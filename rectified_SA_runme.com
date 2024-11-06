#! /bin/tcsh -f
#
#  do SA, but then clean up by taking any atoms that got worse and move them back     -James Holton 10-2-24
#
#
set seed = 1
set temperature = 300
set local_flag = ""
set temperatures = ""

set noise = 3.5
set mapnoise = 3.5

set pdbfile = starthere.pdb
set mtzfile = refme.mtz

set parallel = 0
set debug = 0

set tempfile = /dev/shm/${USER}/temp_SA_$$_
mkdir -p /dev/shm/${USER}

echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set key = `echo $Key | awk '{print tolower($0)}'`
    set val = `echo $Val | awk '{print tolower($0)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`

    if( $assign ) then
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "temp") set temperature = "$Val"
      if("$key" == "parallel") set parallel = "$Val"
    else
      # no equal sign
      if("$key" =~ *.pdb ) set pdbfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = 1
    if("$arg" == "parallel") set parallel = 1
end

if( $debug ) then
    set tempfile = tempfile
endif
set t = "$tempfile"
#touch $logfile

if( ! -e "$pdbfile" ) then
  set BAD = "pdbfile $pdbfile does not exist"
  goto exit
endif

if( $parallel ) goto parallel


set temperature = `echo $temperature | awk '{print $1+0}'`

set pwd = `pwd`
set local_tempdir = /dev/shm/${USER}/temp_SA_$$_

if( $debug < 5 ) then
    mkdir -p $local_tempdir
    cp opts.eff  $local_tempdir >& /dev/null
    cp starthere_*.* $local_tempdir >& /dev/null
    cp *.cif $local_tempdir >& /dev/null
    cp $pdbfile ${local_tempdir}/starthere.pdb
    cp $mtzfile ${local_tempdir}/refme.mtz
    cd $local_tempdir
    echo "tempdir = $local_tempdir"
endif


if(! -e starthere_fullgeo.txt) then

   echo "molprobifying starthere"
   molprobify_runme.com starthere.pdb keepgeo >! molprobify_starthere.log

endif

if(! -e opts.eff) then
  cat << EOF >! opts.eff
refinement {
  refine {
    strategy = *individual_sites *individual_adp
  }
  pdb_interpretation {
    flip_symmetric_amino_acids = False
    automatic_linking {
      amino_acid_bond_cutoff = 2.5
      inter_residue_bond_cutoff = 3.2
    }
  }
}
EOF
endif 

# take any 3-letter-code cif files
set ciffiles = `ls -1rt ???.cif |& grep cif`

echo "runing simulated annealing in phenix"
phenix.refine opts.eff refme.mtz starthere.pdb $ciffiles \
 simulated_annealing.random_seed=$seed \
 prefix=simanneal serial=1 \
 simulated_annealing=True \
 simulated_annealing_torsion=false \
 simulated_annealing.start_temperature=$temperature >! SA.log

if( $status || ! -e simanneal_001.pdb) then
  set BAD = "simulated annealing failed"
  cat SA.log
  goto exit
endif

cp simanneal_001.pdb new.pdb

foreach itr ( 1 2 3 4 5 )

echo "refine $itr"
phenix.refine refme.mtz new.pdb adp.set_b_iso=20 opts.eff $ciffiles \
  wxc_scale=1 \
  main.number_of_macro_cycles=3 \
  prefix=phenix serial=1 >! phenix_${itr}_1.log

echo "refine wxc_scale=0.1"
phenix.refine refme.mtz phenix_001.pdb opts.eff $ciffiles \
  wxc_scale=0.1 \
  main.number_of_macro_cycles=5 \
  prefix=phenix serial=2 >! phenix_${itr}_2.log

echo "refine wxc_scale=1"
phenix.refine refme.mtz phenix_002.pdb opts.eff $ciffiles \
  wxc_scale=1 \
  main.number_of_macro_cycles=5 \
  prefix=phenix serial=3 >! phenix_${itr}_3.log

set num = 003

if(! -e phenix_${num}.pdb) then
  set BAD = "phenix refine failed: no phenix_${num}.pdb"
  goto exit
endif

# look for disimprovements
echo "molprobify..."
molprobify_runme.com phenix_${num}.pdb keepgeo >! molprobify_${itr}.log

awk '{print $0"| OLD"}' starthere_fullgeo.txt |\
cat - phenix_${num}_fullgeo.txt |\
awk '{split($0,w,"|");ent=w[2]}\
 $NF=="OLD"{oldv[ent,$1]=$2;next}\
 {print $1,oldv[ent,$1]-$2,oldv[ent,$1],$2,"|"ent}' |\
sort -k2g |\
awk -v noise=$noise '$2<0 && $4>noise' |\
tee badgeo.txt |\
awk -F "|" '{print $NF}' |\
awk -F "-" '{for(i=1;i<=NF;++i)print substr($i,1,20)}' |\
awk '! seen[$0]{print $0 "| REVERT"} {++seen[$0]}' |\
cat >! revert_geo.txt
# format: atom conf_ typ chain_ resnum | REVERT

egrep "^NONBOND" badgeo.txt |\
awk -F "|" '{print $NF}' |\
awk -F "-" '{for(i=1;i<=NF;++i)print substr($i,1,20)}' |\
awk '{print substr($0,14,1) substr($0,16,5)}' |\
awk '! seen[$0]{print $0 "| REVERTRES"} {++seen[$0]}' |\
cat >! revert_nonbond.txt
# format: chainresnum | REVERTRES



# look for disimproved density
echo "probing density..."
phenix.map_value_at_point phenix_${num}.mtz phenix_${num}.pdb \
          miller_array.labels.name=FOFCWT scale=sigma |\
cat - phenix_${num}.pdb |\
awk '/Map value:/{xyz=$3;rho[xyz]=$NF;n[xyz]=++j;next}\
  ! /^ATOM|^HETAT/{next}\
  {++i;x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
   xyz=sprintf("(%.3f,%.3f,%.3f)",x,y,z);\
  id=substr($0,13,15);\
    a=substr(id,1,4);\
    f=substr(id,5,1);\
    t=substr(id,6,4);\
    c=substr(id,10,1);\
    r=substr(id,11,5);\
    if(f==" ")f="_";if(c==" ")c="_";\
    print "",a,f,t,c,r "|",rho[xyz],"NEW"}' >! new_drho.txt


# look for empties
phenix.map_value_at_point phenix_${num}.mtz starthere.pdb \
          miller_array.labels.name=FOFCWT scale=sigma |\
cat - starthere.pdb |\
awk '/Map value:/{xyz=$3;rho[xyz]=$NF;n[xyz]=++j;next}\
  ! /^ATOM|^HETAT/{next}\
  {++i;x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
   xyz=sprintf("(%.3f,%.3f,%.3f)",x,y,z);\
  id=substr($0,13,15);\
    a=substr(id,1,4);\
    f=substr(id,5,1);\
    t=substr(id,6,4);\
    c=substr(id,10,1);\
    r=substr(id,11,5);\
    if(f==" ")f="_";if(c==" ")c="_";\
    print "",a,f,t,c,r "|",rho[xyz],"OLD"}' >! old_drho.txt


# old not red, new is red, and there is difference > 1 sigma
cat new_drho.txt old_drho.txt |\
 awk '{split($0,w,"|");ent=w[1];rho=$(NF-1)}\
  $NF=="NEW"{newv[ent]=rho+0;next}\
 {print ent "|",rho,newv[ent],rho-newv[ent],"REVERT"}' |\
 sort -k9gr |\
tee sorted_drho.txt |\
awk -v noise=$mapnoise  'BEGIN{red=-sqrt(noise*noise)}\
  $7>red && $8<red && $9>1' |\
cat >! revert_rho.txt
# format: atom conf_ typ chain_ resnum | oldrho newrho deltarho REVERT

cat phenix_${num}_clashes.txt |\
awk -F "|" '{print $2;print $3}' |\
awk '{print substr($0,2,6),"| REVERTRES"}' |\
cat >! revert_clashes.txt
# format: chainresnum | REVERTRES


cat phenix_${num}_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/{energy=($6/0.05)**2;\
   f=$2;c=$4;r=$5;\
   if(f==" ")f="_";if(c==" ")c="_";\
   print energy,f,c,r,"NEW"}' |\
cat >! cbetadevs.txt
# format: E conf resnum NEW

cat starthere_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/{energy=($6/0.05)**2;\
   f=$2;c=$4;r=$5;\
   if(f==" ")f="_";if(c==" ")c="_";\
   id=f" "c" "r;\
   print energy,id,"OLD"}' |\
cat >> cbetadevs.txt
sort -g cbetadevs.txt |\
awk '{f=$2;c=$3;r=$4;\
   if(f==" ")f="_";if(c==" ")c="_";\
  id=f" "c" "r;} \
  ! seen[id] && /OLD/{printf("%s %s%4d | REVERTCONF\n",f,c,r)}\
 {++seen[id]}' |\
cat >! revert_cbetadev.txt
# format: conf chainresnum | REVERTCONF


if(-e always_revert.txt) then
  cp always_revert.txt revert_always.txt
endif

foreach revert ( revert_*.txt )
  echo -n "testing: "
wc -l $revert
cat $revert phenix_${num}.pdb |\
awk '$NF=="REVERT"{++bad[substr($0,1,20)];next}\
  $NF=="REVERTRES"{++bad[substr($0,1,6)];next}\
  $NF=="REVERTCONF"{++bad[substr($0,1,8)];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,13,15);\
  a=substr(id,1,4);\
  f=substr(id,5,1);\
  t=substr(id,6,4);\
  c=substr(id,10,1);\
  r=substr(id,11,5);\
  if(f==" ")f="_";if(c==" ")c="_";\
  eid=" "a" "f" "t" "c" "r;\
  cr=c r;fcr=f" "c r;\
    #print fcr;for(x in bad)if(bad[x])print x;\
  }\
  bad[eid]{print}\
  bad[fcr]{print}\
  bad[cr]{print}\
  ' >! bad.pdb
  wc -l bad.pdb

end


awk -F "|" '{print $1 "|",$NF}' revert_*.txt |\
sort -u >! baddies.txt

set badatoms = `cat baddies.txt | wc -l`
echo "reverting $badatoms bad entities"

cat revert_*.txt phenix_${num}.pdb |\
awk '$NF=="REVERT"{++bad[substr($0,1,20)];next}\
  $NF=="REVERTRES"{++bad[substr($0,1,6)];next}\
  $NF=="REVERCONF"{++bad[substr($0,1,7)];next}\
  ! /^ATOM|^HETAT/{next}\
  {id=substr($0,13,15);\
  a=substr(id,1,4);\
  f=substr(id,5,1);\
  t=substr(id,6,4);\
  c=substr(id,10,1);\
  r=substr(id,11,5);\
  if(f==" ")f="_";if(c==" ")c="_";\
  eid=" "a" "f" "t" "c" "r;\
  cr=c r;fcr=f" "c r}\
  bad[eid]{next}\
  bad[fcr]{next}\
  bad[cr]{next}\
  {print}' >! culled.pdb

wc -l culled.pdb

echo "splicing culled.pdb with starthere.pdb"
combine_pdbs_runme.com culled.pdb starthere.pdb printref=1 outfile=new.pdb >! combine.log

end

echo "finishing..."
mv phenix_${num}.pdb phenix_${temperature}K_${seed}.pdb 
mv phenix_${num}.mtz phenix_${temperature}K_${seed}.mtz 

awk 'substr($0,77,2)!=" H"' starthere.pdb phenix_${temperature}K_${seed}.pdb | rmsd

molprobify_runme.com phenix_${temperature}K_${seed}.pdb 

mkdir -p ${pwd}/pdbs
cp phenix_${temperature}K_${seed}.pdb ${pwd}/pdbs

cd $pwd
if(! $debug ) rm -rf ${local_tempdir}

exit:
if( $?BAD ) then
  echo "ERROR: $BAD"
  exit 9
endif

exit

##########################################################

parallel:

set seed = 0
set round = 0
mkdir -p logs pdbs

if(-e starthere_1.pdb) then
  set round = `ls -1rt starthere_*.pdb | tail -n 1 | awk -F "_" '{print $2+0}'`
endif

if(! -e ./rectified_SA_runme.com) then
  cp $0 ./rectified_SA_runme.com
endif

molprobify_runme.com starthere.pdb keepgeo |& tee molorobify_starthere.log

set wE0 = `awk '/weighted energy/{print $NF;exit}' molorobify_starthere.log`

set since_better = 0
while ( $since_better < 10 )

@ round = ( $round + 1 )

cp starthere.pdb starthere_${round}.pdb
echo "clearing logs"
find logs -name 'trial_*.log' -exec rm -f \{\} \;
echo "clearing old pdbs"
find pdbs -name 'phenix_*.pdb' -exec rm -f \{\} \;

#set temperatures = `awk 'BEGIN{for(T=2000;T<8000;T*=1.1)print int(T)}' | sort -u | sort -g`
if ( $#temperatures <= 1 ) then
  set temperatures = `awk 'BEGIN{for(T=1000;T<=6000;T+=1000)print int(T)}' | sort -u | sort -g`
endif

foreach temperature ( $temperatures )

foreach s ( `seq 1 10` )
@ seed = ( $seed + 1 )

srun ./rectified_SA_runme.com seed=$seed temperature=$temperature noise=$noise mapnoise=$mapnoise debug=$debug >! logs/trial_${temperature}K_${seed}.log &

sleep 0.5

end
end
wait

#append_file_date.com sorted.txt
tail -n 1 logs/trial_* | awk 'NF>10{split($(NF-2),w,"_");print w[2],w[3],$1,$2,$3,$4,$12}' |\
 sort -k3gr | tee sorted.txt
cp sorted.txt sorted_${round}.txt 

set temperature = `tail -n 1 sorted.txt | awk '{print $1}'`
set seed = `tail -n 1 sorted.txt | awk '{print $2}'`

set bestpdb = pdbs/phenix_${temperature}_${seed}.pdb
set bestlog = logs/trial_${temperature}_${seed}.log
cp $bestpdb best_for_round${round}.pdb
cp $bestlog best_for_round${round}.log

echo "tail of best log:"
tail -100 $bestlog

set wE = `awk '/weighted energy/{print $NF;exit}' $bestlog`


set better = `echo $wE $wE0 | awk '{print ( $1 < $2 )}'`
if( $better ) then
  echo "saving new best"
  cp $bestpdb best_so_far.pdb
  cp $bestlog best_so_far.log
  cp best_so_far.pdb starthere.pdb
  molprobify_runme.com starthere.pdb keepgeo |& tee molorobify_starthere.log

  set wE0 = `awk '/weighted energy/{print $NF;exit}' molorobify_starthere.log`
  set since_better = 0
else
  @ since_better = ( $since_better + 1 )
  grep weighted molorobify_starthere.log
endif

echo "$since_better rounds since wE got better"

if(-e exit) exit

end



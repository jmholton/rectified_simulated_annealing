#! /bin/tcsh -f
#
#  do SA, but then clean up by taking any atoms that got worse and move them back     -James Holton 10-2-24
#
#
set seed = 1
set temperature = 300
set local_flag = ""

set noise = 3.0

set pdbfile = starthere.pdb
set mtzfile = refme.mtz

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
    else
      # no equal sign
      if("$key" =~ *.pdb ) set pdbfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = 1
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




set temperature = `echo $temperature | awk '{print $1+0}'`

set pwd = `pwd`
set local_tempdir = /dev/shm/${USER}/temp_SA_$$_

if(! $debug ) then
    mkdir -p $local_tempdir
    cp opts.eff  $local_tempdir >& /dev/null
    cp starthere_*.* $local_tempdir >& /dev/null
    cp *.cif $local_tempdir >& /dev/null
    cp $pdbfile ${local_tempdir}/starthere.pdb
    cp $mtzfile ${local_tempdir}/refme.mtz
    cd $local_tempdir
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

echo "refine"
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

egrep "^NONBOND" badgeo.txt |\
awk -F "|" '{print $NF}' |\
awk -F "-" '{for(i=1;i<=NF;++i)print substr($i,1,20)}' |\
awk '{print substr($0,14,1) substr($0,16,5)}' |\
awk '! seen[$0]{print $0 "| REVERTRES"} {++seen[$0]}' |\
cat >! revert_nonbond.txt



# look for disimproved density
echo "probing density..."
phenix.map_value_at_point phenix_${num}.mtz phenix_${num}.pdb \
          miller_array.labels.name=FOFCWT scale=sigma |\
cat - phenix_${num}.pdb |\
awk '/Map value:/{rho[$3]=$NF;n[$3]=++j;next}\
  /^CRYST1/{print}\
  ! /^ATOM|^HETAT/{next}\
  {++i;x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
   xyz=sprintf("(%.3f,%.3f,%.3f)",x,y,z);\
  id=substr($0,14,15);\
    a=substr(id,1,4);\
    f=substr(id,5,1);\
    t=substr(id,6,4);\
    c=substr(id,10,1);\
    r=substr(id,11,5);\
    print "",a,f,t,c,r "|",rho[xyz],"NEW"}' >! new_drho.txt


# look for empties
phenix.map_value_at_point phenix_${num}.mtz starthere.pdb \
          miller_array.labels.name=FOFCWT scale=sigma |\
cat - starthere.pdb |\
awk '/Map value:/{rho[$3]=$NF;n[$3]=++j;next}\
  /^CRYST1/{print}\
  ! /^ATOM|^HETAT/{next}\
  {++i;x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
   xyz=sprintf("(%.3f,%.3f,%.3f)",x,y,z);\
  id=substr($0,14,15);\
    a=substr(id,1,4);\
    f=substr(id,5,1);\
    t=substr(id,6,4);\
    c=substr(id,10,1);\
    r=substr(id,11,5);\
    print "",a,f,t,c,r "|",rho[xyz],"OLD"}' >! old_drho.txt


cat new_drho.txt old_drho.txt |\
 awk '{split($0,w,"|");ent=w[1];rho=$(NF-1)}\
  $NF=="NEW"{newv[ent]=rho+0;next}\
 {print ent "|",rho,newv[ent],rho-newv[ent],"REVERT"}' |\
 sort -k9gr |\
awk '$7>-3.5 && $8<-3.5 && $9>1' |\
cat >! revert_rho.txt
# old not red, new is red, and there is difference

cat phenix_${num}_clashes.txt |\
awk -F "|" '{print $2;print $3}' |\
awk '{print substr($0,2,6),"| REVERTRES"}' |\
cat >! revert_clashes.txt



cat phenix_${num}_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/{energy=($6/0.05)**2;\
   print energy,$4,$2,$5,"NEW"}' |\
cat >! cbetadevs.txt

cat starthere_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/{energy=($6/0.05)**2;\
   print energy,$4,$2,$5,"OLD"}' |\
cat >> cbetadevs.txt
sort -g cbetadevs.txt |\
awk '{id=$2" "$3" "$4} \
  ! seen[id] && /OLD/{printf("%s %s%4d | REVERTCONF\n",$3,$2,$4)}\
 {++seen[id]}' |\
cat >! revert_cbetadev.txt




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




set temperatures = `awk 'BEGIN{for(T=2000;T<8000;T*=1.1)print int(T)}' | sort -u | sort -g`
foreach temperature ( $temperatures )
#set temperature = 3000
#set temperature = 2500

foreach seed ( `seq 1 20` )

srun ../SA_runme.com $seed $temperature >! logs/trial_${temperature}K_${seed}.log &

sleep 0.5

end
end






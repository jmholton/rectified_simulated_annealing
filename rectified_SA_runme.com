#! /bin/tcsh -f
#
#  do SA, but then clean up by taking any atoms that got worse and move them back     -James Holton 4-1-25
#
#
set seed = 1
set temperature = 300
set final_wxc = 0.5
set local_flag = ""
set temperatures = ""
set nseeds = 10

set noise = 3.5
set mapnoise = 3.5
# C-beta deviations are 3D distances: per-component sigma (A) of the chi(3) deviate, calibrated on correct models (wE v2)
set cbsigma = 0.032

# memory of good models: every accepted model is kept in memory/ (ledger.txt: round wE Rwork Rfree file)
# and the winner of each round is checked against all of them for atoms that crept out of density
set memory = 1
# a trial is eligible only if its Rfree (%) is within this much of the lowest Rfree in memory (negative: no limit)
set rfree_slack = 0.2
# hold every water near its position in the first remembered model (the start of the whole run): before each refinement the
# waters are put back there and restrained by a harmonic restraint of this sigma (A); 0 = off.  wE has no stake in the waters,
# so nothing else holds them, and above ~5000 K the simulated annealing scatters them (1aho: 100 of 154 move > 1 A at 6000 K)
set water_sigma = 0.1

set pdbfile = starthere.pdb
set mtzfile = refme.mtz
set ciffiles = "auto"
set efffiles = opts.eff

set parallel = 0
set debug = 0
set srun = srun

set tempfile = ""

# the helper scripts that ship with this one come first in the path, and the trials (which run as ./rectified_SA_runme.com) inherit it
if(! $?RSA_HOME) then
  set rsa_home = `dirname $0`
  setenv RSA_HOME `cd $rsa_home ; pwd`
endif
set path = ( $RSA_HOME $path )

echo "command-line arguments: $* "

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
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

# tcsh reads its script from disk as it runs, so editing this file (or switching git branches) under a running parent
# changes what it does mid-round.  The parent therefore runs from a private copy in the working directory.
if( $parallel && ! $?RSA_PARENT_COPY ) then
  setenv RSA_PARENT_COPY ./rectified_SA_parent.com
  cp $0 $RSA_PARENT_COPY
  echo "running the parent from a private copy: $RSA_PARENT_COPY"
  exec tcsh -f $RSA_PARENT_COPY $*
endif

set dependencies = ( molprobify_runme.com combine_pdbs_runme.com reset_waters_runme.com )
if( $parallel && $memory ) set dependencies = ( $dependencies memory_repair_runme.com )
foreach dependency ( $dependencies )
  echo -n "using: "
  which $dependency 
  if( $status ) then
      set BAD = "cannot run without $dependency "
      goto exit
  endif
end

if( "$ciffiles" == "auto" ) then
   # take any 3-letter-code cif files
   set ciffiles = `ls -1rt ???.cif |& grep cif`
   echo "selected cif files: $ciffiles"
endif


if ( "$tempfile" == "" ) set tempfile = /dev/shm/${USER}/temp_SA_${seed}_$$_
mkdir -p `dirname $tempfile `

if( $debug ) then
    set tempfile = tempfile
endif
set t = "$tempfile"
#touch $logfile

set temperatures = `echo $temperatures | awk '{gsub("[,:K]"," ");print}'`

if( ! -e "$pdbfile" ) then
  set BAD = "pdbfile $pdbfile does not exist"
  goto exit
endif

set oride
if(-e overrides.txt) set oride = overridefile=overrides.txt
echo $oride

if( $parallel ) goto parallel


set temperature = `echo $temperature | awk '{print $1+0}'`

set pwd = `pwd`
set local_tempdir = /dev/shm/${USER}/temp_SA_$$_

if( $debug < 5 ) then
    mkdir -p $local_tempdir
    cp opts.eff  $local_tempdir >& /dev/null
    cp starthere_*.* $local_tempdir >& /dev/null
    cp *.cif $local_tempdir >& /dev/null
    cp overrides.txt $local_tempdir >& /dev/null
    cp water_reference.pdb $local_tempdir >& /dev/null
    cp $pdbfile ${local_tempdir}/starthere.pdb
    cp $mtzfile ${local_tempdir}/refme.mtz
    cd $local_tempdir
    echo "local_tempdir = $local_tempdir"
endif


if(! -e starthere_fullgeo.txt) then

   echo "molprobifying starthere"
   molprobify_runme.com starthere.pdb keepgeo=1 debug=1 $oride $ciffiles >! molprobify_starthere.log
   if( $status ) then
     set BAD = "cannot score structrues without molprobify_runme.com"
     goto exit
   endif

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

if( "$ciffiles" == "auto" ) then
   # take any 3-letter-code cif files
   set ciffiles = `ls -1rt ???.cif |& grep cif`
   echo "selected cif files: $ciffiles"
endif

# hold the waters: put them back at water_reference.pdb before every refinement and restrain them there
set holdwat = ""
set hold = `echo $water_sigma | awk '{print ($1+0>0)}'`
if( $hold && ! -e water_reference.pdb ) then
  # a trial run on its own (the parallel parent supplies the reference of the whole run) holds the waters of its start
  awk '/^ATOM|^HETATM/ && substr($0,18,3)=="HOH"' starthere.pdb >! water_reference.pdb
  echo "no water_reference.pdb: holding the waters of starthere.pdb"
endif
if( $hold ) then
  cat << EOF >! hold_waters.eff
refinement.pdb_interpretation.reference_coordinate_restraints {
  enabled = True
  exclude_outliers = False
  selection = "water"
  sigma = $water_sigma
  limit = 10
}
EOF
  set holdwat = hold_waters.eff
  reset_waters_runme.com water_reference.pdb starthere.pdb
endif

echo "runing simulated annealing in phenix"
phenix.refine opts.eff $holdwat refme.mtz starthere.pdb $ciffiles \
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
if( $hold ) reset_waters_runme.com water_reference.pdb new.pdb
phenix.refine refme.mtz new.pdb adp.set_b_iso=20 opts.eff $holdwat $ciffiles \
  wxc_scale=1 \
  main.number_of_macro_cycles=3 \
  prefix=phenix serial=1 >! phenix_${itr}_1.log

echo "refine wxc_scale=10"
phenix.refine refme.mtz phenix_001.pdb opts.eff $holdwat $ciffiles \
  wxc_scale=10 \
  main.number_of_macro_cycles=5 \
  prefix=phenix serial=2 >! phenix_${itr}_2.log

echo "refine wxc_scale=0.1"
phenix.refine refme.mtz phenix_002.pdb opts.eff $holdwat $ciffiles \
  wxc_scale=0.1 \
  main.number_of_macro_cycles=5 \
  prefix=phenix serial=3 >! phenix_${itr}_3.log

echo "refine wxc_scale=${final_wxc}"
phenix.refine refme.mtz phenix_003.pdb opts.eff $holdwat $ciffiles \
  wxc_scale=${final_wxc} \
  main.number_of_macro_cycles=5 \
  prefix=phenix serial=4 >! phenix_${itr}_4.log

set num = 004

if(! -e phenix_${num}.pdb) then
  set BAD = "phenix refine failed: no phenix_${num}.pdb"
  goto exit
endif

# look for disimprovements
echo "molprobify..."
set oride
if(-e overrides.txt) set oride = overridefile=overrides.txt
molprobify_runme.com phenix_${num}.pdb keepgeo=1 $oride $ciffiles >! molprobify_${itr}.log
grep wE molprobify_${itr}.log

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


# C-beta deviation: revert the residue conformer only if it got worse AND is beyond noise.  The deviation is a 3D distance,
# so x = dev/cbsigma is a chi(3) deviate; judge it by the Gaussian deviate with the same two-sided tail probability (wE v2):
#   P(>x) = erfc(x/sqrt(2)) + sqrt(2/pi) x exp(-x^2/2) < erfc(noise/sqrt(2)),  i.e. dev > 0.135 A for noise=3.5, cbsigma=0.032
cat phenix_${num}_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/ && NF>5{f=$2;c=$4;r=$5;\
   if(f==" ")f="_";if(c==" ")c="_";\
   print $6+0,f,c,r,"NEW"}' |\
cat >! cbetadevs.txt
cat starthere_cbetadev.log |\
 awk -F ":" '! /^SUMM|^pdb/ && NF>5{f=$2;c=$4;r=$5;\
   if(f==" ")f="_";if(c==" ")c="_";\
   print $6+0,f,c,r,"OLD"}' |\
cat >> cbetadevs.txt
# format: dev conf chain resnum NEW/OLD
cat cbetadevs.txt |\
awk -v noise=$noise -v cbsigma=$cbsigma '\
 function erfc(x,  z,t,a){z=(x<0?-x:x);t=1/(1+0.5*z);\
   a=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));\
   return (x>=0?a:2-a)}\
 function pchi3(x){return erfc(x/sqrt(2))+sqrt(2/3.14159265358979)*x*exp(-x*x/2)}\
 BEGIN{pnoise=erfc(noise/sqrt(2))}\
 {id=$2" "$3" "$4}\
 $NF=="OLD"{old[id]=$1;next}\
 $NF=="NEW"{new[id]=$1}\
 END{for(id in new) if((id in old) && new[id]>old[id] && pchi3(new[id]/cbsigma)<pnoise){\
   split(id,w," ");printf("%s %s%4d | REVERTCONF\n",w[1],w[2],w[3])}}' |\
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
  cr=c r;fcr=f" "c r}\
  bad[eid]{next}\
  bad[fcr]{next}\
  bad[cr]{next}\
  {print}' >! culled.pdb

wc -l culled.pdb

echo "splicing culled.pdb with starthere.pdb"
combine_pdbs_runme.com culled.pdb starthere.pdb printref=1 outfile=new.pdb >! combine.log

end

echo "final refine run..."
if( $hold ) reset_waters_runme.com water_reference.pdb new.pdb
phenix.refine refme.mtz new.pdb opts.eff $holdwat $ciffiles \
  main.number_of_macro_cycles=10 \
  prefix=phenix serial=1 >! phenix_final.log

echo "finishing..."
mv phenix_001.pdb phenix_${temperature}K_${seed}.pdb 
mv phenix_001.log phenix_${temperature}K_${seed}.log 
mv phenix_001.mtz phenix_${temperature}K_${seed}.mtz 

awk 'substr($0,77,2)!=" H"' starthere.pdb phenix_${temperature}K_${seed}.pdb | rmsd

molprobify_runme.com phenix_${temperature}K_${seed}.pdb $oride $ciffiles

mkdir -p ${pwd}/pdbs
foreach ext ( pdb mtz log )
cp phenix_${temperature}K_${seed}.${ext} ${pwd}/pdbs
end

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

molprobify_runme.com starthere.pdb keepgeo=1 $oride $ciffiles |& tee molorobify_starthere.log
if( $status ) then
  set BAD = "cannot score structrues without molprobify_runme.com"
  goto exit
endif

set wE0 = `awk '/weighted energy/{print $NF;exit}' molorobify_starthere.log`

# remember the starting model
mkdir -p memory
if(! -e memory/ledger.txt) then
  cp starthere.pdb memory/round${round}.pdb
  tail -n 1 molorobify_starthere.log |\
  awk -v r=$round -v wE=$wE0 '{print r,wE,$3,$4,"memory/round" r ".pdb"}' >! memory/ledger.txt
endif

# the water reference: waters of the first remembered model (the start of the whole run), kept across restarts
set hold = `echo $water_sigma | awk '{print ($1+0>0)}'`
if( $hold && ! -e water_reference.pdb ) then
  set first = `awk 'NR==1{print $5}' memory/ledger.txt`
  if( ! -e "$first" ) set first = starthere.pdb
  awk '/^ATOM|^HETATM/ && substr($0,18,3)=="HOH"' $first >! water_reference.pdb
  echo "holding `grep -c . water_reference.pdb` water atoms near their positions in $first (sigma $water_sigma A)"
endif
if( $hold ) then
  cat << EOF >! hold_waters.eff
refinement.pdb_interpretation.reference_coordinate_restraints {
  enabled = True
  exclude_outliers = False
  selection = "water"
  sigma = $water_sigma
  limit = 10
}
EOF
endif

set since_better = 0
while ( $since_better < 10 )

@ round = ( $round + 1 )

cp starthere.pdb starthere_${round}.pdb
echo "clearing logs"
find logs -name 'trial_*.log' -exec rm -f \{\} \;
echo "clearing old pdbs"
find pdbs -name 'phenix_*.pdb' -exec rm -f \{\} \;

echo "temperatures: $temperatures"
if ( $#temperatures <= 1 ) then
  set temperatures = `awk 'BEGIN{for(T=1000;T<=8000;T+=1000)print int(T)}' | sort -u | sort -g`
  echo "using default temperatures: $temperatures"
endif

foreach temperature ( $temperatures )

foreach s ( `seq 1 $nseeds` )
@ seed = ( $seed + 1 )

$srun ./rectified_SA_runme.com seed=$seed temperature=$temperature \
   noise=$noise mapnoise=$mapnoise cbsigma=$cbsigma water_sigma=$water_sigma final_wxc=$final_wxc debug=$debug >! logs/trial_${temperature}K_${seed}.log &

sleep 0.5

end
end
wait

#append_file_date.com sorted.txt
tail -n 1 logs/trial_* | awk 'NF>10{split($(NF-2),w,"_");print w[2],w[3],$1,$2,$3,$4,$12}' |\
 sort -k4gr | tee sorted.txt
cp sorted.txt sorted_${round}.txt 

# only trials that fit the data about as well as the best remembered model are eligible
set rmin = `awk '$4+0>0{print $4}' memory/ledger.txt | sort -g | head -n 1`
if( "$rmin" == "" ) set rmin = `awk '$6+0>0{print $6}' sorted.txt | sort -g | head -n 1`
awk -v rmin=$rmin -v slack=$rfree_slack 'slack<0 || ( $6+0>0 && $6+0<=rmin+slack )' sorted.txt >! eligible.txt
cp eligible.txt eligible_${round}.txt
echo "lowest remembered Rfree: $rmin %  eligible trials (Rfree <= $rmin + $rfree_slack): `cat eligible.txt | wc -l` of `cat sorted.txt | wc -l`"

if( -z eligible.txt ) then
  echo "no trial fits the data within $rfree_slack % of the best remembered model"
  @ since_better = ( $since_better + 1 )
else

set temperature = `tail -n 1 eligible.txt | awk '{print $1}'`
set seed = `tail -n 1 eligible.txt | awk '{print $2}'`

set bestphenix = pdbs/phenix_${temperature}_${seed}
set besttrial = logs/trial_${temperature}_${seed}
cp ${bestphenix}.pdb best_for_round${round}.pdb
cp ${bestphenix}.log best_for_round${round}.log
cp ${besttrial}.log best_for_round${round}_trial.log

echo "tail of best log:"
tail -100 ${besttrial}.log

#set wE0 = `awk '/weighted energy/{print $NF}' molorobify_starthere.log | tail -n 1`
set wE = `awk '/weighted energy/{print $NF}' ${besttrial}.log | tail -n 1`
set candidate = ${bestphenix}
set candR = `tail -n 1 eligible.txt | awk '{print $5, $6}'`

if( $memory ) then
  # put back residue conformers that have crept out of density since an earlier, remembered model
  echo "checking the winner against `ls -1 memory/round*.pdb | wc -l` remembered models"
  memory_repair_runme.com ${bestphenix}.pdb ${bestphenix}.mtz memory/round*.pdb \
     mapnoise=$mapnoise outfile=repair_${round}.pdb >&! repair_${round}.log
  grep -E '^restore|^ERROR' repair_${round}.log
  set nrestore = `awk '/^restored residue-conformers:/{print $NF}' repair_${round}.log`
  if( "$nrestore" != "" && "$nrestore" != "0" ) then
    echo "refining the repaired model"
    mkdir -p repairs
    set rnum = `echo $round | awk '{printf("%03d",$1)}'`
    set cifs = ""
    foreach cif ( $ciffiles )
      set cifs = ( $cifs ../$cif )
    end
    set roride = ""
    if(-e overrides.txt) set roride = overridefile=../overrides.txt
    set rhold = ""
    if( $hold ) then
      reset_waters_runme.com water_reference.pdb repair_${round}.pdb
      set rhold = ../hold_waters.eff
    endif
    ( cd repairs ; \
      $srun phenix.refine ../refme.mtz ../repair_${round}.pdb ../opts.eff $rhold $cifs \
        main.number_of_macro_cycles=10 \
        prefix=repaired serial=$round >&! refine_${rnum}.log ; \
      molprobify_runme.com repaired_${rnum}.pdb $roride $cifs >&! molprobify_${rnum}.log )
    if( -e repairs/repaired_${rnum}.pdb ) then
      set rwE = `awk '/weighted energy/{print $NF}' repairs/molprobify_${rnum}.log | tail -n 1`
      set rR = `tail -n 1 repairs/molprobify_${rnum}.log | awk '{print $3, $4}'`
      set rpass = `tail -n 1 repairs/molprobify_${rnum}.log | awk -v rmin=$rmin -v slack=$rfree_slack '{print ( slack<0 || ( $4+0>0 && $4+0<=rmin+slack ) )}'`
      echo "repaired: wE $rwE  R/Rfree $rR   unrepaired: wE $wE  R/Rfree $candR"
      if( "$rwE" != "" && "$rpass" == "1" ) then
        set wE = $rwE
        set candR = ( $rR )
        set candidate = repairs/repaired_${rnum}
      else
        echo "the repaired model does not fit the data well enough, keeping the winner as it is"
      endif
    else
      echo "refinement of the repaired model failed, keeping the winner as it is"
    endif
  endif
endif

set better = `echo $wE $wE0 | awk '{print ( $1 < $2 )}'`
if( $better ) then
  echo "saving new best"
  cp ${candidate}.pdb best_so_far.pdb
  cp ${candidate}.log best_so_far.log
  cp ${candidate}.mtz best_so_far.mtz
  cp ${besttrial}.log best_trial_so_far.log
  cp best_so_far.pdb starthere.pdb
  if( $memory ) then
    cp best_so_far.pdb memory/round${round}.pdb
    echo "$round $wE $candR memory/round${round}.pdb" >> memory/ledger.txt
  endif
  molprobify_runme.com starthere.pdb keepgeo=1 $oride $ciffiles |& tee molorobify_starthere.log

  set wE0 = `awk '/weighted energy/{print $NF;exit}' molorobify_starthere.log`
  set since_better = 0
else
  echo "best run is not better than starting point"
  @ since_better = ( $since_better + 1 )
  grep weighted molorobify_starthere.log
endif

endif

echo "$since_better rounds since wE got better"

if(-e ./exit) then
   echo "exiting because ./exit exists"
   rm ./exit
   goto exit
endif

end


exit:
if($?BAD) then
   echo "ERROR: $BAD"
   echo "tempdir = $tempdir"
   exit 9
endif

# clean up
if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif


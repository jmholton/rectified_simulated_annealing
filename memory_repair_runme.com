#! /bin/tcsh -f
#
#  put back residue conformers that have crept out of density, taking them from the
#  remembered model whose position the current difference map likes best
#
#  usage: memory_repair_runme.com model.pdb model.mtz memory1.pdb [memory2.pdb ...] [mapnoise=3.5] [outfile=repaired.pdb]
#
#  An atom of model.pdb has "crept" if the mFo-DFc map (FOFCWT of model.mtz, in sigma) is below -mapnoise at it,
#  while at the same atom of a remembered model it is above -mapnoise and at least 1 sigma higher.  This is the
#  rectified-SA revert test, but against every remembered model instead of only the start of the current round,
#  so an atom that slid out of density a little at a time is still caught.  The whole residue conformer is copied
#  back from the remembered model the map likes best at its most negative crept atom.  Hydrogens do not trigger.
#
set mapnoise = 3.5
set outfile = repaired.pdb
set pdbfile = ""
set mtzfile = ""
set memories = ""
set debug = 0
set tempfile = /dev/shm/${USER}/temp_memrep_$$_
set nonomatch

foreach Arg ( $* )
    set Key = `echo $Arg | awk -F "=" '{print tolower($1)}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    if( "$Arg" =~ *=* ) then
        if( "$Key" == "mapnoise" ) set mapnoise = "$Val"
        if( "$Key" == "outfile" ) set outfile = "$Val"
        if( "$Key" == "debug" ) set debug = "$Val"
        if( "$Key" == "tempfile" ) set tempfile = "$Val"
        continue
    endif
    if( "$Arg" =~ *.mtz ) set mtzfile = "$Arg"
    if( "$Arg" =~ *.pdb ) then
        if( "$pdbfile" == "" ) then
            set pdbfile = "$Arg"
        else if( "$Arg" != "$pdbfile" ) then
            set memories = ( $memories "$Arg" )
        endif
    endif
end

if( ! -e "$pdbfile" || ! -e "$mtzfile" || "$memories" == "" ) then
    set BAD = "usage: $0 model.pdb model.mtz memory1.pdb [memory2.pdb ...] [mapnoise=$mapnoise] [outfile=$outfile]"
    goto exit
endif

if( $debug ) set tempfile = tempfile_memrep_
set t = "$tempfile"
mkdir -p `dirname $t`
rm -f ${t}rho_*.txt

# FOFCWT at every atom of the model and of each remembered model, all in the model's own map
set n = 0
foreach pdb ( $pdbfile $memories )
    @ n = ( $n + 1 )
    ( phenix.map_value_at_point $mtzfile $pdb miller_array.labels.name=FOFCWT scale=sigma |\
    cat - $pdb |\
    awk -v src=$pdb '/Map value:/{xyz=$3;rho[xyz]=$NF;next}\
      ! /^ATOM|^HETAT/{next}\
      {x=substr($0,31,8);y=substr($0,39,8);z=substr($0,47,8);\
       xyz=sprintf("(%.3f,%.3f,%.3f)",x,y,z);\
       if(! (xyz in rho)) next;\
       id=substr($0,13,15);\
       a=substr(id,1,4);\
       f=substr(id,5,1);\
       t=substr(id,6,4);\
       c=substr(id,10,1);\
       r=substr(id,11,5);\
       if(f==" ")f="_";if(c==" ")c="_";\
       el=substr($0,77,2);gsub(" ","",el);if(el=="")el=substr(a,2,1);\
       print "",a,f,t,c,r "|",rho[xyz],el,src}' >! ${t}rho_${n}.txt ) &
end
wait

# crept atoms, the remembered model the map likes best for each, and one source per residue conformer
cat ${t}rho_*.txt |\
awk -F "|" -v model=$pdbfile -v noise=$mapnoise 'BEGIN{red=-sqrt(noise*noise)}\
  {ent=$1;split($2,v," ");rho=v[1]+0;el=v[2];src=v[3]}\
  src==model{new[ent]=rho;elem[ent]=el;next}\
  ! (ent in best) || rho>best[ent]{best[ent]=rho;from[ent]=src}\
  END{for(ent in new){\
     if(elem[ent]=="H" || ! (ent in best)) continue;\
     if(new[ent]<red && best[ent]>red && best[ent]-new[ent]>1) print new[ent],best[ent],from[ent],"|" ent}}' |\
sort -g |\
awk -F "|" '{ent=$2;key=substr(ent,7,1) substr(ent,14,1) substr(ent,16,5);split($1,v," ")}\
  ! seen[key]++{print key "|",v[3],v[1],v[2],"|" ent}' |\
cat >! ${t}restore.txt
# format: conf chain resnum | source newrho sourcerho | crept atom

set nrestore = `cat ${t}restore.txt | wc -l`
echo "restored residue-conformers: $nrestore"
awk -F "|" '{split($2,v," ");printf("restore %s from %s: crept atom%s at %s sigma, %s sigma at the remembered position\n",$1,v[1],$3,v[2],v[3])}' ${t}restore.txt

if( $nrestore == 0 ) then
    cp $pdbfile $outfile
    goto exit
endif

# copy those residue conformers from their sources and splice them into the model
rm -f ${t}restore.pdb
foreach src ( `awk -F "|" '{split($2,v," ");print v[1]}' ${t}restore.txt | sort -u` )
    awk -F "|" -v src=$src '{split($2,v," ")} v[1]==src{print $1 "| RESTORE"}' ${t}restore.txt |\
    cat - $src |\
    awk '$NF=="RESTORE"{++want[substr($0,1,7)];next}\
      ! /^ATOM|^HETAT/{next}\
      {f=substr($0,17,1);c=substr($0,22,1);r=substr($0,23,5);\
       if(f==" ")f="_";if(c==" ")c="_";}\
      want[f c r]{print}' >> ${t}restore.pdb
end
combine_pdbs_runme.com ${t}restore.pdb $pdbfile printref=1 outfile=$outfile >! ${t}combine.log
echo "restored atoms: `cat ${t}restore.pdb | wc -l`  written to $outfile"

exit:
if( $?BAD ) then
    echo "ERROR: $BAD"
    exit 9
endif
if( ! $debug ) rm -f ${t}rho_*.txt ${t}restore.txt ${t}restore.pdb ${t}combine.log
exit

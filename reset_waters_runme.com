#! /bin/tcsh -f
#
#  put every water of a model back where it is in a reference model                -James Holton 27-10-26
#  usage: reset_waters_runme.com reference.pdb model.pdb [outfile]   (default: overwrite model.pdb)
#  waters are matched by name, altloc, residue name, chain and number; waters missing from the reference are left alone
#
set ref = "$1"
set pdb = "$2"
set outfile = "$3"
if( "$outfile" == "" ) set outfile = "$pdb"
if( ! -e "$ref" || ! -e "$pdb" ) then
  echo "usage: $0 reference.pdb model.pdb [outfile]"
  exit 9
endif
set tmp = ${outfile}.reset_waters_$$
awk 'FNR==NR{if(/^ATOM|^HETATM/ && substr($0,18,3)=="HOH")xyz[substr($0,13,15)]=substr($0,31,24);next}\
  /^ATOM|^HETATM/ && substr($0,18,3)=="HOH" && (substr($0,13,15) in xyz){\
    $0=substr($0,1,30) xyz[substr($0,13,15)] substr($0,55);++n}\
  {print}\
  END{print "reset_waters: " n+0 " waters put back" > "/dev/stderr"}' "$ref" "$pdb" >! $tmp
mv $tmp "$outfile"

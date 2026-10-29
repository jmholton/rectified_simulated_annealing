#! /bin/tcsh -f
#
#  put every water of a model back where it is in a reference model                -James Holton 27-10-26
#  usage: reset_waters_runme.com reference.pdb model.pdb [outfile] [labels=free|same]   (default outfile: overwrite model.pdb)
#
#  labels=free (default): the reference holds the SITES, not the conformer letters.  A water split into A and B goes back
#  onto its two reference sites, each atom onto the site it is nearest (keeping its occupancy and B), and then takes whichever
#  A/B labelling overlaps least with the rest of the model: same-conformer (or blank) atoms closer than a contact distance
#  (2.2 A to water O, 2.4 to O, 2.55 to N, 3.1 to C, 3.0 to S, 2.62 to non-polar H, 1.6 to polar H; the scorer's nonbond and
#  clash limits; symmetry included) cost (d0-d)^2.
#  So a held water follows its neighbours' letters as they untangle, instead of pinning a clash between conformers.
#  labels=same: waters are matched by name, altloc, residue name, chain and number (the original behaviour).
#  Waters missing from the reference are left alone.
#
set ref = ""
set pdb = ""
set outfile = ""
set labels = free
foreach arg ( $* )
  if( "$arg" =~ *=* ) then
    set key = `echo $arg | awk -F "=" '{print $1}'`
    set val = `echo $arg | awk -F "=" '{print $2}'`
    if( "$key" == "labels" ) set labels = "$val"
    continue
  endif
  if( "$ref" == "" ) then
    set ref = "$arg"
  else if( "$pdb" == "" ) then
    set pdb = "$arg"
  else
    set outfile = "$arg"
  endif
end
if( "$outfile" == "" ) set outfile = "$pdb"
if( ! -e "$ref" || ! -e "$pdb" ) then
  echo "usage: $0 reference.pdb model.pdb [outfile] [labels=free|same]"
  exit 9
endif
set tmp = ${outfile}.reset_waters_$$

if( "$labels" == "free" ) then
  set python = ""
  foreach try ( phenix.python libtbx.python cctbx.python )
    which $try >& /dev/null
    if( ! $status ) then
      set python = $try
      break
    endif
  end
  if( "$python" == "" ) then
    echo "reset_waters: no phenix.python, keeping the reference labels" > /dev/stderr
    set labels = same
  endif
endif

if( "$labels" == "free" ) then
  cat << 'EOF' >! ${tmp}.py
import sys
from iotbx import pdb
from cctbx import crystal
from scitbx.array_family import flex
ref_file, model_file, out_file = sys.argv[1:4]
D0 = {'O': 2.4, 'N': 2.55, 'C': 3.1, 'S': 3.0, 'H': 2.62}   # contact distance of a water O to each element
D0_WATER = 2.2
D0_POLAR_H = 1.6
CUT = 3.2

def is_atom(l): return l[:6] in ('ATOM  ', 'HETATM')
def is_water(l): return is_atom(l) and l[17:20] == 'HOH'
def site(l): return l[12:16] + l[17:27]               # name, residue name, chain, number, insertion code
def xyz(l): return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
def element(l):
    e = l[76:78].strip().upper()
    if not e:
        e = l[12:14].strip().upper().lstrip('0123456789') or l[13].upper()
    return e[:1] if e[:1] in 'HCNOS' else e
def d2(a, b): return sum((p - q) ** 2 for p, q in zip(a, b))

ref = {}
for l in open(ref_file):
    if is_water(l): ref.setdefault(site(l), {})[l[16]] = xyz(l)
lines = open(model_file).read().splitlines()
pi = pdb.input(file_name=model_file)
cs = pi.crystal_symmetry()
if cs is None or cs.unit_cell() is None:
    cs = pdb.input(file_name=ref_file).crystal_symmetry()

# the model's water records, by site and altloc
recs = {}
for i, l in enumerate(lines):
    if is_water(l) and site(l) in ref: recs.setdefault(site(l), {})[l[16]] = i
# split waters: two altloc records in the model and the same two in the reference
split = [k for k in recs if len(recs[k]) == 2 and set(recs[k]) == set(ref[k])]
out = list(lines)
# every other water goes back by its own label
nback = 0
for k in recs:
    if k in split: continue
    for al, i in recs[k].items():
        if al in ref[k]:
            x, y, z = ref[k][al]
            out[i] = lines[i][:30] + '%8.3f%8.3f%8.3f' % (x, y, z) + lines[i][54:]
            nback += 1
# split waters: each model record onto its nearest reference site
w_sites = []                  # per split water: ((letter, xyz, model record), (letter, xyz, model record)) as they come in
for k in split:
    (a1, i1), (a2, i2) = sorted(recs[k].items())
    s1, s2 = ref[k][a1], ref[k][a2]
    p1, p2 = xyz(lines[i1]), xyz(lines[i2])
    if d2(p1, s2) + d2(p2, s1) < d2(p1, s1) + d2(p2, s2):
        s1, s2 = s2, s1
    w_sites.append(((a1, s1, i1), (a2, s2, i2)))
    nback += 2

# atoms the waters can touch: everything but the split waters, plus the split-water sites (index >= nfix)
fix = [i for i, l in enumerate(lines) if is_atom(l) and not (is_water(l) and site(l) in split)]
pos = [xyz(out[i]) for i in fix]
alt = [lines[i][16] for i in fix]
ele = [element(lines[i]) for i in fix]
nfix = len(fix)
wat = set(j for j, i in enumerate(fix) if is_water(lines[i]))
for w in w_sites:
    for (a, s, i) in w:
        wat.add(len(pos)); pos.append(s); alt.append(None); ele.append('O')
sites_frac = cs.unit_cell().fractionalize(flex.vec3_double(pos))
asu = cs.asu_mappings(buffer_thickness=CUT)
asu.process_sites_frac(original_sites=sites_frac, min_distance_sym_equiv=0.01)
pairs = crystal.neighbors_fast_pair_generator(asu_mappings=asu, distance_cutoff=CUT)
close = {}
hbond = {}
for p in pairs:
    i, j, d = p.i_seq, p.j_seq, p.dist_sq ** 0.5
    if i == j: continue
    close.setdefault(i, []).append((j, d)); close.setdefault(j, []).append((i, d))
# polar H: an H whose nearest heavy atom (within 1.25 A) is N or O
polar = set()
for i in range(nfix):
    if ele[i] != 'H': continue
    nb = [(d, j) for j, d in close.get(i, []) if j < nfix and ele[j] != 'H' and d < 1.25]
    if nb and ele[min(nb)[1]] in ('N', 'O'): polar.add(i)

# label of each split-water site: lab[w] = 0 keeps the incoming letters, 1 swaps them
lab = [0] * len(w_sites)
def site_index(w, n): return nfix + 2 * w + n
def letter(idx):
    if idx < nfix: return alt[idx]
    w, n = divmod(idx - nfix, 2)
    return w_sites[w][n ^ lab[w]][0]
def sees(a, b): return a == ' ' or b == ' ' or a == b
def cost(w):
    e = 0.0
    for n in (0, 1):
        si = site_index(w, n)
        for j, d in close.get(si, []):
            if j >= nfix and (j - nfix) // 2 == w: continue
            if not sees(letter(si), letter(j)): continue
            d0 = D0_WATER if j in wat else D0_POLAR_H if j in polar else D0.get(ele[j], 3.0)
            if d < d0: e += (d0 - d) ** 2
    return e
for sweep in range(20):
    flips = 0
    for w in range(len(w_sites)):
        e0 = cost(w); lab[w] ^= 1; e1 = cost(w)
        if e1 < e0 - 1e-6: flips += 1
        else: lab[w] ^= 1
    if not flips: break

nswap = 0
nref = 0
for w, ((a1, s1, i1), (a2, s2, i2)) in enumerate(w_sites):
    # site n keeps its model record (occupancy, B) and takes the letter lab chooses
    for n, (a, s, i) in enumerate(((a1, s1, i1), (a2, s2, i2))):
        newa = w_sites[w][n ^ lab[w]][0]
        slot = recs[split[w]][newa]                    # write it where that letter's record was
        x, y, z = s
        out[slot] = lines[slot][:6] + lines[slot][6:11] + lines[i][11:16] + newa + lines[i][17:30] \
                    + '%8.3f%8.3f%8.3f' % (x, y, z) + lines[i][54:]
    k = split[w]
    if lab[w]: nswap += 1
    if (ref[k][a1] == s1) == bool(lab[w]): nref += 1
open(out_file, 'w').write('\n'.join(out) + '\n')
sys.stderr.write('reset_waters: %d waters put back, %d of %d split waters relabelled, %d of them labelled opposite to the reference\n'
                 % (nback, nswap, len(w_sites), nref))
'EOF'
  $python ${tmp}.py "$ref" "$pdb" $tmp
  if( $status || ! -s $tmp ) then
    echo "reset_waters: label-free reset failed, keeping the reference labels" > /dev/stderr
    set labels = same
  endif
  rm -f ${tmp}.py
endif

if( "$labels" != "free" ) then
  awk 'FNR==NR{if(/^ATOM|^HETATM/ && substr($0,18,3)=="HOH")xyz[substr($0,13,15)]=substr($0,31,24);next}\
    /^ATOM|^HETATM/ && substr($0,18,3)=="HOH" && (substr($0,13,15) in xyz){\
      $0=substr($0,1,30) xyz[substr($0,13,15)] substr($0,55);++n}\
    {print}\
    END{print "reset_waters: " n+0 " waters put back" > "/dev/stderr"}' "$ref" "$pdb" >! $tmp
endif
mv $tmp "$outfile"

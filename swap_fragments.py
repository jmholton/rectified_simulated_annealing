# swap_fragments.py: break a two-conformer (A/B) model into swap-together fragments with good internal geometry, then deduce which
# relative A/B flips of those fragments (and of the waters) the coordinates alone allow -- no refinement, no ground truth.
# JH 5 Nov 2026: "write an analyzer for the model, given geometry, to try and break it up into swap-together fragments that have good
# internal geometry. There will have to be a threshold of some sort for splitting. then apply principles of
# ../claude_numbers/swapwater_2026-11-05/water_deduce.py before finally putting it into phenix."
#
#   phenix.python swap_fragments.py MODEL.pdb [cut=25] [hinge=40] [clash=0.6] [sym=0] [singles=1] [pairs=1] [nsample=N seed=S] prefix=OUT
#   (phenix.pdb_interpretation from PATH; numpy)
#
# 1. Fragments.  phenix.pdb_interpretation's restraints for MODEL (bonds, angles, dihedrals, chirality, planes; CDL; both letters)
#    are summed onto the covalent bonds they contain ("edge strain" = sum of the residuals (delta/sigma)^2 of every restraint that
#    contains both atoms of the bond, A and B copies together).  A label crossover between two atoms strains exactly those
#    restraints.  Bonds with strain > CUT are cut (bonds to H never: H follows its heavy atom), and the connected pieces of the
#    remaining covalent graph (disulfides included) are the fragments: inside a fragment every restraint is satisfied under the
#    model's own A/B pairing, so the fragment swaps as a unit.  Each water residue is its own unit.
# 2. Deduction (the water_deduce.py principle, for protein fragments too).  A split copy of MODEL (A copy = chain A, B copy = chain
#    Z, water B copies resseq+1000, all altlocs blank) gives phenix's nonbonded pairs between EVERY two copies, crystal symmetry
#    included.  For a pair of units X != Y with overlap (vdw - distance) > T, atoms more than 3 bonds apart:
#       same letters      -> they clash now; a relative flip of X and Y removes it    (parity 1);
#       different letters -> a relative flip of X and Y would make them clash         (parity 0).
#    Union-find with parity, strongest overlap first; a constraint contradicting stronger ones is a conflict (reported, dropped).
#    A component of units is fixed relative to itself; each component has one free bit (flip as a unit); the whole model's A<->B
#    exchange is a symmetry, so the component holding the largest fragment is the anchor.
#    Bonded links across the cuts are reported (straight vs crossed restraint cost) but are NOT used as constraints: on a refined
#    model the straight pairing has been refined into place, so the crossed cost is large whether the pairing is right or wrong.
#    singles=1: base = the deduced relative flips applied, and one model per free component holding protein flipped alone
#    (OUT_base.pdb, OUT_c<k>.pdb); pairs=1 adds one model per two such components that meet at a cut bond (OUT_p<k>.pdb): a right
#    flip of one fragment is penalised when its neighbour is also wrong (1aho round 15, 5 Nov 2026).
#    Defaults cut=25 hinge=40 clash=0.6 sym=0 from a scan on 1aho (swapfrag_2026-11-05/README.md in the 1aho project).
# 3. Draws (nsample): uniform over the allowed assignments = every free component flipped as a unit with p 0.5 (the anchor never),
#    internal parities as deduced.  Writes OUT_<k>.pdb (letters exchanged in the flipped units) and OUT_<k>.txt.
import sys, os, math, random, subprocess
import numpy as np

other = {'A': 'B', 'B': 'A'}

def opts(argv):
    o = {'cut': 25.0, 'hinge': 40.0, 'clash': 0.6, 'sym': 0, 'singles': 0, 'pairs': 0, 'nsample': 0, 'seed': 0, 'prefix': None, 'geo': None, 'splitgeo': None, 'verbose': 0}
    for a in argv:
        k, v = a.split('=', 1)
        o[k] = type(o[k])(v) if o.get(k) is not None and not isinstance(o[k], str) else v
    return o

def is_water(l): return l[17:20] == 'HOH'

def read_model(path):
    xyz, lines = {}, open(path).readlines()
    for l in lines:
        if l[:6] in ('ATOM  ', 'HETATM'):
            xyz[(site_of(l[12:16], l[17:20], int(l[22:26])), l[16])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return xyz, lines

def site_of(name, resn, resseq):
    return ('W', resseq) if resn == 'HOH' else (resseq, name.strip())

def is_H(site): return site[0] != 'W' and site[1][0] in 'HD' and not site[1].startswith('HG') or (site[0] != 'W' and site[1] in ('HG', 'HG1', 'HG2', 'HG3') or (site[0] != 'W' and site[1].startswith('HG') and len(site[1]) > 2))

def H_test(name, el):  # by element column when present
    return el == 'H' or el == 'D'

def run_interp(pdb, extra=()):
    geo = pdb + '.geo'
    if not os.path.exists(geo):
        cmd = ['phenix.pdb_interpretation', pdb, 'write_geo_files=True', 'const_shrink_donor_acceptor=0.6'] + list(extra)
        with open(pdb + '_interp.log', 'w') as f: subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
    return geo

def pid(s, split=False):
    # pdb=" N  AVAL A   1 " -> (site, letter)
    name, alt, resn, ch, rs = s[0:4], s[4], s[5:8], s[9], int(s[10:14])
    if split:
        if resn == 'HOH': return (('W', rs % 1000), 'B' if rs >= 1000 else 'A')
        return (site_of(name, resn, rs), 'B' if ch == 'Z' else 'A')
    return (site_of(name, resn, rs), alt)

def parse_geo(path, split=False):
    L = open(path).read().splitlines(); R = []; NB = []
    i = 0
    while i < len(L):
        l = L[i]; w = l.split(' ', 1)[0]
        if w in ('bond', 'angle', 'dihedral', 'chirality'):
            n = {'bond': 2, 'angle': 3, 'dihedral': 4, 'chirality': 4}[w]
            ids = [pid(L[i + k].split('"')[1], split) for k in range(n)]
            v = L[i + n + 1].split()
            if w == 'chirality':
                R.append(dict(kind=w, atoms=ids, both=(v[0] == 'True'), ideal=float(v[1]), model=float(v[2]), sigma=float(v[4]), res=float(v[6])))
            elif w == 'dihedral':
                R.append(dict(kind=w, atoms=ids, ideal=float(v[0]), model=float(v[1]), per=int(v[3]), sigma=float(v[4]), res=float(v[6])))
            else:
                R.append(dict(kind=w, atoms=ids, ideal=float(v[0]), model=float(v[1]), sigma=float(v[3]), res=float(v[5])))
            i += n + 2; continue
        if w == 'plane':
            ids, sig = [], []
            k = i
            while L[k].lstrip().startswith(('plane pdb=', 'pdb=')):
                ids.append(pid(L[k].split('"')[1], split)); f = L[k].split('"')[2].split(); sig.append(float(f[1]))
                if len(f) > 3: res = float(f[4])
                k += 1
            R.append(dict(kind='plane', atoms=ids, sigmas=sig, res=res)); i = k + 1; continue
        if w == 'nonbonded':
            a = pid(L[i].split('"')[1], split); b = pid(L[i + 1].split('"')[1], split); v = L[i + 3].split()
            NB.append((a, b, float(v[0]), float(v[1]), v[2] if len(v) > 2 else ''))
            i += 4; continue
        i += 1
    return R, NB

def dihed(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))

def residual(r, X):
    # residual of restraint r with coordinates X (list, same order as r['atoms']); harmonic forms (periodic for dihedrals)
    k = r['kind']
    if k == 'bond': return ((np.linalg.norm(X[0] - X[1]) - r['ideal']) / r['sigma']) ** 2
    if k == 'angle':
        u, v = X[0] - X[1], X[2] - X[1]
        a = math.degrees(math.acos(max(-1, min(1, np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))))
        return ((a - r['ideal']) / r['sigma']) ** 2
    if k == 'dihedral':
        p = max(r['per'], 1); d = (dihed(*X) - r['ideal']) % (360.0 / p)
        d = min(d, 360.0 / p - d); return (d / r['sigma']) ** 2
    if k == 'chirality':
        vol = np.dot(X[1] - X[0], np.cross(X[2] - X[0], X[3] - X[0]))
        if r['both']: vol = abs(vol); ideal = abs(r['ideal'])
        else: ideal = r['ideal']
        return ((vol - ideal) / r['sigma']) ** 2
    if k == 'plane':
        P = np.array(X); w = 1.0 / np.array(r['sigmas']) ** 2
        c = (P * w[:, None]).sum(0) / w.sum(); Q = (P - c) * np.sqrt(w)[:, None]
        n = np.linalg.svd(Q)[2][-1]
        return float((((P - c) @ n) ** 2 * w).sum())

class UF:  # union-find with parity (from water_deduce.py)
    def __init__(self): self.p, self.x = {}, {}
    def find(self, a):
        if a not in self.p: self.p[a], self.x[a] = a, 0
        if self.p[a] == a: return a, 0
        r, px = self.find(self.p[a]); self.p[a] = r; self.x[a] ^= px; return r, self.x[a]
    def union(self, a, b, par):
        (ra, xa), (rb, xb) = self.find(a), self.find(b)
        if ra == rb: return (xa ^ xb) == par
        self.p[rb] = ra; self.x[rb] = xa ^ xb ^ par; return True

def write_split(lines, out):
    o = []
    for l in lines:
        if l[:6] in ('ATOM  ', 'HETATM'):
            alt = l[16]
            if is_water(l):
                if alt == 'B': l = l[:16] + ' ' + l[17:22] + '%4d' % (int(l[22:26]) + 1000) + l[26:]
                else: l = l[:16] + ' ' + l[17:]
                o.append(l)
            else:
                for c, ch in (('A', 'A'), ('B', 'Z')):
                    if alt in (c, ' '): o.append(l[:16] + ' ' + l[17:21] + ch + l[22:])
        elif l[:6] in ('ANISOU', 'TER   ', 'MODEL ', 'ENDMDL', 'CONECT'):
            continue
        else: o.append(l)
    # chain A first, then chain Z, then waters, each followed by TER
    head = [l for l in o if l[:6] not in ('ATOM  ', 'HETATM', 'END   ', 'END\n')]
    atoms = [l for l in o if l[:6] in ('ATOM  ', 'HETATM')]
    A = [l for l in atoms if l[21] == 'A' and not is_water(l)]; Z = [l for l in atoms if l[21] == 'Z']; W = [l for l in atoms if is_water(l)]
    open(out, 'w').writelines(head + A + ['TER\n'] + Z + ['TER\n'] + W + ['END\n'])

def main():
    model = sys.argv[1]; o = opts(sys.argv[2:])
    xyz, lines = read_model(model)
    el = {}
    for l in lines:
        if l[:6] in ('ATOM  ', 'HETATM'): el[site_of(l[12:16], l[17:20], int(l[22:26]))] = l[76:78].strip()
    isH = lambda s: s[0] != 'W' and el.get(s) in ('H', 'D')
    R, NB0 = parse_geo(o['geo'] or run_interp(model))
    # covalent graph over sites (letter-independent)
    adj = {}
    for r in R:
        if r['kind'] == 'bond':
            a, b = r['atoms'][0][0], r['atoms'][1][0]
            adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    sites = sorted({s for (s, L) in xyz if s[0] != 'W'}, key=lambda s: (s[0], s[1]))
    for s in sites: adj.setdefault(s, set())
    edges = {tuple(sorted((a, b))) for a in adj for b in adj[a]}
    strain = dict.fromkeys(edges, 0.0)
    for r in R:
        S = {s for s, L in r['atoms']}
        for e in edges:
            pass
    # faster: for each restraint, the edges among its sites
    for r in R:
        S = list({s for s, L in r['atoms']})
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                e = tuple(sorted((S[i], S[j])))
                if e in strain: strain[e] += r['res']
    heavy_edges = [e for e in edges if not isH(e[0]) and not isH(e[1])]
    # crossed cost of each bond: the same restraints with the letters of the atoms on one side of the bond exchanged (the side is
    # taken within each restraint's own atoms; a restraint that stays connected without the bond, a ring, is left out)
    crossed = dict.fromkeys(edges, 0.0); straight_open = dict.fromkeys(edges, 0.0)
    for r in R:
        S = [s for s, L in r['atoms']]; Sset = set(S)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                e = tuple(sorted((S[i], S[j])))
                if e not in crossed or S[j] not in adj[S[i]]: continue
                side = {S[j]}; st = [S[j]]
                while st:
                    a = st.pop()
                    for b in adj[a]:
                        if b in Sset and b not in side and not (a == S[j] and b == S[i]):
                            side.add(b); st.append(b)
                if S[i] in side: continue
                keys = [(s, other[L] if s in side else L) for s, L in r['atoms']]
                if any(k not in xyz for k in keys): continue      # a lone atom (His54 tautomer H) has no partner copy
                crossed[e] += residual(r, [xyz[k] for k in keys]); straight_open[e] += r['res']
    bind = {e: crossed[e] - straight_open[e] for e in edges}
    def nbr(e): return [f for f in heavy_edges if f != e and (e[0] in f or e[1] in f)]
    # strained bonds: strain > cut and the largest among the bonds sharing an atom with it (one cut per crossover, not a smear)
    cut = {e for e in heavy_edges if strain[e] > o['cut'] and all(strain[e] >= strain[f] for f in nbr(e))}
    # hinges: the crossed pairing costs less than hinge more than the straight one (A and B nearly coincide around the bond, so
    # geometry cannot tell the pairings apart); smallest among its neighbours
    hinge = set()
    if o['hinge'] > 0:
        hinge = {e for e in heavy_edges if bind[e] < o['hinge'] and all(bind[e] <= bind[f] for f in nbr(e))}
    cut |= hinge
    # fragments
    frag = {}; n = 0
    for s in sites:
        if s in frag: continue
        n += 1; stack = [s]; frag[s] = n
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if b not in frag and tuple(sorted((a, b))) not in cut:
                    frag[b] = n; stack.append(b)
    unit = dict(frag)
    for (s, L) in xyz:
        if s[0] == 'W': unit[s] = 'W%d' % s[1]
    members = {}
    for s, f in frag.items(): members.setdefault(f, []).append(s)
    def desc(f):
        rs = sorted({s[0] for s in members[f]}); heavy = [s for s in members[f] if not isH(s)]
        runs, st = [], rs[0]
        for a, b in zip(rs, rs[1:] + [None]):
            if b != a + 1: runs.append('%d' % st if st == a else '%d-%d' % (st, a)); st = b
        part = '' if len(heavy) > 12 else ' [' + ' '.join('%s%d' % (s[1], s[0]) for s in sorted(heavy)) + ']'
        return '%s (%d heavy)%s' % (','.join(runs), len(heavy), part)
    order = sorted(members, key=lambda f: -len(members[f]))
    print('MODEL %s  cut=%.1f  hinge=%.1f  clash=%.2f' % (model, o['cut'], o['hinge'], o['clash']))
    print('edge strain (heavy-atom bonds, sum of residuals of the restraints containing the bond, A+B): top 15')
    for e in sorted(heavy_edges, key=lambda e: -strain[e])[:15]:
        print('   %7.1f  %s%d-%s%d  crossed-straight %8.1f%s' % (strain[e], e[0][1], e[0][0], e[1][1], e[1][0], bind[e], '  CUT' if e in cut else ''))
    print('bonds with the smallest crossed-straight (hinge candidates): top 25')
    for e in sorted(heavy_edges, key=lambda e: bind[e])[:25]:
        print('   %7.1f  %s%d-%s%d  strain %.1f%s' % (bind[e], e[0][1], e[0][0], e[1][1], e[1][0], strain[e], '  HINGE' if e in hinge else ''))
    print('fragments: %d (%d cut bonds); largest first:' % (len(members), len(cut)))
    for f in order: print('   F%-3d %s' % (f, desc(f)))
    # bonded links across cuts: straight vs crossed cost
    print('cut links (restraints spanning two fragments; residual sum A+B): straight -> crossed')
    pairs = {}
    for r in R:
        F = {frag[s] for s, L in r['atoms'] if s in frag}
        if len(F) == 2: pairs.setdefault(tuple(sorted(F)), []).append(r)
    for (f1, f2), rs in sorted(pairs.items()):
        st = sum(r['res'] for r in rs); cr = 0.0
        for r in rs:
            keys = [(s, other[L] if frag[s] == f2 else L) for s, L in r['atoms']]
            cr += residual(r, [xyz[k] for k in keys]) if all(k in xyz for k in keys) else r['res']
        print('   F%d|F%d  %d restraints  %.1f -> %.1f' % (f1, f2, len(rs), st, cr))
    # nonbonded constraints from the split model
    sp = model[:-4] + '_split.pdb'
    if o['splitgeo'] is None:
        if not os.path.exists(sp): write_split(lines, sp)
        sgeo = run_interp(sp, ['clash_guard.nonbonded_distance_threshold=None', 'link_distance_cutoff=0', 'disulfide_distance_cutoff=0'])
    else: sgeo = o['splitgeo']
    _, NB = parse_geo(sgeo, split=True)
    # topological distance <= 3 between sites: excluded (bonded geometry, not contacts)
    near = {}
    for s in adj:
        seen = {s: 0}; fr = [s]
        for d in (1, 2, 3):
            nf = []
            for a in fr:
                for b in adj[a]:
                    if b not in seen: seen[b] = d; nf.append(b)
            fr = nf
        near[s] = seen
    cons, internal, selfsym = [], [], []
    for (sa, La), (sb, Lb), d, vdw, sym in NB:
        ov = vdw - d
        if ov <= o['clash'] or sa == sb: continue
        if sa[0] != 'W' and sb[0] != 'W' and (sb in near.get(sa, {}) or sa[0] == sb[0]) and not sym: continue  # bonded or same residue
        ua, ub = unit[sa], unit[sb]
        if ua == ub:
            if La == Lb: (selfsym if sym else internal).append((ov, sa, La, sb, Lb, sym))
            continue
        if sym and not o['sym']: continue   # crystal contacts: off by default (r15 test: 4 of 5 pointed the wrong way)
        cons.append((ov, ua, ub, 1 if La == Lb else 0, sa, La, sb, Lb, sym))
    cons.sort(key=lambda c: -c[0])
    uf = UF(); conflicts = []; used = 0; viol = []
    allunits = sorted(set(unit.values()), key=str)
    for u in allunits: uf.find(u)
    status = []
    for c in cons:
        if c[3] == 1: viol.append(c)
        if uf.union(c[1], c[2], c[3]): used += 1; status.append('used')
        else: conflicts.append(c); status.append('conflict')
    comp = {}
    for u in allunits: comp.setdefault(uf.find(u)[0], []).append(u)
    def lab(u): return ('S%s' % u[1:]) if isinstance(u, str) else 'F%d' % u
    def atom(s, L): return ('S%d%s' % (s[1], L)) if s[0] == 'W' else '%s%d%s' % (s[1], s[0], L)
    print('nonbonded constraints (overlap > %.2f, >3 bonds apart): %d (clash now -> flip: %d, flip would clash -> keep: %d); '
          'contradicted by stronger ones: %d' % (o['clash'], len(cons), len(viol), len(cons) - len(viol), len(conflicts)))
    for c in viol[:12]: print('   clash now: %.2f %s %s  (%s|%s) %s' % (c[0], atom(c[4], c[5]), atom(c[6], c[7]), lab(c[1]), lab(c[2]), c[8]))
    for c in conflicts[:12]: print('   conflict:  %.2f %s %s  (%s|%s) parity %d %s' % (c[0], atom(c[4], c[5]), atom(c[6], c[7]), lab(c[1]), lab(c[2]), c[3], c[8]))
    for c in sorted(internal, reverse=True)[:8]: print('   inside one unit: %.2f %s %s %s' % (c[0], atom(c[1], c[2]), atom(c[3], c[4]), c[5]))
    for c in sorted(selfsym, reverse=True)[:4]: print('   unit vs its own symmetry mate: %.2f %s %s %s' % (c[0], atom(c[1], c[2]), atom(c[3], c[4]), c[5]))
    anchor = uf.find(order[0])[0]
    comps = sorted(comp.values(), key=lambda c: -len(c))
    prot = [c for c in comps if any(not isinstance(u, str) for u in c)]
    print('components: %d (%d holding protein); free bits (anchor = component of F%d): %d' % (len(comps), len(prot), order[0], len(comps) - 1))
    for c in prot:
        fl = [u for u in c if uf.find(u)[1] != uf.find(c[0])[1]]
        print('   {%s}  relative flips inside: %s' % (' '.join(lab(u) for u in sorted(c, key=str)), ' '.join(lab(u) for u in sorted(fl, key=str)) or '-'))
    if o['prefix']:
        with open(o['prefix'] + '.cons', 'w') as F:
            for c, st_ in zip(cons, status):
                F.write('%.3f %s %s %d %s %s %s %s\n' % (c[0], lab(c[1]), lab(c[2]), c[3], atom(c[4], c[5]), atom(c[6], c[7]), st_, c[8] or '-'))
    # per-site table: resseq name fragment component parity(deduced flip relative to the component's root) anchor|free
    if o['prefix']:
        with open(o['prefix'] + '.frag', 'w') as F:
            for s_ in sites:
                u = unit[s_]; r_, x_ = uf.find(u)
                F.write('%d %s F%d C%s %d %s\n' % (s_[0], s_[1], u, lab(r_), x_, 'anchor' if r_ == anchor else 'free'))
    # the deduced state relative to the model: within each component, flip the side opposite to its largest member
    size = lambda u: len(members[u]) if not isinstance(u, str) else 0
    def write_flipped(flip, path, note):
        out = []
        for l in lines:
            if l[:6] in ('ATOM  ', 'HETATM', 'ANISOU') and l[16] in 'AB' and flip.get(unit.get(site_of(l[12:16], l[17:20], int(l[22:26])))):
                l = l[:16] + other[l[16]] + l[17:]
            out.append(l)
        open(path, 'w').writelines(out)
        fl = sorted((u for u in flip if flip[u]), key=str)
        open(path[:-4] + '.txt', 'w').write('%s | %d units flipped: %s\n' % (note, len(fl), ' '.join(lab(u) for u in fl)))
    if o['singles']:
        # base = the deduced relative flips applied (inside every component, the units whose parity differs from its largest unit's);
        # then one model per component holding protein (not the anchor's): base + that whole component flipped
        base = {}
        for c in comps:
            big = max(c, key=size)
            for u in c: base[u] = uf.find(u)[1] ^ uf.find(big)[1]
        write_flipped(base, o['prefix'] + '_base.pdb', 'base (deduced relative flips)')
        k = 0
        for c in comps:
            if not any(not isinstance(u, str) for u in c) or uf.find(c[0])[0] == anchor: continue
            k += 1; fl = dict(base)
            for u in c: fl[u] = 1 - base[u]
            write_flipped(fl, '%s_c%d.pdb' % (o['prefix'], k), 'component {%s} flipped' % ' '.join(lab(u) for u in sorted(c, key=str)))
        print('singles: %s_base.pdb + %d single-component flips %s_c<k>.pdb' % (o['prefix'], k, o['prefix']))
        if o['pairs']:
            croot = {u: uf.find(u)[0] for u in allunits}
            adjc = set()
            for a, b in cut:
                ca, cb = croot[frag[a]], croot[frag[b]]
                if ca != cb and anchor not in (ca, cb): adjc.add(tuple(sorted((ca, cb), key=str)))
            k = 0
            for ca, cb in sorted(adjc, key=str):
                k += 1; fl = dict(base)
                for u in comp[ca] + comp[cb]: fl[u] = 1 - base[u]
                write_flipped(fl, '%s_p%d.pdb' % (o['prefix'], k), 'components {%s} + {%s} flipped' % (
                    ' '.join(lab(u) for u in sorted(comp[ca], key=str)), ' '.join(lab(u) for u in sorted(comp[cb], key=str))))
            print('pairs: %d flips of two adjacent components %s_p<k>.pdb' % (k, o['prefix']))
    if o['nsample']:
        rng = random.Random(o['seed'])
        for k in range(1, o['nsample'] + 1):
            flip = {}
            for c in comps:
                big = max(c, key=size); f = 0 if uf.find(big)[0] == anchor else int(rng.random() < 0.5)
                for u in c: flip[u] = uf.find(u)[1] ^ uf.find(big)[1] ^ f
            out = []
            for l in lines:
                if l[:6] in ('ATOM  ', 'HETATM', 'ANISOU') and l[16] in 'AB' and flip.get(unit.get(site_of(l[12:16], l[17:20], int(l[22:26])))):
                    l = l[:16] + other[l[16]] + l[17:]
                out.append(l)
            open('%s_%d.pdb' % (o['prefix'], k), 'w').writelines(out)
            fl = sorted((u for u in flip if flip[u]), key=str)
            open('%s_%d.txt' % (o['prefix'], k), 'w').write('%d %s\n' % (len(fl), ' '.join(lab(u) for u in fl)))

if __name__ == '__main__':
    main()

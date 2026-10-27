# rectified_simulated_annealing

its regular simulated annealing, except any atom that ends up in a worse place gets put back to where it started

## Motivation

Simulated annealing (SA) refinement is a great way to help your model escape from a local minimum, but it is also a great way to put your model into a different local minimum at the same time. Fortunately, local minima are often "local", and characterized by strained chemical geometry restraints and/or poor fit to the density. This allows each atom to be evaluated before and after the SA run to see if the operation made it fit better or worse. If it was worse, why not just "cancel" the SA for that particular atom by moving it back to where it was before the run? And if an atom fits better, leave it where it is. 

## Description

 This script
 
<br>



## Getting Started

### Dependencies

* These are linux c-shell scripts, so you will need a working `/bin/tcsh`
* the SA and other operations are performed using the Phenix Suite, which must be installed.
* The helper scripts `molprobify_runme.com`, `combine_pdbs_runme.com`, `rmsd`, `memory_repair_runme.com` and `reset_waters_runme.com` ship with this repo.

### Installing

* git clone this repo
* put the repo directory in your shell `$PATH` (or copy all the scripts, `rmsd` included, into one directory that is), and make them executable:
```
    chmod u+x *.com
```
Yes, I know the extension says `*.com`, but these are not Windows executables. The use of `.com` to denote shell scripts pre-dates Windows.

The parent puts its own directory at the front of the path it hands to the trials (`$RSA_HOME`), so the helpers in the repo are used even if other copies are in your `$PATH`.

### Executing program

#### test just one seed
For best results you will want to try ~1000 random number seeds, but to start, just run it with one to make sure it works on your system. Now, this script relies on well-known file names, so you need to copy your input pdb and mtz files:
```
cp refined_001.pdb starthere.pdb
cp my_data.mtz refme.mtz
rectified_SA_runme.com seed=1 debug=1 temperature=300
```
This will serve as a good test to see if everything is working. Success is indicated by the output of some files called smooth_vs_X.txt, which is the smoothed version of the X-coordinate of atom number 123 in all the files. The parallel parent script below uses these text files to make the new PDBs. Other command-line options:
- seed    random number seed to give to the phenix.refine SA job. default: 1
- temperature     simulation temperature to use, 2000 to 5000K is typical, but default is 300
- noise   the sigma deviate in density or geometry score that is considered small enough to be noise. default: 3
- cbsigma  per-component sigma (A) of a C-beta deviation. The deviation is a 3D distance, so dev/cbsigma is scored as a chi deviate with 3 degrees of freedom, converted to the Gaussian deviate with the same tail probability; a residue conformer whose C-beta deviation got worse and is beyond `noise` sigma is put back (0.135 A at the defaults). default: 0.032, the rms/sqrt(3) of the C-beta deviations of correct models
- water_sigma  hold the waters: before every refinement each water is put back at its position in `water_reference.pdb` (written by the parallel parent from the first remembered model; a trial run on its own uses its own `starthere.pdb`) and held there by a harmonic restraint of this sigma (A) through the simulated annealing and every refinement. wE has no term that holds a water, and above ~5000 K the annealing scatters them (1aho: ~100 of 154 waters > 1 A at 6000 K, none with the restraint). 0 = off. default: 0.1
- debug    keep all intermediate files for debugging purposes

#### the parallel parent script
Run `rectified_SA_runme.com parallel=1` in a directory containing `starthere.pdb`, `refme.mtz` and `opts.eff`. Each round runs 80 trials (1000 to 8000 K, 10 seeds each, via `srun`), takes the best eligible trial as the next `starthere.pdb` if its wE is lower, and stops after 10 rounds without improvement or when a file called `exit` appears. Options for the parent:
- noise, mapnoise, cbsigma, water_sigma  passed to the trials (defaults: 3.5, 3.5, 0.032, 0.1); the repair refinement holds the waters too
- nseeds    seeds per temperature (default: 10)
- temperatures  comma-separated list (default: 1000,2000,...,8000)
- memory    keep every accepted model in `memory/` and check each round's winner against all of them (default: 1)
- rfree_slack   a trial is eligible only if its Rfree (in %) is no more than this above the lowest Rfree of any remembered model; negative means no limit (default: 0.2)

### Memory of good models
Each trial only compares atoms with the start of the current round, so two things can slip through one small step at a time:
- **the fit to the data:** ranking trials by wE alone will happily trade R factor for geometry, round after round. In one 12-round run from a model with Rfree 4.4 %, wE went from 105 to 19 while Rfree went to 9.6 %.
- **single atoms:** an atom can slide out of its density a little each round, because each round is tested against an already crept start. In that run a lysine NZ ended 2.4 A from where it belongs, sitting in -4.5 sigma difference density with +16.8 sigma at its old position.

With `memory=1` every accepted model is copied into `memory/round<N>.pdb` and listed in `memory/ledger.txt` (round, wE, Rwork, Rfree, file). The Rfree gate (`rfree_slack`) is measured from the lowest Rfree in memory, not from the current start, so the R factor cannot creep. Then `memory_repair_runme.com` looks at the winner's mFo-DFc map: any non-hydrogen atom below -mapnoise sigma, where the same atom in a remembered model sits above -mapnoise and at least 1 sigma higher, has its whole residue conformer copied back from the remembered model the map likes best. The repaired model gets the same 10-cycle final refinement as a trial, and replaces the winner if it still passes the Rfree gate. Log: `repair_<round>.log`; refinement files: `repairs/`.

`memory_repair_runme.com` also runs on its own:
```
memory_repair_runme.com model.pdb model.mtz memory/round*.pdb mapnoise=3.5 outfile=repaired.pdb
```



## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.


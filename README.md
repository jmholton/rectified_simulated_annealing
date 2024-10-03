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
* The molprobify_runme.com script must also be installed.

### Installing

* git clone this repo
* copy the *.com files into somewhere in your shell `$PATH`, and make them executable:
```
    chmod u+x *.com
```
Yes, I know the extension says `*.com`, but these are not Windows executables. The use of `.com` to denote shell scripts pre-dates Windows.

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
- debug    keep all intermediate files for debugging purposes



## Help

Let me know if you have any questions or find any bugs.  In general, debugging information is provided by adding the command-line option: `debug=1`, and possibly `tempfile=temp` to change the default temporary file prefix.


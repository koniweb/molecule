molecule
=========

molecule is a python class to read in/write out xyz,lammps,pwscf files.  
Additionally it allows you to modify the input data or use the data to
do analysis.  

It is in principle the basis for analysis codes.  
  
It contains:  

file              | comment
:-----------------|:-----------------------------------------
class_molecule.py | main module file with molecule class
mod_calc.py       | module with calculation functions
mod_extend.py	  | module to extend the molecule
mod_lmp.py	      | modules to read and write files (LAMMPS)
mod_pwscf.py	  | modules to read and write files (pwscf input and output)
mod_xyz.py        | modules to read and write files (xyz)

additional module files
=========================

file                      | comment
:-------------------------|:-----------------------------------------
class_molecule_crystal.py | extention to class_molecule.py to cut nanocrystals

To Do
=====
- [X] add additional read write flags for lammps in and output
- [X] allow for extended xyz output
- [X] allow extended xyz input
- [X] combine extended and normal xyz read
- [X] combine extended and normal xyz write
- [X] extended xyz write also take stored data of atoms in molecule
- [X] add rewrapping function for atoms outside the box
- [ ] pwscf read output -- redo the write of coordinates - Final Coordinates
- [ ] correct additional numbering in xyz output
- [S] change to numpy arrays
- [ ] use f2py to use fortran code for expensive calculations

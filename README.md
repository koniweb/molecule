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
class_molecule_crystal.py -- extention to class_molecule.py to cut nanocrystals

To Do
=====
- [ ] correct additional numbering in xyz output
- [X] add additional read write flags for lammps in and output
- [ ] allow for extended xyz output

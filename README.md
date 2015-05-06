molecule
=========

molecule is a python class to read in/write out xyz,lammps,pwscf files.  
Additionally it allows you to modify the input data or use the data to
do analysis.  

It is in principle the basis for analysis codes.  
  
It contains:  

file                | comment
:-------------------|:-----------------------------------------
class_molecule.py   | main module file with molecule class
mod_calc.py         | module with calculation functions
mod_extend.py	    | module to extend the molecule
mod_lmp.py	        | modules to read and write files (LAMMPS)
mod_pwscf.py	    | modules to read and write files (pwscf input and output)
mod_xyz.py          | modules to read and write files (xyz)
fortran_modules.f90 | modules to calculate bonds

additional module files
=========================

file                      | comment
:-------------------------|:-----------------------------------------
class_molecule_crystal.py | extention to class_molecule.py to cut nanocrystals

Compile fortran module
======================
> f2py -c fortran_modules.f90 -m fortran_modules

To Do
=====
- [x] add additional read write flags for lammps in and output
- [x] allow for extended xyz output
- [x] allow extended xyz input
- [x] combine extended and normal xyz read
- [x] combine extended and normal xyz write
- [x] extended xyz write also take stored data of atoms in molecule
- [x] add rewrapping function for atoms outside the box
- [x] pwscf read output -- redo the write of coordinates - Final Coordinates
- [ ] correct additional numbering in xyz output
- [ ] change to numpy arrays
- [x] use f2py to use fortran code for expensive calculations
- [x] use fortran with f2py for bond calculation
- [ ] change typelist to namelist
- [ ] check offset in wrapping

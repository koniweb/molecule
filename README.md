molecule
========

molecule is a python class to read in/write out xyz,lammps,pwscf files.
Additionally it allows you to modify the input data or use the data to do analysis.

It is in principle the basis for analysis codes.

It contains:
class_molecule.py -- main module file with molecule class
mod_calc.py       -- module with calculation functions
mod_extend.py	  -- module to extend the molecule
mod_lmp.py	  -- modules to read and write files (LAMMPS)
mod_pwscf.py	  --  &&                             (pwscf input and output)
mod_xyz.py        --  &&                             (xyz)

class_molecule.py -- extention to class_molecule.py to cut nanocrystals
 
# TODO
change vec to __vec
change offset to __offset
change Rvec() to vec()
change Roffset() to offset()
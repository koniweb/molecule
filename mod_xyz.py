#!/usr/bin/env python
##########################################################################
# writexyz
# readfxyz
##########################################################################
version=3.3
versiontext='# mod_xyz.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import copy

#----------------------------------------------------------------------
# module mod_xyz
#----------------------------------------------------------------------
class molecule_rw:
    # write xyz file
    def writexyz(self,filename="",status='w'):
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        mol=self
        print >>f, mol.natoms()
        print >>f, mol.comment()
        for cntat in range(0,mol.natoms()):
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                    mol.at()[cntat].type()[0], 
                    mol.at()[cntat].coord()[0], 
                    mol.at()[cntat].coord()[1], 
                    mol.at()[cntat].coord()[2]
                    )
                )
        if filename != "": f.close()
        return
    
    # write extended xyz file
    def writeexyz(self,filename="",status='w',data=[]):
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
            
        mol=self
        # data contains a list of data for each atom and the text
        # -> do more convenient
        if len(data)==0:
            data=["",[""]*self.natoms()]
        # add vectors to commentline
        a=""
        b=""
        c=""
        if mol.vec()[0]!=[0.0,0.0,0.0]:
            a="a {:15.10f} {:15.10f} {:15.10f}".format(
                mol.vec()[0][0],mol.vec()[0][1],mol.vec()[0][2])
        if mol.vec()[1]!=[0.0,0.0,0.0]:
            b="b {:15.10f} {:15.10f} {:15.10f}".format(
                mol.vec()[1][0],mol.vec()[1][1],mol.vec()[1][2])
        if mol.vec()[2]!=[0.0,0.0,0.0]:
            c="c {:15.10f} {:15.10f} {:15.10f}".format(
                mol.vec()[2][0],mol.vec()[2][1],mol.vec()[2][2])
        # add energy to comment line
        E=""
        if mol.energy()!=0.0:E="E {:.15f}".format(mol.energy())
        # print output
        print >>f, mol.natoms()
        print >>f, "{:s} {:s} {:s} {:s} {:s}".format(a,b,c,E,data[0]) 
        for cntat in range(0,mol.natoms()):
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f} {:s}'.format(
                    mol.at()[cntat].type()[0], 
                    mol.at()[cntat].coord()[0], 
                    mol.at()[cntat].coord()[1], 
                    mol.at()[cntat].coord()[2],
                    data[1][cntat]
                    )
                )
        if filename != "": f.close()
        return

    # read molecules in xyz file
    def readxyz(self,filename,start=1,end=-1):
        # only last molecule via start=-1 and end=-1
        # set molecule
        molecules=[]
        # check read file
        try: 
            file=open(filename, 'r')
        except IOError:
            print >> sys.stderr, "... input file not found"
            exit()
        # read file
        cntline=0
        oldline=0
        natoms=0
        cntmol=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            # create new molecule
            if (cntline-oldline)%(natoms+2)==1:
                mol=self.__class__()
                mol.clear_atoms()
                cntat=0
                natoms=int(linesplit[0])
            # save comment
            if (cntline-oldline)%(natoms+2)==2:
                comment=line.strip('\n')
            # read atoms
            if ((cntline-oldline)%(natoms+2) >= 3 or 
                (cntline-oldline)%(natoms+2) == 0 ):
                # check if number or atomtype given
                if linesplit[0].isdigit():
                    number=int(linesplit[0])
                    name=linesplit[0]
                else:
                    name=linesplit[0]
                    number=self.name2element(name)[1]

                # append atoms
                mol.append_atom(
                    self.__class__.atom(
                        mol,
                        cntat,
                        name,
                        number,
                        float(linesplit[1]),
                        float(linesplit[2]),
                        float(linesplit[3])
                        )
                    )
                cntat+=1
            # finish molecule and append to list
            if (cntline-oldline)%(natoms+2)==0: 
                cntmol+=1
                oldline=cntline
                if start!=-1 and (  
                    (cntmol>=start and (cntmol<=end or end==-1))  ):
                    mol.set(filename,cntmol,comment)
                    molecules.append(copy.copy(mol))
        # if start==-1 add last frame only
        if start==-1:
            mol.set(filename,cntmol,comment)
            molecules.append(copy.copy(mol))
        # close file
        file.close()
        # return molecules
        return molecules
    

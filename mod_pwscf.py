#!/usr/bin/env python
##########################################################################
# writepwscf
# readfpwscfin
# readfpwscfout
##########################################################################
version=1.1
versiontext='# mod_pwscf.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import string
import copy
import mod_calc as calc
from numpy import matrix
from numpy import linalg

#----------------------------------------------------------------------
# module mod_xyz
#----------------------------------------------------------------------
class molecule_rw:
    # write pwscf output format
    def writepwscf(self,filename="",status='w'):
        # calculate relative coordinates
        vecM = matrix( [ [self.vec()[0][0],self.vec()[1][0],self.vec()[2][0]],
                         [self.vec()[0][1],self.vec()[1][1],self.vec()[2][1]],
                         [self.vec()[0][2],self.vec()[1][2],self.vec()[2][2]] ])
        for atom in self.at():
            tmp = matrix([ [atom.coord()[0]],
                           [atom.coord()[1]],
                           [atom.coord()[2]] ])
            res=linalg.solve(vecM, tmp)
            atom.coord_rel=[float(res[0]),float(res[1]),float(res[2])]
        # set options if not set already
        if not hasattr(self, 'setup_pwscf'):
            self.setup_pwscf=self.SETUP_PWSCF()
        # set option celldm in bohr
        self.setup_pwscf.set_celldm(self.celldm_vec()[0]/calc.b2A)
        self.setup_pwscf.set_nat(self.natoms())
        self.setup_pwscf.set_ntyp(self.ntypes())
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        # print options if present
        if hasattr(self.setup_pwscf, 'control'):
            print >>f, "&control"
            for i in range(len(self.setup_pwscf.control)):
                print >>f, ('{:s}={:s}').format(
                    self.setup_pwscf.control[i][0],
                    self.setup_pwscf.control[i][1])
            print >>f, "/"
        if hasattr(self.setup_pwscf, 'system'):
            print >>f, "&system"
            for i in range(len(self.setup_pwscf.system)):
                print >>f, ('{:s}={:s}').format(
                    self.setup_pwscf.system[i][0],
                    self.setup_pwscf.system[i][1])
            print >>f, "/"
        if hasattr(self.setup_pwscf, 'electrons'):
            print >>f, "&electrons"
            for i in range(len(self.setup_pwscf.electrons)):
                print >>f, ('{:s}={:s}').format(
                    self.setup_pwscf.electrons[i][0],
                    self.setup_pwscf.electrons[i][1])
            print >>f, "/"
        if hasattr(self.setup_pwscf, 'ions'):
            print >>f, "&ions"
            for i in range(len(self.setup_pwscf.ions)):
                print >>f, ('{:s}={:s}').format(
                    self.setup_pwscf.ions[i][0],
                    self.setup_pwscf.ions[i][1])
            print >>f, "/"
        # print atomtypes
        print >>f
        print >>f, "ATOMIC_SPECIES"
        for type in self.typelist():
            # check for local name
            if type[1]=="": name=self.pse()[type[0]][0]
            else:           name=type[1]
            print >>f, ('{:5s} {:11.8f} {:s}').format(
                 name,
                 self.pse()[type[0]][2],
                 self.pse()[type[0]][3])
        print >>f
        # print relative coordinates
        print >>f
        print >>f, "ATOMIC_POSITIONS crystal"
        for atom in self.at():
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                    atom.type()[0],
                    atom.coord_rel[0], 
                    atom.coord_rel[1], 
                    atom.coord_rel[2]
                    )
                )
        print >>f
        print >>f, "CELL_PARAMETERS"
        for cntvec in range(3):
            print >>f , (
                '{:15.10f} {:15.10f} {:15.10f}'.format(
                    self.celldm_vec()[1][cntvec][0], 
                    self.celldm_vec()[1][cntvec][1], 
                    self.celldm_vec()[1][cntvec][2]
                    )
                )
        if (hasattr(self.setup_pwscf,"kpoints")):
            print >>f
            print >>f, "K_POINTS"
            print >>f , ('{:s}'.format(self.setup_pwscf.kpoints))
        f.close()
        return
    
    # read molecules in pwscf input file format
    def readpwscfin(self,filename):
        # set molecule
        molecules=[]
        natoms=0
        ntypes=0
        vec=[]
        # read file
        file=open(filename, 'r')
        opt=""
        cntline=0
        cntvec=0
        cntat=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            #
            # create new molecule
            #
            if cntline==1:
                mol=self.__class__()
                mol.clear_atoms()
                # additional optional fields
                mol.setup_pwscf=mol.SETUP_PWSCF()
            #
            # Do READ IN
            #
            # read vectors
            if  opt=="readvec":
                vec.append([float(linesplit[0]),
                            float(linesplit[1]),
                            float(linesplit[2])])
                cntvec+=1
                if cntvec==3: opt=""
            # read species
            if  opt=="readspecies":
                if not len(linesplit)==0:
                    name=linesplit[0]
                    mass=float(linesplit[1])
                    pot=linesplit[2]
                    mol.set_element(name,mass,pot)
                    ntypes+=1
                else: opt=""
            # read atoms
            elif opt=="readcoord":
                if not len(linesplit)==0:
                    mol.append_atom(
                        self.__class__.atom(
                            mol,
                            cntat,
                            linesplit[0],-1,
                            float(linesplit[1]),
                            float(linesplit[2]),
                            float(linesplit[3])
                            )
                        )
                    # cntat
                    cntat+=1
                    natoms=cntat
                else: opt=""
            # read input options
            opt=mol.readpwscfin_opt(linesplit,opt)
            # read main options
            opt=mol.readpwscfin_blocks(linesplit,opt)
        #
        # set data
        #
        # set celldm
        set=False
        for i in range(len(mol.setup_pwscf.system)):
            if mol.setup_pwscf.system[i][0]=="celldm(1)": 
                celldm=float(mol.setup_pwscf.system[i][1])
                set=True
        if set==False:
            celldm=float(1.0)
        mol.set_celldm(celldm*calc.b2A)    # in A
        mol.setup_pwscf.set_celldm(celldm) # in bohr
        # set vector
        mol.set_vecs(vec[0],vec[1],vec[2])
        # set real coordinates
        mol.rel2real()
        # set molecule and append
        mol.set(filename,1,"")
        molecules.append(copy.copy(mol))
        # close file
        file.close()
        # return molecules
        return molecules

    # read molecules in pwscf output file format
    def readpwscfout(self,filename):
        # set molecule
        molecules=[]
        natoms=0
        vec=[]
        # read file
        file=open(filename, 'r')
        opt=""
        cntline=0
        oldline=0
        cntat=0
        cntvec=0
        cnttypes=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            #
            # create new molecule
            #
            if  cntline==1:
                mol=self.__class__()
                mol.clear_atoms()
                # additional optional fields
                mol.setup_pwscf=mol.SETUP_PWSCF()        
            if  opt=="readcoord" and cntat==0:
                if not mol.natoms()==0:
                    mol=self.__class__()
                    mol.clear_atoms()
                    # additional optional fields
                    mol.setup_pwscf=mol.SETUP_PWSCF()
            #
            # Do READ IN
            #
            # read vectors after CELL_PARAMETERS
            if  opt=="readvec":
                # delete old vector
                if cntvec==0: vec=[]
                # read vector
                vec.append([float(linesplit[0]),
                            float(linesplit[1]),
                            float(linesplit[2])])
                cntvec+=1
                # set option and cntvec to zero after vectors are read
                if cntvec==3: 
                    opt=""
                    cntvec=0
            # read vectors after crystal axes
            if  opt=="readcrys":
                # delete old vector
                if cntvec==0: vec=[]
                # read vector
                vec.append([float(linesplit[3]),
                            float(linesplit[4]),
                            float(linesplit[5])])
                cntvec+=1
                # set option and cntvec to zero after vectors are read
                if cntvec==3: 
                    opt=""
                    cntvec=0
            # read coordinates
            if  opt=="readspecies":
                name=linesplit[0]
                mass=float(linesplit[2])
                pot=name+".uspp736.pbe.UPF"
                mol.set_element(name,mass,pot)
                cnttypes+=1
                if cnttypes==ntypes: opt=""
            # read unitvector
            elif opt=="readcoord":
                # append atom
                mol.append_atom(
                    self.__class__.atom(
                        mol,
                        cntat,
                        linesplit[0],-1,
                        float(linesplit[1]),
                        float(linesplit[2]),
                        float(linesplit[3])
                        )
                    )
                # cntat
                cntat+=1
                if int(cntat)==int(natoms): 
                    opt=""
                    cntat=0                   
                    # set celldm
                    mol.set_celldm(celldm*calc.b2A)    # in A
                    mol.setup_pwscf.set_celldm(celldm) # in bohr
                    # set periodicity
                    mol.set_vecs(vec[0],vec[1],vec[2])
                    # set real coordinates
                    mol.rel2real()
                    # set molecule and append
                    mol.set(filename,1,"")
                    # append mol to molecules
                    molecules.append(copy.copy(mol))
            # read input options
            elif opt=="":
                if   len(linesplit)>1 and linesplit[0:3]==["number","of","atoms/cell"]:
                    natoms=int(linesplit[4])
                elif len(linesplit)>3 and linesplit[0:4]==["number","of","atomic","types"]:
                    ntypes=int(linesplit[5])                        
                elif len(linesplit)>1 and linesplit[0]=="celldm(1)=":
                    celldm=(float(linesplit[1]))
            #
            # read main options
            #
            if len(linesplit)>0:
                option=linesplit[0]
                if   option=="atomic" and linesplit[1]=="species":
                    opt="readspecies"
                elif option=="ATOMIC_POSITIONS":
                    opt="readcoord"
                elif option=="CELL_PARAMETERS":
                    opt="readvec"
                elif option=="crystal" and linesplit[1]=="axes:":
                    opt="readcrys"

        # close file
        file.close()
        # return molecules
        return molecules

    # read input file options
    def readpwscfin_opt(self,linesplit,opt):
        if   opt=="system":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.system.append(linesplit[i].split("="))
        elif opt=="control":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.control.append(linesplit[i].split("="))
        elif opt=="electrons":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.electrons.append(linesplit[i].split("="))
        elif opt=="ions":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.ions.append(linesplit[i].split("="))
        elif opt=="readkpoints":
            line=" ".join(linesplit)
            self.setup_pwscf.kpoints=line
            opt=""
        return opt

    # find option blocks in pwscf
    # return opt
    def readpwscfin_blocks(self,linesplit,opt):
        # find keywords for blocks
        if len(linesplit)>0:
            option=linesplit[0]
            if  option=="&control":
                opt="control"
            elif option=="&system":
                opt="system"
            elif option=="&electrons":
                opt="electrons"
            elif option=="&ions":
                opt="ions"
            elif option=="ATOMIC_SPECIES":
                opt="readspecies"
            elif option=="ATOMIC_POSITIONS":
                opt="readcoord"
            elif option=="CELL_PARAMETERS":
                opt="readvec"
            elif option=="K_POINTS":
                opt="readkpoints"                    
        # return internal keywords
        return opt
        


    def rel2real(self):
        for atom in self.at():
            coo=[float(0.0),float(0.0),float(0.0)]
            # calculate coordinates
            for dim in range(len(atom.coord())):
                coo=calc.vecadd(coo,calc.scal_vecmult(
                    atom.coord()[dim],self.vec()[dim])
                )
            # set relative and real coordinates
            atom.coord_rel=atom.coord()
            atom.set_pos(coo)

    def read_setup_pwscf(self,filename):
        file=open(filename,"r")
        opt=""
        cntline=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            # read input options
            opt=self.readpwscfin_opt(linesplit,opt)
            # read main options
            opt=self.readpwscfin_blocks(linesplit,opt)
        file.close()
    #
    # PWSCF input file option class
    #
    class SETUP_PWSCF():
        def __init__(self):
            self.control=[]
            self.system=[]
            self.electrons=[]
            self.ions=[]
            # kpoints Na Nb Nc Sa Sb Sc
            self.kpoints="1 1 1 0 0 0"
            
        #############################################################
        # return functions
        #############################################################                      
        def Rcelldm(self):
            celldm=1.0
            for i in range(len(system)):
                if self.system[0]=="celldm(1)":
                    celldm=float(self.system[1])
            return celldm
        #############################################################
        # modify functions
        #############################################################    
        def set_celldm(self,celldm):
            # IN BOHR
            set=False
            # find celldm in system
            for i in range(len(self.system)):
                if self.system[i][0]=="celldm(1)": 
                    self.system[i][1]=str(celldm)
                    set=True
            # if it is not in there add it
            if set==False:
                self.system.append(["celldm(1)",str(celldm)])
            return
        
        def set_nat(self,nat):
            set=False
            # find nat in system
            for i in range(len(self.system)):
                if self.system[i][0]=="nat": 
                    self.system[i][1]=str(nat)
                    set=True
            # if it is not in there add it
            if set==False:
                self.system.append(["nat",str(nat)])
            return

        def set_ntyp(self,ntypes):
            set=False
            # find ntypes in system
            for i in range(len(self.system)):
                if self.system[i][0]=="ntyp": 
                    self.system[i][1]=str(ntypes)
                    set=True
            # if it is not in there add it
            if set==False:
                self.system.append(["ntyp",str(ntypes)])
            return
        

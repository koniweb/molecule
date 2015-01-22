#!/usr/bin/env python
##########################################################################
# writepwscf
# readfpwscfin
# readfpwscfout
##########################################################################
version=1.2
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
            atom.coord_crystal=[float(res[0]),float(res[1]),float(res[2])]
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
        ## CONTROL
        pstr="    "
        if hasattr(self.setup_pwscf, 'control'):
            if len(self.setup_pwscf.control)>0:
                print >>f, " {:s}".format("&control")
                for i in range(len(self.setup_pwscf.control)):
                    if self.setup_pwscf.control[i]=="LB":
                        print >>f , pstr[0:-2]
                        pstr="    "
                    else: 
                        pstr=('{:s}{:s}={:s}, ').format(
                            pstr,
                            self.setup_pwscf.control[i][0],
                            self.setup_pwscf.control[i][1])
                if self.setup_pwscf.control[-1]!="LB": print >>f , pstr[0:-2]
                print >>f, "/"
        ## SYSTEM
        pstr="    "
        if hasattr(self.setup_pwscf, 'system'):
            if len(self.setup_pwscf.system)>0:
                print >>f, " {:s}".format("&system")
                for i in range(len(self.setup_pwscf.system)):
                    if self.setup_pwscf.system[i]=="LB":
                        print >>f , pstr[0:-2]
                        pstr="    "
                    else: 
                        pstr=('{:s}{:s}={:s}, ').format(
                            pstr,
                            self.setup_pwscf.system[i][0],
                            self.setup_pwscf.system[i][1])
                print >>f, "/"
                if self.setup_pwscf.system[-1]!="LB": print >>f , pstr[0:-2]
        ## ELECTRONS
        pstr="    "
        if hasattr(self.setup_pwscf, 'electrons'):
            if len(self.setup_pwscf.electrons)>0:
                print >>f, " {:s}".format("&electrons")
                for i in range(len(self.setup_pwscf.electrons)):
                    if self.setup_pwscf.electrons[i]=="LB":
                        print >>f , pstr[0:-2]
                        pstr="    "
                    else: 
                        pstr=('{:s}{:s}={:s}, ').format(
                            pstr,
                            self.setup_pwscf.electrons[i][0],
                            self.setup_pwscf.electrons[i][1])
                print >>f, "/"
                if self.setup_pwscf.electrons[-1]!="LB": print >>f , pstr[0:-2]
        ## IONS
        pstr="    "
        if hasattr(self.setup_pwscf, 'ions'):
            if len(self.setup_pwscf.ions)>0:
                print >>f, " {:s}".format("&ions")
                for i in range(len(self.setup_pwscf.ions)):
                    if self.setup_pwscf.ions[i]=="LB":
                        print >>f , pstr[0:-2]
                        pstr="    "
                    else: 
                        pstr=('{:s}{:s}={:s}, ').format(
                            pstr,
                            self.setup_pwscf.ions[i][0],
                            self.setup_pwscf.ions[i][1])
                print >>f, "/"
                if self.setup_pwscf.ions[-1]!="LB": print >>f , pstr[0:-2]
        ## CELL
        if hasattr(self.setup_pwscf, 'cell'):
            print >>f, "&CELL"
            for i in range(len(self.setup_pwscf.cell)):
                print >>f, ('{:s}={:s}').format(
                    self.setup_pwscf.cell[i][0],
                    self.setup_pwscf.cell[i][1])
            print >>f, "/"
        ## ATOMIC SPECIES
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
        ## ATOMIC POSITIONS in crystal
        print >>f, "ATOMIC_POSITIONS  {:s}".format(self.setup_pwscf.atomic_positions_info)
        for atom in self.at():
            # define fixes for atoms
            fixatom="  "
            for i in range(len(atom.fixes())): fixatom+=" {:s}".format(str(atom.fixes()[i]))
            # print atoms
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f}{:s}'.format(
                    atom.type()[0],
                    atom.coord_crystal[0], 
                    atom.coord_crystal[1], 
                    atom.coord_crystal[2],
                    fixatom
                    )
                )
        print >>f
        ## CELL PARAMETERS
        print >>f, "CELL_PARAMETERS  {:s}".format(self.setup_pwscf.cell_parameters_info)
        for cntvec in range(3):
            print >>f , (
                '{:15.10f} {:15.10f} {:15.10f}'.format(
                    self.celldm_vec()[1][cntvec][0], 
                    self.celldm_vec()[1][cntvec][1], 
                    self.celldm_vec()[1][cntvec][2]                    
                    )
                )
        ## KPOINTS
        if (hasattr(self.setup_pwscf,"kpoints")):
            k=self.setup_pwscf.kpoints
            print >>f
            print >>f, "K_POINTS  {:s}".format(self.setup_pwscf.kpoints_info)
            # automatic
            if self.setup_pwscf.kpoints_info=="automatic":
                print >>f , ('   {:d} {:d} {:d}   {:d} {:d} {:d}'.format(k[0][0],k[0][1],k[0][2],
                                                                         k[0][3],k[0][4],k[0][5]))
            # rest without gamma
            elif self.setup_pwscf.kpoints_info!="gamma":
                print >> f, ("{:d}").format(k[0][0])
                for i in range(1,len(k)):
                    print >>f , ('{:f} {:f} {:f}   {:f}'.format(k[i][0],k[i][1],k[i][2],
                                                                   k[i][3]))
        # close filename
        if filename != "": f.close()
        return
    
    # read molecules in pwscf input file format
    def readpwscfin(self,filename):
        # set molecule
        molecules=[]
        natoms=0
        ntypes=0
        vec=[]
        # check read file
        try: 
            file=open(filename, 'r')
        except IOError:
            print >> sys.stderr, "... input file not found"
            exit()
        # read file
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
                    for i in range(4,len(linesplit)): linesplit[i]=int(linesplit[i])
                    mol.at()[-1].set_fixes(linesplit[4:])
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
        mol.crystal2angstroem()
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
        cntmol=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            #
            # create new molecule
            #
            # first vector readin or standart vector readin
            if cntline==1:
                mol=self.__class__()
                mol.clear_atoms()
                # additional optional fields
                mol.setup_pwscf=mol.SETUP_PWSCF()
            if opt=="endmol":
                opt=""
                if cntline==1 or not mol.natoms()==0:
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
            # read vectors after crystal axes (0. step)
            if  opt=="readvec0":
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
            # read species
            if  opt=="readspecies":
                name=linesplit[0]
                mass=float(linesplit[2])
                pot=name+".uspp736.pbe.UPF"
                mol.set_element(name,mass,pot)
                cnttypes+=1
                if cnttypes==ntypes: 
                    opt=""
                    cnttypes=0
            # read coordinates of atoms ATOMIC_POSITIONS
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
                    cntmol+=1
            # read coordinates of atoms "atom positions" (0. step)
            elif opt=="readcoord0":
                # append atom
                mol.append_atom(
                    self.__class__.atom(
                        mol,
                        cntat,
                        linesplit[1],-1,
                        float(linesplit[6]),
                        float(linesplit[7]),
                        float(linesplit[8])
                        )
                    )
                # cntat
                cntat+=1
                if int(cntat)==int(natoms): 
                    opt=""
                    cntat=0                   
            # read energy
            elif opt=="" and len(linesplit)>0 and linesplit[0]=="!":
                # read energy and set it
                E=float(linesplit[4])
                mol.set_energy(E)
                #
                # set all molecule features
                #
                # set celldm
                mol.set_celldm(celldm*calc.b2A)    # in A
                mol.setup_pwscf.set_celldm(celldm) # in bohr
                # set periodicity
                mol.set_vecs(vec[0],vec[1],vec[2])
                # set real coordinates
                if cntmol==0: mol.alat2angstroem()
                else:         mol.crystal2angstroem()

                # set molecule and append
                mol.set(filename,1,"")
                # append mol to molecules
                molecules.append(copy.copy(mol))
                #
                # print info that molecule is finished
                #
                opt="endmol"
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
                if  (len(linesplit)>3 and 
                     option=="atomic" and linesplit[1]=="species" and linesplit[2]=="valence"):
                    opt="readspecies"
                elif option=="ATOMIC_POSITIONS":
                    opt="readcoord"
                elif len(linesplit)>3 and linesplit[2]=="atom" and linesplit[3]=="positions":
                    opt="readcoord0"
                elif option=="CELL_PARAMETERS":
                    opt="readvec"
                elif option=="crystal" and linesplit[1]=="axes:":
                    opt="readvec0"
            #
            # if Final enthalpy/ energy keyword --> exit loop # RECHECK
            #
            if len(linesplit)>0 and linesplit[0]=="Final":
                break
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
                self.setup_pwscf.system.append("LB") #TEST
        elif opt=="control":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.control.append(linesplit[i].split("="))
                self.setup_pwscf.control.append("LB") #TEST
        elif opt=="electrons":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.electrons.append(linesplit[i].split("="))
                self.setup_pwscf.electrons.append("LB") #TEST
        elif opt=="ions":
            if linesplit[0]=="/": opt=""
            else:
                for i in range(len(linesplit)):
                    linesplit[i]=linesplit[i].replace(",","")
                    self.setup_pwscf.ions.append(linesplit[i].split("="))
                self.setup_pwscf.ions.append("LB") #TEST
        elif opt=="readkpoints":
            # read info for automatic -- only one line
            if self.setup_pwscf.kpoints_info=="automatic":
                for i in range(len(linesplit)): linesplit[i]=int(linesplit[i])
                self.setup_pwscf.kpoints.append(linesplit)
                opt=""
            # read info for automatic -- no line
            elif self.setup_pwscf.kpoints_info=="gamma":
                opt=""
            # read info for multiple kpoint lines with positions
            else:
                print "yeah"
                if len(linesplit)==1: 
                    self.setup_pwscf.kpoints.append([int(linesplit[0]),0])
                else:
                    self.setup_pwscf.kpoints[0][1]+=1
                    if self.setup_pwscf.kpoints[0][1]<=self.setup_pwscf.kpoints[0][0]:
                        for i in range(len(linesplit)): linesplit[i]=float(linesplit[i])
                        self.setup_pwscf.kpoints.append(linesplit)
                    else: opt=""
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
                # additional argument for atomic positions
                if len(linesplit)==2 : self.setup_pwscf.atomic_positions_info=linesplit[1]
                else:                  self.setup_pwscf.atomic_positions_info=""                
            elif option=="CELL_PARAMETERS":
                opt="readvec"
                if len(linesplit)==2 : self.setup_pwscf.cell_parameters_info=linesplit[1]
                else:                  self.setup_pwscf.cell_parameters_info=""
            elif option=="K_POINTS":
                opt="readkpoints"
                self.setup_pwscf.kpoints=[]
                if len(linesplit)==2 : self.setup_pwscf.kpoints_info=linesplit[1]
                else:                  self.setup_pwscf.kpoints_info=""
        # return internal keywords
        return opt
        
    def crystal2angstroem(self):
        for atom in self.at():
            coo=[float(0.0),float(0.0),float(0.0)]
            # calculate coordinates
            for dim in range(len(atom.coord())):
                coo=calc.vecadd(coo,calc.scal_vecmult(
                    atom.coord()[dim],self.vec()[dim])
                )
            # set relative and real coordinates
            atom.coord_crystal=atom.coord()
            atom.set_pos(coo)

    def alat2angstroem(self):
        for atom in self.at():
            coo=[float(0.0),float(0.0),float(0.0)]
            # calculate coordinates
            coo=calc.scal_vecmult(self.celldm_vec()[0],atom.coord())
            # set relative and real coordinates
            atom.coord_crystal=atom.coord()
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
            self.kpoints=[[1, 1, 1, 0, 0, 0]]
            self.kpoints_info="automatic"
            # cell parameters
            self.cell_parameters_info=""
            # atomic positions
            self.atomic_positions_info="crystal"
            self.atomic_position_fixes=[]

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
        

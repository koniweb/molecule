#!/usr/bin/env python
##########################################################################
# Molecule Class
# subclasses:
#   Atom Class
#   Bond Class
##########################################################################
version=3.6
versiontext='# class_molecule.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import copy
import re
from operator import itemgetter, attrgetter
import mod_calc as calc # several functions
import numpy
#from numpy import matrix
#from numpy import linalg

ndim=calc.ndim
pse=[
    # pse[element][0] name
    # pse[element][1] type
    # pse[element][2] ~ weight
    # pse[element][3] pseudopotential
        ["X",  0, 0.0,"pseudopotential"],
        # first row
        ["H",  1, 1.0079, "PSEUDO"], ["He",2,4.0026,   "PSEUDO"],

        # second row
        ["Li", 3, 6.941,  "PSEUDO"], ["Be", 4, 9.0122, "PSEUDO"], 
        ["B" , 5,10.811,  "PSEUDO"], ["C" , 6,12.0107, "PSEUDO"], ["N" , 7,14.007, "PSEUDO"], 
        ["O" , 8,15.999,  "PSEUDO"], ["F" , 9,18.998,  "PSEUDO"], ["Ne",10,20.18,  "PSEUDO"],
        
        # third row
        ["Na",11, 22.99,   "PSEUDO"], ["Mg",12, 24.305,  "PSEUDO"],  
        ["Al",13, 26.982,  "PSEUDO"], ["Si",14, 28.086,  "PSEUDO"], ["P" ,15, 30.974, "PSEUDO"], 
        ["S", 16, 32.065,  "PSEUDO"], ["Cl",17, 35.453,  "PSEUDO"], ["Ar",18, 39.948, "PSEUDO"],
        
        # fourth row
        ["K", 19, 39.098,  "PSEUDO"], ["Ca",20, 40.078,  "PSEUDO"], 
        # transition metals
        ["Sc",21, 44.956,  "PSEUDO"], ["Ti",22, 47.867,  "PSEUDO"], ["V" ,23, 50.942, "PSEUDO"], 
        ["Cr",24, 51.996,  "PSEUDO"], ["Mn",25, 54.938,  "PSEUDO"], ["Fe",26, 55.845, "PSEUDO"], 
        ["Co",27, 58.933,  "PSEUDO"], ["Ni",28, 58.693,  "PSEUDO"], ["Cu",29, 63.546, "PSEUDO"], 
        ["Zn",30, 65.380,  "PSEUDO"], 
        # transition metals end
        ["Ga",31, 69.723,  "PSEUDO"], ["Ge",32, 72.640,  "PSEUDO"], ["As",33, 74.922, "PSEUDO"], 
        ["Se",34, 78.96,   "PSEUDO"], ["Br",35, 79.904,  "PSEUDO"], ["Kr",36, 83.798, "PSEUDO"],

        # fifth row
        ["Rb",37, 85.468,  "PSEUDO"], ["Sr",38, 87.620,  "PSEUDO"], 
        # transition metals
        ["Y" ,39, 88.906,  "PSEUDO"], ["Zr",40, 91.224,  "PSEUDO"], ["Nb",41, 92.906, "PSEUDO"], 
        ["Mo",42, 95.960,  "PSEUDO"], ["Tc",43, 97.900,  "PSEUDO"], ["Ru",44,101.070, "PSEUDO"], 
        ["Rh",45,102.910,  "PSEUDO"], ["Pd",46,106.420,  "PSEUDO"], ["Ag",47,107.870, "PSEUDO"], 
        ["Cd",48,112.410,  "PSEUDO"], 
        # transition metals end
        ["In",49,114.820,  "PSEUDO"], ["Sn",50,118.710,  "PSEUDO"], ["Sb",51,121.760, "PSEUDO"], 
        ["Te",52,127.600,  "PSEUDO"], ["I" ,53,126.900,  "PSEUDO"], ["Xe",54,131.290, "PSEUDO"]        
        ]
#----------------------------------------------------------------------
# classes
#----------------------------------------------------------------------
import mod_xyz    as mxyz  # molecule subclass
import mod_pwscf  as mpw   # molecule subclass
import mod_lmp    as mlmp  # molecule subclass
import mod_extend as mext  # molecule subclass

######################################################################
# MOLECULE CLASS
######################################################################
class molecule(mxyz.molecule_rw,mpw.molecule_rw,mlmp.molecule_rw,
               mext.molecule_extend):
    # initialize
    def __init__(self):           
        # pse
        self.__pse=copy.deepcopy(pse)
        # set molecule info
        self.__file=""
        self.__filemolnumber=0
        self.__comment=""
        self.__typelist=[] # which element by name
        self.__id=id(self)
        # set atom list
        self.__at=[]
        self.__celldm=1.0
        self.__vec=[[0.0 for x in xrange(0,ndim)]for x in xrange(0,ndim)]
        # additional data
        self.__data=[] # [[ info, [data atom0,data atom1,...]],...]
        # vec[0]: a
        # vec[1]: b
        # vec[2]: c
        self.__offset=[0.0,0.0,0.0]
        # energy
        self.__energy=0.0

    #############################################################
    # return functions
    #############################################################
    # return id
    def id(self):
        return self.__id

    # return at
    def at(self):
        return self.__at

    # return periodicity vectors
    def vec(self):
        return (calc.scal_matmult(self.__celldm,self.__vec))
    def celldm_vec(self):
        return [self.__celldm,self.__vec]
        
    # return offset vector
    def offset(self):
        return self.__offset
        
    # return natoms
    def natoms(self):
        return len(self.__at)
        
    # return typelist
    def typelist(self):
        return self.__typelist

    # return ntypes
    def ntypes(self):
        return len(self.typelist())
        
    # return comment
    def comment(self):
        return self.__comment

    # return file
    def file(self):
        return self.__file

    # return filemolnumber
    def filemolnumber(self):
        return self.__filemolnumber

    # return pse
    def pse(self):
        return self.__pse

    # converts number, name or weight of atom to element
    def name2element(self,name):
        localname=re.sub("[^a-zA-Z]","", name)
        returndata=[]
        for ielement in self.pse():
            if  localname==ielement[0]: returndata=ielement
        # if not found remove last character
        if len(returndata)==0:
            for ielement in self.pse():
                if  localname[0:-1]==ielement[0]: returndata=ielement
        return returndata
    def number2element(self,number):
        returndata=[]
        for ielement in self.pse():
            if number==ielement[1]: returndata=ielement
        return returndata  
    def weight2element(self,weight):
        returndata=[]
        weightdiff=0.01
        for ielement in self.pse():
            if ( float(weight)<=float(ielement[2])+weightdiff and
                 float(weight)>=float(ielement[2])-weightdiff ): 
                returndata=ielement
        return returndata

    # return energy
    def energy(self):
        return self.__energy

    # return array bondlenghts
    def bondlengths(self,types=[]):
        list=[]
        for at in self.at():
            list.append(at.bondlengths(types))
        return list
    
    # def return additional data
    def data(self):
        return self.__data

    #############################################################
    # set functions
    #############################################################    
    # set id
    def set_id(self,id):
        self.__id=id
        
    # append to typelist
    def typelist_append(self,tid,name=""):
        self.__typelist.append([int(tid),name])
        return

    def typelist_append_check(self,name="",number=-1):
        type=[]
        if bool(type)==False: type=self.name2element(name)
        if bool(type)==False: type=self.name2element(re.sub("[a-zA-Z]","",name))
        if bool(type)==False: type=self.number2element(number)
        if bool(type):
            if [type[1],name] not in self.typelist():
                self.typelist_append(type[1],name)
                tid=self.ntypes()-1
            else: 
                tid=self.typelist().index([type[1],name])
        return tid

    # set typelist
    def set_typelist2list(self,list):
        self.__typelist=list
        return

    # set atomlist
    def set_atomlist(self,list):
        self.__at=list
        return
        
    # clear atom list
    def clear_atoms(self):
        self.__at=[]
        return

    # delete atoms
    def delete_atoms(self,idlist=[]):
        # sort in reverse order
        idlist.sort(reverse=True)
        # delete from last to first
        for item in idlist:
            del self.__at[item]
        return

    # append atom
    def append_atom_coo(self,type,x,y,z,pos=-1):
        # if position is -1 append atom at the end
        if pos==-1:pos=self.natoms()
        # append atom
        i=len(self.at())
        if type.isdigit(): self.__at.insert(pos, self.atom(
                self,i,""  ,int(type),
                float(x),float(y),float(z)))
        else:              self.__at.insert(pos,self.atom(
                self,i,type,-1       ,
                float(x),float(y),float(z)))
        return

    # append an instance atom
    def append_atom_cp(self,addat,pos=-1):
        # if position is -1 append atom
        if pos==-1:pos=self.natoms()
        # insert atom at position
        self.__at.insert(pos,self.atom(addat.mol(),
                                   addat.id(),
                                   addat.type()[0],
                                   addat.type()[1],
                                   addat.coord()[0],
                                   addat.coord()[1],
                                   addat.coord()[2],
                                   addat.mult()[0],
                                   addat.mult()[1],
                                   addat.mult()[2],
                                   atomcharge=addat.charge())
                       )
        return

    # append atom
    def append_atom(self,atom,pos=-1):
        # if position is -1 append atom
        if pos==-1:pos=self.natoms()
        # insert atom at position
        self.__at.insert(pos,atom)
        return

    # set type in pse
    def set_element(self,identification,mass=-1.0,pseudopotential=""):
        # if name calculate id
        if type(identification)==type(''): 
            tid=self.name2element(str(identification))[1]
        else: tid=int(identification)
        # set pseudopotential
        if not pseudopotential=="":
            self.__pse[tid][3]=pseudopotential
        # correct mass
        if mass>0.0:
            self.__pse[tid][2]=float(mass)
        return

    # set periodicity vectors separately
    def set_vecs(self,a=[0.0,0.0,0.0],b=[0.0,0.0,0.0],
                 c=[0.0,0.0,0.0],off=[0.0,0.0,0.0]):
        if not (a[0]==0.0 and a[1]==0.0 and a[2]==0.0):
            self.__vec[0]=a
        if not (b[0]==0.0 and b[1]==0.0 and b[2]==0.0):
            self.__vec[1]=b
        if not (c[0]==0.0 and c[1]==0.0 and c[2]==0.0):
            self.__vec[2]=c
        if not (off[0]==0.0 and off[1]==0.0 and off[2]==0.0):
            self.__offset=off
        return

    # set celldm
    def set_celldm(self,length):
        self.__celldm=float(length)

    # set comment
    def set_comment(self,comment):
        self.__comment=comment
        return    

    # set natoms ntypes and filename
    def set(self,filename="",filemolnumber=0,comment=""):
        self.__file=filename
        self.__filemolnumber=filemolnumber
        self.set_comment(comment)
        return

    # set energy
    def set_energy(self,energy):
        self.__energy=float(energy)
        return

    # def set additional data
    def set_data(self,data):
        self.__data=data
        return 

    #############################################################
    # modify functions
    #############################################################    
    # append molecule
    def append_mol(self,molecule,shiftv=[0.0,0.0,0.0],
                rotangle=0.0,rota=[1.0,0.0,0.0],rotp=[0.0,0.0,0.0]):
         # rotate and shift molecules
        molecule.rot(rotangle,
                     rota[0],rota[1],rota[2],
                     rotp[0],rotp[1],rotp[2])
        molecule.shift(shiftv[0],shiftv[1],shiftv[2])
        # add to new molecule
        for iat in range(0,molecule.natoms()):
            self.append_atom_cp(molecule.at()[iat])
        # shift and rotate molecules back
        molecule.shift(-shiftv[0],-shiftv[1],-shiftv[2])
        molecule.rot(-rotangle,rota[0],rota[1],rota[2],rotp[0],rotp[1],rotp[2])
        return

    # set number range
    def set_number(self,type,st=0,end=0):
        if end <= 0: 
            end = self.natoms() + end
            if end <=0: end = self.natoms()
        print >> sys.stderr, '... setting atom {:d} - {:d} type {:d}'.format(st,end,type)
        for cntat in range(st,end):
            self.at()[cntat].set_number(int(type))
        return

    # set name range
    def set_name(self,name,st=0,end=0):
        if end <= 0: 
            end = self.natoms() + end
            if end <=0: end = self.natoms()
        print >> sys.stderr, '... setting atom {:d} - {:d} name {:s}'.format(st,end,name)
        for cntat in range(st,end):
            self.at()[cntat].set_name(name)
        return
    
    # shift molecule
    def shift(self,x,y,z):
        for cnt in range(0,self.natoms()):
            self.at()[cnt].shift(x,y,z)
        return
            
    # rotate molecule
    def rot(self,angle,ax,ay,az,px,py,pz):
        # abbreviations and angles
        cos=math.cos
        sin=math.sin
        angle=angle*2.*math.pi/360.0
        # move molecule prior to rotation
        self.shift(-px,-py,-pz)
        # normalize vector
        la=math.sqrt(ax*ax + ay*ay + az*az)
        ax=ax/la
        ay=ay/la
        az=az/la
        # rotation matrix
        mat=[[0.0,0.0,0.0],
             [0.0,0.0,0.0],
             [0.0,0.0,0.0]]
        # first column
        mat[0][0]=ax*ax*(1.-cos(angle))+   cos(angle)
        mat[1][0]=ax*ay*(1.-cos(angle))-az*sin(angle)
        mat[2][0]=ax*az*(1.-cos(angle))+ay*sin(angle)
        # second column
        mat[0][1]=ay*ax*(1.-cos(angle))+az*sin(angle)
        mat[1][1]=ay*ay*(1.-cos(angle))+   cos(angle)
        mat[2][1]=ay*az*(1.-cos(angle))-ax*sin(angle)
        # third column
        mat[0][2]=az*ax*(1.-cos(angle))-ay*sin(angle)
        mat[1][2]=az*ay*(1.-cos(angle))+ax*sin(angle)
        mat[2][2]=az*az*(1.-cos(angle))+   cos(angle)
        # do rotation on all atoms
        for atom in self.at():
            coo=atom.coord()
            x= mat[0][0]*coo[0] + mat[1][0]*coo[1] + mat[2][0]*coo[2]
            y= mat[0][1]*coo[0] + mat[1][1]*coo[1] + mat[2][1]*coo[2]
            z= mat[0][2]*coo[0] + mat[1][2]*coo[1] + mat[2][2]*coo[2]
            coo[0]=x
            coo[1]=y
            coo[2]=z
        # rotate vectors
        vecs=self.celldm_vec()[1]
        vecx=[ (mat[0][0]*vecs[0][0]+mat[1][0]*vecs[0][1]+mat[2][0]*vecs[0][2]),
               (mat[0][1]*vecs[0][0]+mat[1][1]*vecs[0][1]+mat[2][1]*vecs[0][2]),
               (mat[0][2]*vecs[0][0]+mat[1][2]*vecs[0][1]+mat[2][2]*vecs[0][2])]
        vecy=[ (mat[0][0]*vecs[1][0]+mat[1][0]*vecs[1][1]+mat[2][0]*vecs[1][2]),
               (mat[0][1]*vecs[1][0]+mat[1][1]*vecs[1][1]+mat[2][1]*vecs[1][2]),
               (mat[0][2]*vecs[1][0]+mat[1][2]*vecs[1][1]+mat[2][2]*vecs[1][2])]
        vecz=[ (mat[0][0]*vecs[2][0]+mat[1][0]*vecs[2][1]+mat[2][0]*vecs[2][2]),
               (mat[0][1]*vecs[2][0]+mat[1][1]*vecs[2][1]+mat[2][1]*vecs[2][2]),
               (mat[0][2]*vecs[2][0]+mat[1][2]*vecs[2][1]+mat[2][2]*vecs[2][2])]
        for i in range(ndim):
            if vecx[i]<10**-10 and vecx[i]>-10**-10: vecx[i]=0.0
            if vecy[i]<10**-10 and vecy[i]>-10**-10: vecy[i]=0.0
            if vecz[i]<10**-10 and vecz[i]>-10**-10: vecz[i]=0.0
        self.set_vecs(vecx,vecy,vecz)
        # move molecule back after rotation
        self.shift(px,py,pz)

    # stretch structure
    def stretch(self,factor):
        # stretch atoms
        for atom in self.at():
            atom.set_pos( calc.vec_vecmult(factor,atom.coord()) )
        # stretch vectors
        if not (factor[0]==factor[1] and factor[1]==factor[2]):
            self.set_vecs(calc.vec_vecmult(factor,self.celldm_vec()[1][0]),
                          calc.vec_vecmult(factor,self.celldm_vec()[1][1]),
                          calc.vec_vecmult(factor,self.celldm_vec()[1][2]))
        else: self.set_celldm(self.celldm_vec()[0]*factor[0])
        # stretch offset
        self.set_vecs(off=calc.vec_vecmult(factor,self.offset()) )
        return

    # wrap atoms into unitcell
    def wrapmol(self):
        if self.vec()==[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]:
            print >> sys.stderr, "... wrapping not possible"
            print >> sys.stderr, "... vectors have to be set"
        else:
            nwrap=self.func_on_alat(self.wrap_atom)
            print >> sys.stderr, "... {:d} atoms in molecule wrapped".format(nwrap)
        return

    # do wrap on coordinate of atoms
    def wrap_atom(self,coord_alat,atom,nwrap):
        # calculate shift to box
        shift=[0.0,0.0,0.0]
        outofbox=True
        while outofbox==True:
            for i in range(ndim):
                if coord_alat[i]>1.0: 
                    for idim in range(ndim): shift[idim]-=self.vec()[i][idim]
                    coord_alat[i]-=1
                elif coord_alat[i]<0.0: 
                    for idim in range(ndim): shift[idim]+=self.vec()[i][idim]
                    coord_alat[i]+=1
            outofbox=False
            for i in range(ndim):
                if coord_alat[i]>1.0 or coord_alat[i]<0.0:  outofbox=True
                
        # reset coordinates
        if shift!=[0.0,0.0,0.0]:
            nwrap+=1
            atom.set_pos(calc.vecadd(atom.coord(),shift))
        return nwrap

    # calculate relative coordinates and operate function on those
    def func_on_alat(self,function):
        counter=0
        vecM = numpy.matrix( [ [self.vec()[0][0],self.vec()[1][0],self.vec()[2][0]],
                         [self.vec()[0][1],self.vec()[1][1],self.vec()[2][1]],
                         [self.vec()[0][2],self.vec()[1][2],self.vec()[2][2]] ])
        for atom in self.at():
            tmp = numpy.matrix([ [atom.coord()[0]-self.offset()[0]],
                           [atom.coord()[1]-self.offset()[1]],
                           [atom.coord()[2]-self.offset()[2]] ])
            coord_alat=numpy.linalg.solve(vecM, tmp)
            counter=function(coord_alat,atom,counter)
        return counter
    
    # rotate vectors a->b b->c c->a
    def rotate_vecs(self):
        vecs=self.celldm_vec()[1]
        self.set_vecs(vecs[1],vecs[2],vecs[0])
        return

    # exchange vectors a->b b->a
    def exchange_vecs(self):
        vecs=self.celldm_vec()[1]
        self.set_vecs(vecs[1],vecs[0],vecs[2])
        return
    
    # set vectors part to all plus
    def repositiv_vecs(self):
        vecs=self.celldm_vec()[1]
        for i in range(ndim):
            for j in range(ndim):
                if vecs[i][j]<0.0: vecs[i][j]=-vecs[i][j]
        self.set_vecs(vecs[0],vecs[1],vecs[2])
        return

    # sorting atoms via x,y,z
    def sortatoms(self,dir,factor=1):
        if dir>=0 and dir<=2:
            self.set_atomlist(sorted(self.at(), key=lambda atom:factor*atom.coord()[dir]))
        return


    # find bonds for atoms with length smaller than cutoff
    def define_bonds(self,cutoff,cutmin=10**(-10),atomrange=[0,-1],atomlist=[],
                     periodicity=False,nbondmax=10):
        # create atomlist
        if len(atomlist)==0:
            if atomrange[1]==-1: atomrange[1]=self.natoms()
            atomlist=range(atomrange[0],atomrange[1])
        # fortran bondcnt
        bndcnt=self.Fdefine_bonds(cutoff,cutmin,atomlist,periodicity,nbondmax)
        # if there is an error try a larger nbondmax
        if bndcnt<0:
            bndcnt=self.Fdefine_bonds(cutoff,cutmin,atomlist,periodicity,nbondmax*50)
            # if it still does not work use the python function
            if bnccnt<0:
                print >> sys.stderr, "...use python to define bonds"
                bndcnt=self.Pdefine_bonds(cutoff,cutmin,atomlist,periodicity)
        return bndcnt

    # find bonds for atoms with length smaller than cutoff
    def Pdefine_bonds(self,cutoff,cutmin=10**(-10),atomlist=[],periodicity=False):
        bndcnt=0
        for at in self.at():
            # print if 
            if (at.id()%100)==0 or at.id()==self.natoms()-1: 
                print >> sys.stderr, "...bonding for atom {:d} of {:d} calculated".format(at.id()+1,self.natoms())
            # only check atoms in range
            if at.id() in atomlist: 
                # do bond calculation for periodic structures
                if periodicity:
                    perx=[-1,0,1]
                    pery=[-1,0,1]
                    perz=[-1,0,1]
                else:
                    perx=[0]
                    pery=[0]
                    perz=[0]
                for x in perx:
                    for y in pery:
                        for z in perz:
                            for neigh in self.at():             
                                l=calc.a_dist(at, neigh, per=[x,y,z])
                                if (l>cutmin and l<cutoff):
                                    #print at.id(),neigh.id(),x,y,z,l #infos
                                    at.add_bond(neigh,per=[x,y,z]) 
                                    bndcnt+=1
        return bndcnt

    def Fdefine_bonds(self,cutoff,cutmin=10**(-10),atomlist=[],periodicity=False,nbondmax=10):
        from fortran_modules import fortran_modules as f
        # only check atoms in range
        alist=numpy.array(copy.deepcopy(atomlist))
        # bonding -> change
        bonding=numpy.array([[int(-1) for x in xrange(5)]for x in xrange(self.natoms()*nbondmax)],order='F')
        error=False
        # call fortran code
        #print self.vec()[0] # DEBUG
        f.define_bonds(
            cutoff,
            cutmin,
            numpy.array([i.coord() for i in self.at()]),
            alist,
            periodicity,
            numpy.array(self.vec()),
            bonding,
            error
        )
        # add bonds
        bndcnt=0
        for bond in bonding:
            if bond[0]==-1: 
                break
            else:
                if bond[0]!=-1: self.at()[bond[0]].add_bond(self.at()[bond[1]],per=[bond[2],bond[3],bond[4]])
                bndcnt+=1
        # return
        if error:
            print >> sys.stderr, "test" #DEBUG
            return -1
        else:
            return bndcnt


######################################################################
# ATOM CLASS
######################################################################
    class atom:   
        # initialize
        def __init__(self,mol,id,name,number,
                     x,y,z,mx=0,my=0,mz=0,atomcharge=0.0,tid=-1):
            self.__parentmol=mol
            self.__id=id
            if name=="" and number>0:     name   = self.__parentmol.number2element(number)[0] 
            if number<0 and not name=="": number = self.__parentmol.name2element(name)[1] 
            self.__tid=self.__parentmol.typelist_append_check(name,number)
            # element is self.parentmol().pse()[self.parentmol().typelist[tid][0]]
            #   you can get elementinfos by using self.type() 
            #   !!! Attention number is still local number
            self.__coord=[x,y,z]
            self.__mult=[mx,my,mz]
            self.__charge=atomcharge
            # bonds
            self.__bonds=[]
            # fixes
            self.__fixes=[] # fixed in x y or z direction: 1 for moveable 0 for fixed

        #############################################################
        # return functions
        #############################################################
        def parentmol(self):
            return self.__parentmol

        def id(self):
            return self.__id

        def type(self):
            if self.parentmol().typelist()[self.tid()][1]=="": 
                name=self.parentmol().pse()[self.parentmol().typelist()[self.tid()][0]][0]
            else: name   = self.parentmol().typelist()[self.tid()][1]
            number = self.parentmol().typelist()[self.tid()][0]
            weight = self.parentmol().pse()[self.parentmol().typelist()[self.tid()][0]][2]
            pot    = self.parentmol().pse()[self.parentmol().typelist()[self.tid()][0]][3]
            return [name,number,weight,pot]

        def tid(self):
            return self.__tid
            
        def mult(self):
            return self.__mult

        def coord(self):
            return self.__coord

        def charge(self):
            return self.__charge
        
        def bonds(self):
            return self.__bonds

        def nbonds(self,cutoff=-1.0,types=[]):
            if cutoff<0.0: N=len(self.bonds())
            else:
                N=0
                for b in self.bonds(): 
                    # if types are defined
                    if len(types)>0: 
                        if ( (types[0]==b.atom().type()[0]    and 
                              types[1]==b.neighbor().type()[0])  or
                             (types[1]==b.atom().type()[0]    and 
                              types[0]==b.neighbor().type()[0])  ):
                            if b.bondlength()<cutoff:N+=1
                    # for all types
                    else:
                        if b.bondlength()<cutoff:N+=1                    
            return N
      
        def bondlengths(self,types=[]):
            bondlengths=[]
            for bond in self.bonds():
                # if types are defined
                if len(types)>0: 
                    if ( (types[0]==bond.atom().type()[0]    and 
                          types[1]==bond.neighbor().type()[0])  or
                         (types[1]==bond.atom().type()[0]    and 
                          types[0]==bond.neighbor().type()[0])  ):
                        bondlengths.append(bond.bondlength())
                # for all types
                else:
                    bondlengths.append(bond.bondlength())
            return bondlengths

        # return bondangles
        def bondangles(self,cutoff=[],types=[]):
            #types[0]   center of angle (atom)
            #types[1:2]  sides of angle (neighbor)
            bondangles=[]
            # check cutoffs -> too expensive??
            if True is not ((len(types)==3 and len(cutoff)==2) or len(types)==0):
                print >> sys.stderr, "...cutoffs not defined in bondangles function"
                exit()
            # calculate angle with AxB=|A|x|B|*cos(theta)
            for bond_i in self.bonds():
                for bond_j in self.bonds():
                    # if types are defined
                    if len(types)>0: 
                        if ( (types[0]==bond_i.atom().type()[0]     and 
                              types[1]==bond_i.neighbor().type()[0] and
                              types[2]==bond_j.neighbor().type()[0] ) ):
                            self.calc_bondangle_append(bondangles,bond_i,bond_j,cutoff)
                    # for all types
                    else:
                        self.calc_bondangle_append(bondangles,bond_i,bond_j,cutoff)
            return bondangles
        # calculation externalized
        def calc_bondangle_append(self,bondangles,bond1,bond2,cutoff):
            acc=10**-6 # accuracy for angle calculation
            # if both bonds are smaller than cutoff
            if bond1.bondlength()<cutoff[0] and bond2.bondlength()<cutoff[1]:
                preangle=( calc.vecprod(bond1.bondvector(),bond2.bondvector())
                           /calc.length(bond1.bondvector())
                           /calc.length(bond2.bondvector()) )
                # if angle is not 360 or 0
                if preangle<1-acc: 
                    bondangles.append(numpy.arccos(preangle)*180.0/numpy.pi)

        def fixes(self):
            return self.__fixes

        #############################################################
        # set functions
        #############################################################
        # set position of atom
        def set_pos(self,vec):
            self.__coord=vec
    
        # set typeid
        def set_tid(self,tid):
            self.__tid=tid

        # set molecule
        def set_parentmol(self,mol):
            self.__parentmol=mol

        # set charge of atom
        def set_charge(self,charge):
            self.__charge==charge

        # set bond of atom to neighbor atom in periodicity box
        def add_bond(self,neighbor,per=[0,0,0]):
            self.__bonds.append(self.parentmol().bond(self,neighbor,per))

        # set fixes for atoms
        def set_fixes(self,fixes):
            if len(fixes)<=3: self.__fixes=fixes           

        # set fixes for atoms
        def set_bonds(self,newbonds):
            self.__bonds=newbonds

        #############################################################
        # modify functions
        #############################################################   
        # shift atom
        def shift(self,x,y,z):
            self.set_pos([self.coord()[0]+x,self.coord()[1]+y,self.coord()[2]+z])

######################################################################
# BOND CLASS
######################################################################       
    class bond:
        def __init__(self,atom,neighbor,per=[0,0,0]):
            # define bond start and end
            self.__atom=atom
            self.__neighbor=neighbor
            # define periodicity
            self.__per=per

        #############################################################
        # return functions
        #############################################################
        def atom(self):
            return self.__atom

        def neighbor(self):
            return self.__neighbor
            
        def per(self):
            return self.__per
        
        def bondlength(self):
            return calc.a_dist(self.atom(),self.neighbor(),self.per())

        def bondvector(self):
            return calc.a_vec(self.atom(),self.neighbor(),self.per())
        #############################################################
        # set functions
        #############################################################
        def reset_bond(self,neighbor,per=[0.0,0.0,0.0]):
            self.__neighbor=neighbor
            if not per==[0.0,0.0,0.0]: self.__per=per

        def reset_periodicity(self,per):
            self.__per=per

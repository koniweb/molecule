#!/usr/bin/env python
##########################################################################
# Molecule Class
# subclasses:
#   Atom Class
#   Bond Class
##########################################################################
version=3.3
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

ndim=calc.ndim
pse=[
    # pse[element][0] name
    # pse[element][1] type
    # pse[element][2] ~ weight
    # pse[element][3] pseudopotential
        ["X",  0, 0.0,"pseudopotential"],
        ["H",  1, 1.0079, "PSEUDO"], ["He",2,4.0026,   "PSEUDO"],
        
        ["Li", 3, 6.941,  "PSEUDO"], ["Be", 4, 9.0122, "PSEUDO"], 
        ["B" , 5,10.811,  "PSEUDO"], ["C" , 6,12.0107, "PSEUDO"], ["N" , 7,14.007, "PSEUDO"], 
        ["O" , 8,15.999,  "PSEUDO"], ["F" , 9,18.998,  "PSEUDO"], ["Ne",10,20.18,  "PSEUDO"],
        
        ["Na",11,22.99,   "PSEUDO"], ["Mg",12,24.305,  "PSEUDO"],  
        ["Al",13,26.982,  "PSEUDO"], ["Si",14,28.086,  "PSEUDO"], ["P" ,15,30.974, "PSEUDO"], 
        ["S", 16,32.065,  "PSEUDO"], ["Cl",17,35.453,  "PSEUDO"], ["Ar",18,39.948, "PSEUDO"],

        ["K", 19,39.098,  "PSEUDO"], ["Ca",20,40.078,  "PSEUDO"], 
        # row transition metals
        ["Ga",31,69.723,  "PSEUDO"], ["Ge",32,72.64,   "PSEUDO"], ["As",33,74.922, "PSEUDO"], 
        ["Se",34,78.96,   "PSEUDO"], ["Br",35,79.904,  "PSEUDO"], ["Kr",36,83.798, "PSEUDO"]
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
        self.__typelist=[] # which element
        self.__id=id(self)
        # set atom list
        self.__at=[]
        self.__celldm=1.0
        self.__vec=[[0.0 for x in xrange(0,ndim)]for x in xrange(0,ndim)] 
        # vec[0]: a
        # vec[1]: b
        # vec[2]: c
        self.__offset=[0.0,0.0,0.0]

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
        return returndata
    def number2element(self,number):
        returndata=[]
        for ielement in self.pse():
            if number==ielement[1]: returndata=ielement
        return returndata  
    def weight2element(self,weight):
        returndata=[]
        for ielement in self.pse():
            if float(weight)==float(ielement[2]): returndata=ielement
        return returndata
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
            else: tid=self.typelist().index([type[1],name])
        return tid

    # set typelist
    def set_typelist2list(self,list):
        self.__typelist=list

    # set atomlist
    def set_atomlist(self,list):
        self.__at=list

    # clear atom list
    def clear_atoms(self):
        self.__at=[]

    # append atom
    def append_atom_coo(self,type,x,y,z):
        i=len(self.at())
        if type.isdigit(): self.__at.append(self.atom(self,i,""  ,int(type),
                                                      float(x),float(y),float(z)))
        else:              self.__at.append(self.atom(self,i,type,-1       ,
                                                      float(x),float(y),float(z)))
        return

    # append an instance atom
    def append_atom_cp(self,addat):
        self.__at.append(self.atom(addat.mol(),
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
    def append_atom(self,atom):
        self.__at.append(atom)

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

    # set natoms ntypes and filename
    def set(self,filename="",filemolnumber=0,comment=""):
        self.__file=filename
        self.__filemolnumber=filemolnumber
        self.__comment=comment
        # get ntypes
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
        print 'setting atom {:d} - {:d} type {:d} ...'.format(st,end,type)
        for cntat in range(st,end):
            self.at()[cntat].set_number(int(type))
        return

    # set name range
    def set_name(self,name,st=0,end=0):
        if end <= 0: 
            end = self.natoms() + end
            if end <=0: end = self.natoms()
        print 'setting atom {:d} - {:d} name {:s} ...'.format(st,end,name)
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
        angle=angle*2*math.pi/360.0
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
            self.__neighbor=neibg
            # define periodicty
            self.__per=per

        #############################################################
        # return functions
        #############################################################

        def atom(self):
            return self.__atom

        def neighbor(self):
            return self.__neighbor
            
        def bondlength(self):
            return calc.a_dist(self.atom(),self.neighbor(),per)

        def bondvector(self):
            return calc.a_vec(self.atom(),self.neighbor(),per)
            
        #############################################################
        # set functions
        #############################################################
        def reset_bond(self,neighbor,per=[0.0,0.0,0.0]):
            self.__neighbor=neighbor
            if not per==[0.0,0.0,0.0]: self.__per=per

        def reset_periodicity(self,per):
            self.__per=per
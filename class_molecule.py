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
from operator import itemgetter, attrgetter
import mod_calc as calc # several functions
import copy

ndim=calc.ndim
types=[
    # types[type][0] name
    # types[type][1] id
    # types[type][2] ~ weight
    # types[type][3] pseudopotential
        ["X",  0, 0.0,"pseudopotential"],
        ["H",  1, 1.0079, "P"], ["He",2,4.0026,   "P"],
        
        ["Li", 3, 6.941,  "P"], ["Be", 4, 9.0122, "P"], 
        ["B" , 5,10.811,  "P"], ["C" , 6,12.0107, "P"], ["N" , 7,14.007, "P"], 
        ["O" , 8,15.999,  "P"], ["F" , 9,18.998,  "P"], ["Ne",10,20.18,  "P"],
        
        ["Na",11,22.99,   "P"], ["Mg",12,24.305,  "P"],  
        ["Al",13,26.982,  "P"], ["Si",14,28.086,  "P"], ["P" ,15,30.974, "P"], 
        ["S", 16,32.065,  "P"], ["Cl",17,35.453,  "P"], ["Ar",18,39.948, "P"],

        ["K", 19,39.098,  "P"], ["Ca",20,40.078,  "P"], 
        # row transition metals
        ["Ga",31,69.723,  "P"], ["Ge",32,72.64,   "P"], ["As",33,74.922, "P"], 
        ["Se",34,78.96,   "P"], ["Br",35,79.904,  "P"], ["Kr",36,83.798, "P"]
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
class molecule(mxyz.molecule_rw,mpw.molecule_rw,mlmp.molecule_rw,mext.molecule_extend):

    # initialize
    def __init__(self):           
        # other inits
        mpw.molecule_rw.__pwscfinit__(self)
        # pse
        self.__types=copy.deepcopy(types)
        # set molecule info
        self.__typelist=[]
        self.file=""
        self.filemolnumber=0
        self.comment=""
        self.__id=id(self)
        # set atom list
        self.__at=[]
        self.vec=[[0.0 for x in xrange(0,ndim)]for x in xrange(0,ndim)]
        self.offset=[0.0,0.0,0.0]

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
    def Rvec(self):
        return self.vec       
        
    # return offset vector
    def Roffset(self):
        return self.offset
        
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
    def Rcomment(self):
        return self.comment

    # return comment
    def types(self):
        return self.__types

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
    
    # set typelist
    def set_typelist(self,list):
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
        self.__at.append(self.atom(self,i,"test",-1,"0.0","0.0","0.0"))
        self.__at[i].set(type,x,y,z)
        return
    # append an instance atom
    def append_atom_cp(self,addat):
        self.__at.append(self.atom(addat.mol,addat.id,addat.name,addat.number,
                                   addat.coord[0],addat.coord[1],addat.coord[2],
                                   addat.mult[0],addat.mult[1],addat.mult[2])
                       )
        self.__at[self.natoms()].charge=addat.charge
        return
    def append_atom(self,atom):
        self.__at.append(atom)

    # set pseudopotential
    def set_type(self,identification,mass=-1.0,pseudopotential=""):
        if mass>0.0:
            self.__types[tid][2]=float(mass)
        # set pseudopotential
        if not pseudopotential=="":
            # if name calculate id
            if type(identification)==type(''): tid=self.type_name2number(name)
            else: tid=identification
            self.__types[tid][3]=pseudopotential

    #############################################################
    # modify functions
    #############################################################    
    # append molecule
    def append_mol(self,molecule,shiftv=[0.0,0.0,0.0],
                rotangle=0.0,rota=[1.0,0.0,0.0],rotp=[0.0,0.0,0.0]):
        # rotate and shift molecules
        molecule.rot(rotangle,rota[0],rota[1],rota[2],rotp[0],rotp[1],rotp[2])
        molecule.shift(shiftv[0],shiftv[1],shiftv[2])
        # add to new molecule
        for iat in range(0,molecule.natoms()):
            self.append_atom_cp(molecule.at()[iat])
        # shift and rotate molecules back
        molecule.shift(-shiftv[0],-shiftv[1],-shiftv[2])
        molecule.rot(-rotangle,rota[0],rota[1],rota[2],rotp[0],rotp[1],rotp[2])
        return

    # set natoms ntypes and filename
    def set(self,filename="",filemolnumber=0,comment=""):
        self.file=filename
        self.filemolnumber=filemolnumber
        self.comment=comment
        # get ntypes
        if any(x.tid == -1 for x in self.at()): self.set_ntypes()
        else: self.set_ntypes_tid()
        return

    # set ntypes
    def set_ntypes(self):
        for cnt in range(self.natoms()):
            for cnttype in range(self.ntypes()):
                if self.at()[cnt].number==self.typelist()[cnttype][0]:
                    self.at()[cnt].tid=cnttype
                    break
            else:
                self.typelist_append(self.at()[cnt].number)
                self.at()[cnt].tid=self.ntypes()-1
        return
    # set ntypes if tid already defined
    def set_ntypes_tid(self):
        for cnt in range(self.natoms()):
            for cnttype in range(self.ntypes()):
                if (self.at()[cnt].number==self.typelist()[cnttype][0] and 
                    self.at()[cnt].tid==cnttype):
                    break
            else:
                self.typelist_append(self.at()[cnt].number)
        return
                           
    # set periodicity
    def set_periodicity(self,vec0,vec1,vec2,off=[0.0,0.0,0.0]):
        self.vec[0]=vec0
        self.vec[1]=vec1
        self.vec[2]=vec2
        self.offset=off
        return
    
    # set periodicity vectors separately
    def set_vecs(self,a=[0.0,0.0,0.0],b=[0.0,0.0,0.0],c=[0.0,0.0,0.0],offset=[0.0,0.0,0.0]):
        if not (a[0]==0.0 and a[1]==0.0 and a[2]==0.0):
            self.vec[0]=a
        if not (b[0]==0.0 and b[1]==0.0 and b[2]==0.0):
            self.vec[1]=b
        if not (c[0]==0.0 and c[1]==0.0 and c[2]==0.0):
            self.vec[2]=c
        o=offset
        if not (o[0]==0.0 and o[1]==0.0 and o[2]==0.0):
            self.offset=offset
        return

    # set number range
    def set_number(self,type,st=0,end=0):
        if end <= 0: 
            end = self.natoms() + end
            if end <=0: end = self.natoms()
        print 'setting atom {:d} - {:d} type {:d} ...'.format(st,end,type)
        for cntat in range(st,end):
            self.at()[cntat].snumber(int(type))
        return

    # set name range
    def set_name(self,name,st=0,end=0):
        if end <= 0: 
            end = self.natoms() + end
            if end <=0: end = self.natoms()
        print 'setting atom {:d} - {:d} name {:s} ...'.format(st,end,name)
        for cntat in range(st,end):
            self.at()[cntat].sname(name)
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
            coo=atom.coord
            x= mat[0][0]*coo[0] + mat[1][0]*coo[1] + mat[2][0]*coo[2]
            y= mat[0][1]*coo[0] + mat[1][1]*coo[1] + mat[2][1]*coo[2]
            z= mat[0][2]*coo[0] + mat[1][2]*coo[1] + mat[2][2]*coo[2]
            coo[0]=x
            coo[1]=y
            coo[2]=z
        # move molecule back after rotation
        self.shift(px,py,pz)

    # return typedata
    def type_number2weight(self,number):
        return types[number][2]
    
    # return typenumber
    def type_name2number(self,name):
        # types[type][0] name
        # types[type][1] id
        # types[type][2] ~ weight
        # types[type][3] pseudopotential
        number=0
        for cnt in range(len(types)):
            if name==types[cnt][0]:
                number=types[cnt][1]
        return number

    # return typenumber
    def type_weight2name(self,weight):
        # types[type][0] name
        # types[type][1] id
        # types[type][2] ~ weight
        # types[type][3] pseudopotential
        name=""
        for cnt in range(len(types)):
            if weight==types[cnt][2]:
                name=types[cnt][0]
        return name

    # stretch structure
    def stretch(self,factor):
        # stretch atoms
        for atom in self.at():
            for d in range(ndim): 
                atom.coord[d]=atom.coord[d]*factor[d]
        for d in range(ndim):
            for vcnt in range(ndim):
                self.vec[vcnt][d]=self.vec[vcnt][d]*factor[d]
            self.offset[d]=self.offset[d]*factor[d]    
        return

######################################################################
# ATOM CLASS
######################################################################
    class atom:   
        # initialize
        def __init__(self,mol,id,name,number,x,y,z,mx=0,my=0,mz=0,atomcharge=0.0,tid=-1):
            self.mol=mol
            self.id=id
            self.name=name
            self.number=int(number)
            self.tid=int(tid)
            self.coord=[x,y,z]
            self.mult=[mx,my,mz]
            self.charge=atomcharge

        #############################################################
        # return functions
        #############################################################                      
        def Rid(self):
            return self.id

        def Rname(self):
            return self.name

        def Rnumber(self):
            return self.number

        def Rtid(self):
            return self.tid

        def Rcoord(self):
            return self.coord
        #############################################################
        # return functions
        #############################################################                      
        # set positon and name/number
        def set(self,type,x,y,z):
            if type.isdigit():
                self.sname(type)
                self.snumber(int(type))
                self.spos(x,y,z)
            else:
                self.sname(type)
                self.snumber(self.name2number(type))
                self.spos(x,y,z)
    
        # set position of atom
        def spos(self,x,y,z):
            self.coord=[x,y,z]
    
        # set name of atom
        def sname(self,nameinput):
            self.name=nameinput
    
        # def number of atom
        def snumber(self,numberinput):
            self.number=numberinput
                
        # converts name to a number
        def name2number(self,name):
            # types[type][0] name
            # types[type][1] id
            # types[type][2] ~ weight
            number=0
            for cnt in range(len(types)):
                if name==types[cnt][0]:
                    number=types[cnt][1]
            return number

        # converts number to name
        def number2name(self,number):
            # types[type][0] name
            # types[type][1] id
            # types[type][2] ~ weight
            name=""
            for cnt in range(len(types)):
                if number==types[cnt][1]:
                    name=types[cnt][0]
            return name
    
        # converts number to name
        def weight2number(self,weigth):
            # types[type][0] name
            # types[type][1] id
            # types[type][2] ~ weight
            number=""
            for cnt in range(len(types)):
                if weight==types[cnt][2]:
                    number=types[cnt][1]
            return number
    
        # shift atom
        def shift(self,x,y,z):
            self.coord=[self.coord[0]+x,self.coord[1]+y,self.coord[2]+z]

       

        
######################################################################
# BOND CLASS
######################################################################       
    class bond:
        def __init__(self,atom,neighbor,per=[0,0,0],v=[0.0,0.0,0.0]):
            if ( v == [0.0,0.0,0.0] ): 
                v=calc.a_vec(atom,neighbor,per)
                d=calc.a_dist(atom,neighbor,per)
            else:
                d=calc.length(v)
            self.atom=int(atom.id)
            self.neighbor=int(neighbor.id)
            self.per=per
            self.dist=d
            self.vec=v
            self.vec_trans=[0.0,0.0,0.0]

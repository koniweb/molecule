#!/usr/bin/env python
##########################################################################
# writelmp
# readlmp
# readlmpcustomout
##########################################################################
version=1.1
versiontext='# mod_lmp.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import copy
import re
import mod_calc as calc

#----------------------------------------------------------------------
# module mod_xyz
#----------------------------------------------------------------------
class molecule_rw:
    # write lammps file
    def writelmp(self,filename="",status='w',lcharge=False,lmoltype=False):
        ndim=calc.ndim
        # define vec and offset
        v=self.vec()
        o=self.offset()
        # check if vector is correct for lammps
        if not (v[0][1]==0.0 and v[0][2]==0.0 and v[1][2]==0.0): 
            print >> sys.stderr, "vectors are not defined correctly for a lammps file"
            quit()
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        mol=self
        print >>f
        print >>f
        # numbers
        print >>f, ("{:d} atoms".format(mol.natoms()))
        print >>f
        # types
        print >>f, ("{:d} atom types".format(mol.ntypes()))
        print >>f
        # vectors
        if not (v[1][0]==0.0 and v[2][0]==0.0 and v[2][1]==0.0):
            print >>f, ("{:f} {:f} {:f} xy xz yz".format(v[1][0],v[2][0],v[2][1]))
            print >>f
        # offset
        print >>f, ("{:20.10f} {:20.10f} xlo xhi".format(o[0],v[0][0]+o[0]))
        print >>f, ("{:20.10f} {:20.10f} ylo yhi".format(o[1],v[1][1]+o[1]))
        print >>f, ("{:20.10f} {:20.10f} zlo zhi".format(o[2],v[2][2]+o[2]))
        print >>f
        print >>f, "Masses"
        print >>f
        for cnttype in range(mol.ntypes()):
            if not self.typelist()[cnttype][1]=="":
                name=self.typelist()[cnttype][1] # name from typelist              
            else:
                name=self.number2element(mol.typelist()[cnttype][0])[0] # name from pse
            print >>f ,(
                '{:6d} {:f} #{:s}'.format(
                    cnttype+1,
                    self.number2element(mol.typelist()[cnttype][0])[2], # weight
                    name,
                    )
                )       
        print >>f
        # atoms
        print >> f, "Atoms"
        print >> f
        for cntat in range(mol.natoms()):
            a=mol.at()[cntat]
            moltype=""
            charge=""
            if lcharge:  charge="{:15.10f}".format(a.charge()) 
            if lmoltype: moltype='{:4d}'.format(1)
            print >>f ,(
                '{:6d} {:s} {:4d} {:s} {:15.10f} {:15.10f} {:15.10f}'.format(
                    cntat+1, moltype, a.tid()+1, charge, a.coord()[0], a.coord()[1], a.coord()[2])
                )
            
        print >>f
        if filename != "": f.close()
        return
    
    # read molecules in lammps file
    def readlmp(self,filename,lcharge=False,lmoltype=False):
        ndim=calc.ndim
        # set molecule
        moltype=0
        charge=0.0
        molecules=[]
        tilt=[0.0 for j in range(ndim)]
        # set columns to read in
        ctid=1
        ccoord=2
        if lcharge: 
            ccharge=2
            ccoord+=1
        if lmoltype: 
            cmoltype=1
            ccoord+=1
            ctid+=1
            ccharge+=1
        # check read file
        try: 
            file=open(filename, 'r')
        except IOError:
            print >> sys.stderr, "... input file not found"
            exit()
        # read file
        cntline=0
        natoms=0
        opt=""
        for line in file:
            cntline+=1
            linesplit=line.split()
            # create new molecule
            if (cntline)==1:
                mol=self.__class__()
                mol.clear_atoms()
                cntat=0
                cnttype=0
            # check commented lines
            if len(linesplit)>0 and linesplit[0][0]=="#":
                continue
            # check for empty lines
            if len(linesplit)==0:
                continue
            # read atoms
            if opt=="atoms":
                # find atomtype and so on
                tid    = int(linesplit[ctid])
                # find name and number out of typelist
                for t in range(len(mol.typelist())):
                    if int(re.sub("[a-zA-Z]","", mol.typelist()[t][1]))==tid:
                        name=str(mol.pse()[mol.typelist()[t][0]][0])+str(t+1)
                        number=t
                        break
                # read charge
                if lcharge: charge=float(linesplit[ccharge])
                # read molid
                if lmoltype: moltype=int(linesplit[cmoltype])
                # append atoms
                mol.append_atom(
                    self.__class__.atom(
                        mol,
                        cntat,
                        name,
                        number,
                        float(linesplit[ccoord+0]),
                        float(linesplit[ccoord+1]),
                        float(linesplit[ccoord+2]),
                        tid=tid-1,
                        atomcharge=charge
                        )
                    )
                cntat+=1
                if cntat==natoms: opt=""
            # read types
            if opt=="masses":
                # append types
                type=mol.weight2element(float(linesplit[1]))
                if len(type)==0: 
                    print >> sys.stderr, "ERROR: type could not be found"
                    exit()
                mol.typelist_append(
                    type[1],type[0]+linesplit[0])
                # until cnttype==ntypes
                cnttype+=1
                if cnttype==ntypes: opt=""
            #
            # read main options
            #
            # opt for blocks
            if len(linesplit)>0:
                option=linesplit[0]
                if option=="Masses":
                    opt="masses"
                elif option=="Atoms":
                    opt="atoms"
                elif option=="Bonds":
                    opt="bonds"
            # NUMBERS
            if len(linesplit)>=2:                
                option=linesplit[1]
                if   option=="atoms":
                    natoms=int(linesplit[0])
                elif option=="bonds":
                    nbonds=int(linesplit[0])
                elif option=="angles":
                    nangles=int(linesplit[0])
                elif option=="dihredrals":
                    ndihredrals=int(linesplit[0])
                elif option=="impropers":
                    nimpropers=int(linesplit[0])
            # TYPES
            if len(linesplit)>=3:                
                option=linesplit[1:3]
                if   option==["atom","types"]:
                    ntypes=int(linesplit[0])
                elif option==["bond","types"]:
                    tbonds=int(linesplit[0])
                elif option==["angle","types"]:
                    tangle=int(linesplit[0])
                elif option=="dihredrals":
                    tdihredrals=int(linesplit[0])
                elif option=="impropers":
                    timpropers=int(linesplit[0])
            # Vectors
            if len(linesplit)>=4:
                option=linesplit[2:4]
                if   option==["xlo","xhi"]:
                    x=[float(linesplit[0]),float(linesplit[1])]
                elif option==["ylo","yhi"]:
                    y=[float(linesplit[0]),float(linesplit[1])]
                elif option==["zlo","zhi"]:
                    z=[float(linesplit[0]),float(linesplit[1])]
            if len(linesplit)==6:
                if linesplit[3:6]==["xy","xz","yz"]:
                    tilt=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2])]
        # calculate vectors a,b,c out of the box and tilt
        vec=self.calc_vec(x,y,z,tilt)
        # set periodicity
        mol.set_vecs(vec[0],vec[1],vec[2],[x[0],y[0],z[0]])
        # set molecule and append
        mol.set(filename,1,"")
        molecules.append(copy.copy(mol))
        # close file
        file.close()
        # return molecules
        return molecules

    # read lmp custom
    def readlmpcustomout(self,filename,start=1,end=-1):
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
        cntat=0
        natoms=0
        cntmol=0
        #block
        block=[]
        box=[]
        customdata=[]
        adddata=[]
        adddata_info=[]
        for line in file:
            linesplit=line.split()
            # go through blocks
            if linesplit[0]!="ITEM:":
                if   len(block)==3 and block[0:3]==["NUMBER","OF","ATOMS"]:
                    natoms=int(linesplit[0])
                elif len(block)>=2 and block[0:2]==["BOX","BOUNDS"]:
                    for i in range(len(linesplit)):linesplit[i]=float(linesplit[i])
                    box.append(linesplit)
                elif len(block)>=1 and block[0]=="ATOMS":
                    number=self.name2element(linesplit[data[0]])[0]
                    mol.append_atom(
                        self.__class__.atom(
                            mol,
                            cntat,
                            linesplit[data[0]],
                            number,
                            float(linesplit[data[1]]),
                            float(linesplit[data[2]]),
                            float(linesplit[data[3]])
                            )
                        )
                    # additional data
                    for i in range(4,len(data)):
                        adddata[i-4][1].append(linesplit[data[i]])
                    # finish atom read in
                    cntat+=1
                    if len(mol.at())==natoms:
                        block=[]
            # get block information
            else:
                block=linesplit[1:]
                # get custom data information
                if linesplit[1]=="ATOMS":
                    data=[-1,-1,-1,-1] # name,x,y,z
                    customdata=linesplit[2:]
                    for i in range(len(customdata)):
                        if   customdata[i]=="element":  data[0]=i
                        elif customdata[i]=="x":        data[1]=i
                        elif customdata[i]=="y":        data[2]=i
                        elif customdata[i]=="z":        data[3]=i
                        elif customdata[i]=="id":       nid=i
                        else: data.append(i)
                    # create molecules
                    cntat=0
                    if cntmol==0:
                        mol=self.__class__()
                        mol.clear_atoms()
                    else:
                        mol.set(filename,cntmol)
                        mol.set_data(adddata)
                        # calculate vectors a,b,c out of the box and tilt
                        if len(box[0])==2 and len(box[1])==2 and len(box[2])==2: tilt=[0.0,0.0,0.0]
                        else: tilt=[box[0][2],box[1][2],box[2][2]]
                        vec=self.calc_vec(box[0][0:2],box[1][0:2],box[2][0:2],tilt)
                        mol.set_vecs(vec[0],vec[1],vec[2],[box[0][0],box[1][0],box[2][0]])
                        # append molecules
                        molecules.append(copy.copy(mol))
                        adddata=[]
                    cntmol+=1
                    # additional data
                    for i in range(4,len(data)):
                        adddata.append([customdata[data[i]],[]])
        # append last frame
        mol.set(filename,cntmol)
        mol.set_data(adddata)
        # calculate vectors a,b,c out of the box and tilt
        if len(box[0])==2 and len(box[1])==2 and len(box[2])==2: tilt=[0.0,0.0,0.0]
        else: tilt=[box[0][2],box[1][2],box[2][2]]
        vec=self.calc_vec(box[0][0:2],box[1][0:2],box[2][0:2],tilt)
        mol.set_vecs(vec[0],vec[1],vec[2],[box[0][0],box[1][0],box[2][0]])
        # append molecule
        molecules.append(copy.copy(mol))
        # close file 
        file.close()
        # return molecules                                                                                         
        return molecules

    # calculate vectors out of box and tilt
    def calc_vec(self,x,y,z,tilt):
        print "NOW",x,y,z # DEBUG
        vec=[]
        vec.append([ x[1]-x[0], 0.0,       0.0])
        vec.append([ tilt[0],   y[1]-y[0], 0.0])
        vec.append([ tilt[1],   tilt[2],   z[1]-z[0] ])
        return vec

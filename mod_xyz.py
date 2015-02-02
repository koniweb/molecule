#!/usr/bin/env python
##########################################################################
# writexyz
# readfxyz
##########################################################################
version=3.4
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
    def writexyz(self,filename="",status='w',extended=False,data=[]):
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        mol=self
        if extended==False: data.append(["",[""]*self.natoms()])
        # EXTENDED
        else:
            # empty data field
            if len(data)==0:
                # if additional data present
                if len(mol.data())!=0:
                    data=mol.data()
                else:
                    data.append(["",[""]*self.natoms()])
            # add data to data array if not already done
            for val in data:
                if len(val)==1:
                    val.append([""]*self.natoms())
                    # loop over atoms get data
                    for iat in range(self.natoms()):
                        val[1][iat]=str(  getattr(self.at()[iat],str(val[0]))()  )
            # create commentline
            mol.set_comment(mol.exyz_writecomment(data))
        # print output
        print >>f, mol.natoms()
        print >>f, mol.comment()
        for cntat in range(0,mol.natoms()):
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f} {:s}'.format(
                    mol.at()[cntat].type()[0], 
                    mol.at()[cntat].coord()[0], 
                    mol.at()[cntat].coord()[1], 
                    mol.at()[cntat].coord()[2],
                    " ".join([row[1][cntat] for row in data])
                    )
                )
        if filename != "": f.close()
        return

    # read molecules in xyz file
    def readxyz(self,filename,start=1,end=-1,extended=False):
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
                comment=line.strip()
                # EXTENDED
                # strip comment, set molecule and append it
                if extended: exyzdata=mol.exyz_readcomment(comment)
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
                # EXTENDED
                # read exyz data
                if extended: mol.exyz_readdata(linesplit,exyzdata)
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

    def exyz_readcomment(self,comment):
        commentsplit=comment.strip().split(" ")
        commentsplit=filter(bool, commentsplit)
        exyzdata=[]
        i=0
        while i<len(commentsplit):
            str=commentsplit[i]
            # select data for molecule
            if str=="E":
                self.set_energy(float(commentsplit[i+1]))
                i+=1
            elif str=="a":
                self.set_vecs(a=[float(commentsplit[i+1]), 
                                 float(commentsplit[i+2]), 
                                 float(commentsplit[i+3])]
                              )
                i+=3
            elif str=="b":
                self.set_vecs(b=[float(commentsplit[i+1]), 
                                 float(commentsplit[i+2]), 
                                 float(commentsplit[i+3])]
                              )
                i+=3
            elif str=="c":
                self.set_vecs(c=[float(commentsplit[i+1]), 
                                 float(commentsplit[i+2]), 
                                 float(commentsplit[i+3])]
                              )
                i+=3
            elif str=="off":
                self.set_vecs(off=[float(commentsplit[i+1]), 
                                   float(commentsplit[i+2]), 
                                   float(commentsplit[i+3])]
                              )
                i+=3
            # select data per atom
            else:
                exyzdata.append(str)
            i+=1
        return exyzdata

    def exyz_readdata(self,linesplit,exyzdata):
        # 0:   Atomtype
        # 1-3: Coordinates
        # get datafields to fill
        Dsplit=[row[0] for row in exyzdata]
        # check if all data is there
        if len(Dsplit)==len(linesplit)-4:           
            for i in range(len(Dsplit)):
                if Dsplit[i]=="charge":
                    self.at()[self.natoms()-1].set_charge(float(linesplit[i+4]))
                #elif Dsplit[i]=="datavaluename""
                #   self.at()[self.natoms()-1].set_datavaluename(float(linesplit[i+4]))
        else:
            print >> sys.stderr, "... extended xyz data not present"
            exit()
        return

    def exyz_writecomment(self,exyzdata):
        # add vectors to commentline
        a=""
        b=""
        c=""
        off=""
        if self.vec()[0]!=[0.0,0.0,0.0]:
            a="a {:15.10f} {:15.10f} {:15.10f}".format(
                self.vec()[0][0],self.vec()[0][1],self.vec()[0][2])
        if self.vec()[1]!=[0.0,0.0,0.0]:
            b="b {:15.10f} {:15.10f} {:15.10f}".format(
                self.vec()[1][0],self.vec()[1][1],self.vec()[1][2])
        if self.vec()[2]!=[0.0,0.0,0.0]:
            c="c {:15.10f} {:15.10f} {:15.10f}".format(
                self.vec()[2][0],self.vec()[2][1],self.vec()[2][2])
        if self.offset()!=[0.0,0.0,0.0]:
            off="off {:15.10f} {:15.10f} {:15.10f}".format(
                self.offset()[0],self.offset()[1],self.offset()[2])
        # add energy to comment line
        E=""
        if self.energy()!=0.0:E="E {:.15f}".format(self.energy())
        # return comment line
        return "{:s} {:s} {:s} {:s} {:s} {:s}".format(a,b,c,off,E," ".join([row[0] for row in exyzdata]) ).strip()


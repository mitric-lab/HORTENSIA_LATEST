#!/usr/bin/env python

import string

import numpy as np

import hortensia_latest.wigner.atom as atom
import hortensia_latest.wigner.HarmonicBase as HarmonicBase

class HarmonicMolpro(HarmonicBase.HarmonicBase):
    def __init__(self,linear,filename1,filename2):
        """
        This is a class handling Molpro hessian and normal mode calculations
        """
        HarmonicBase.HarmonicBase.__init__(self,linear)
        self.readFile(filename1)

    def setStructure(self):
        self.sym = []
        for i in range(len(self.lines)):
            if "geometry={" in self.lines[i]:
                j = i
                break
                
        self.nat = int(self.lines[j+1].split()[0])
        self.coordinates = []
        for i in range(j+2, j+2+self.nat):
            h = self.lines[i].split()
            self.sym.append(self.lines[i].split()[0])
            
            for k in [2,4,6]:
                self.coordinates.append(float(h[k])/self.autoangs)
                
        self.coordinates = np.array(self.coordinates)
        self.ndim = len(self.coordinates)

    def setMasses(self):
        self.masses = []
        for i in range(self.ndim):
            self.masses.append(atom.atomicmasses[self.sym[i/3]])
        self.masses = np.array(self.masses)

    def setHessian(self):
        def getInd(s):
            #converts indices GX1,GY1,GZ1,GX2 etc. to integers 0,1,2,3...
            n1 = int(s[2:]) - 1
            if s[1] == "X":
                n2 = 0
            elif s[1] == "Y":
                n2 = 1
            else:
                n2 = 2
                
            return 3*n1 + n2

        self.hessian = np.zeros((self.ndim, self.ndim))
        # reading hessian from Molpro output
        for i in range(len(self.lines)):
            if "Force Constants (Second Derivatives of the Energy)" in self.lines[i]:
                ind1 = i
            if "Atomic Masses" in self.lines[i+3]:
                ind2 = i+1
                break
                
        flag = 0
        for i in range(ind1+1, ind2):
            h = self.lines[i].split()
            if ("G" in h[0]) and (len(h) == 1):
                flag = 1
            elif ("G" in h[0]) and ("G" in h[1]):
                flag = 1
            if flag == 1:
                hessind2 = []
                for ii in range(len(h)):
                    hessind2.append(getInd(h[ii]))
                hessind2 = np.array(hessind2)
                flag = 0
            else:
                hessind1 = getInd(h[0])
                for ii in range(1, len(h)):
                    self.hessian[hessind1][hessind2[ii-1]] = \
                        float(string.replace(h[ii],"D","E"))
                    self.hessian[hessind2[ii-1]][hessind1] = \
                        self.hessian[hessind1][hessind2[ii-1]]
#!/usr/bin/env python

import numpy as np

import hortensia_latest.wigner.atom as atom
import hortensia_latest.wigner.HarmonicBase as HarmonicBase

class HarmonicMndo(HarmonicBase.HarmonicBase):
    def __init__(self,linear,filename1,filename2):
        """
        This is a class handling mndo99 hessian and normal mode calculations
        """
        HarmonicBase.HarmonicBase.__init__(self,linear)
        self.linear       = linear
        self.atomicmasses = atom.atomicmasses
        self.readFile(filename2)
        with open(filename1) as f:
            self.lines2 = f.readlines()

    def setStructure(self):
        # changed so as to be compatible with molden.dat and fort.13
        # changed 1.5.14 for molden.dat with coordinates in a.u. 
        # after [FR-COORD]
        self.coordinates = []
        self.sym = []

        flag=0
        for i in range(len(self.lines)):
            if flag == 1:
                if "[" in self.lines[i]:
                    break
                h = self.lines[i].split()
                self.sym.append(h[0])

                for j in range(1,4):
                    self.coordinates.append(float(h[j]))

            if "[FR-COORD]" in self.lines[i]:
                    flag = 1

        self.coordinates = np.array(self.coordinates)
        self.ndim = len(self.coordinates)
        self.nat  = self.ndim/3

    def setMasses(self):
        self.masses = []
        for i in range(len(self.sym)):
            for ii in range(3):
                self.masses.append(self.atomicmasses[self.sym[i]])
        self.masses = np.array(self.masses)

    def setHessian(self):
        flag  = 0
        flag2 = 0
        self.hessian = np.zeros((self.ndim, self.ndim))
        for i in range(len(self.lines2)):
            if flag2 == 1:
                h = self.lines2[i].split()
                firstrow = int(h[0])-1
                flag2 = 0

            if flag == 2:
                h = self.lines2[i].split()
                rowindex = int(h[0]) - 1
                for j in range(1, len(h)):
                    self.hessian[rowindex][firstrow+j-1] = \
                        float(h[j])
                    self.hessian[firstrow+j-1][rowindex] = \
                        self.hessian[rowindex][firstrow+j-1]

                if rowindex == self.ndim-1:
                    flag = 1

            if "CARTESIAN FORCE CONSTANT MATRIX IN ATOMIC UNITS" in self.lines2[i]:
                flag = 1

            if (flag == 1) and ("-----" in self.lines2[i]):
                flag  = 2
                flag2 = 1

    def readTransformationMatrix(self):
        # read eigenvectors directly from molden.dat, but which units are used??
        flag = 0
        ind  = 0
        modeindex = -1
        if self.linear == 0:
            nrows = self.ndim - 6
        else:
            nrows = self.ndim - 5

        self.lmat = np.zeros((nrows, self.ndim))
        for i in range(len(self.lines)-1):
            if flag == 1:
                h = map(float, self.lines[i].split())
                for ii in range(len(h)):
                    self.lmat[modeindex][ind] = h[ii]
                    ind += 1
            if "vib" in self.lines[i]:
                flag = 1
                ind  = 0
                modeindex += 1
            if "vib" in self.lines[i+1]:
                flag = 0
            if "[" in self.lines[i+1]:
                flag = 0

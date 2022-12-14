#!/usr/bin/env python3

import string
import numpy as np

import hortensia_latest.wigner.HarmonicBase as HarmonicBase

class HarmonicGaussian(HarmonicBase.HarmonicBase):
    def __init__(self, linear, filename1, filename2):
        """
        This is a class handling Gaussian hessian and
        normal mode calculations
        """
        HarmonicBase.HarmonicBase.__init__(self, linear)
        self.readFile(filename1)
        self.suffix    = filename1[-4:]
        self.filename2 = filename2


    def setStructure(self):
        if self.suffix == "fchk":
            self.setStructureChk()
        else:
            self.setStructureLog()


    def setStructureLog(self):
        flag = 0
        self.sym = []
        for i in range(2, len(self.lines)-1):
            if flag == 1:
                self.sym.append(self.lines[i].split()[0])
                if len(self.lines[i+1].split()) != 4:
                    flag = 0
            if "Symbolic Z-matrix" in self.lines[i-1]:
                flag = 1
            if "Standard orientation" in self.lines[i]:
                startind = i

        flag = 1
        self.coordinates = []
        for i in range(startind+5, len(self.lines)):
            if "-------------" in self.lines[i]:
                flag = 0
            if flag == 1:
                h = self.lines[i].split()
                self.coordinates.append(float(h[3])/self.autoangs)
                self.coordinates.append(float(h[4])/self.autoangs)
                self.coordinates.append(float(h[5])/self.autoangs)

        self.coordinates = np.array(self.coordinates)
        self.ndim = len(self.coordinates)
        self.nat  = self.ndim/3


    def setStructureChk(self):
        flag = 0
        f = open(self.filename2)
        li = f.readlines()
        f.close()
        self.sym = []
        for i in range(2,len(li)-1):
            if flag == 1:
                self.sym.append(li[i].split()[0])
                if len(li[i+1].split()) != 4:
                    flag = 0

            if flag == 2:
                self.sym.append(li[i].split(",")[0])
                if len(li[i+1].split()) != 1:
                    flag = 0

            if "Symbolic Z-matrix" in li[i-1]:
                flag = 1

            if "Redundant internal coordinates found in file." in li[i]:
                # this added 18/9/2017 for the case of geom=check
                flag=2

        self.coordinates=[]
        for i in range(len(self.lines)):
            if "Current" in self.lines[i]:
                ncoord = int(self.lines[i].split()[5])

                n = 0
                while n < ncoord:
                    i += 1
                    for k in self.lines[i].split():
                        self.coordinates.append(k)
                    n += len(self.lines[i].split())
                break

        self.coordinates = np.array(self.coordinates, dtype=float)
        self.ndim = len(self.coordinates)
        self.nat  = self.ndim // 3


    def setMasses(self):
        if self.filename2 == "":
            # this option if everything is read from a log file
            li = self.lines
        else:
            # this option if hessian is read from fchk file,
            # the file here should then be the log file
            f  = open(self.filename2)
            li = f.readlines()
            f.close()

        self.masses = []
        for i in range(len(li)):
            if ("has atomic number" in li[i]) and ("and mass" in li[i]):
                self.masses.append(self.toaumass*float(li[i].split()[8]))
                self.masses.append(self.toaumass*float(li[i].split()[8]))
                self.masses.append(self.toaumass*float(li[i].split()[8]))
        self.masses = np.array(self.masses)


    def setHessian(self):
        if self.suffix == "fchk":
            self.setHessianChk()
        else:
            self.setHessianLog()


    def setHessianLog(self):
        self.hessian = np.zeros((self.ndim, self.ndim))

        # reading hessian from Gaussian output
        for i in range(len(self.lines)):
            if "Force constants in Cartesian coordinates:" in self.lines[i]:
                ind1 = i
            if "Force constants in internal coordinates:" in self.lines[i]:
                ind2 = i
        flag = 0
        for i in range(ind1+1, ind2):
            h = self.lines[i].split()
            if (len(h) == 5) and (flag == 0):
                hessind2 = np.array(map(int, h)) - 1
                flag = 1
            else:
                hessind1 = int(h[0]) - 1
                for ii in range(1, len(h)):
                    self.hessian[hessind1][hessind2[ii-1]] = \
                        float(string.replace(h[ii],"D","E"))
                    self.hessian[hessind2[ii-1]][hessind1] = \
                        self.hessian[hessind1][hessind2[ii-1]]

                if (len(h) == 6) and (len(self.lines[i+1].split()) == 5):
                    flag = 0


    def setHessianChk(self):
        self.hessian = np.zeros((self.ndim, self.ndim))

        for i in range(len(self.lines)):
            if "Cartesian Force Constants" in self.lines[i]:
                startind = i+1

        endind = startind + (self.ndim*self.ndim + self.ndim)//10

        # only half the force constant matrix is printed in fchk file
        # (N**2 - N)/2 + N = (N**2 + N)/2 distinct values,
        # arranged in 5 columns per line
        # if this number cannot be divided by 5,
        # one incomplete line has to be added
        if ((self.ndim**2 + self.ndim)/2)%5 != 0:
            endind += 1
        h = []
        for i in range(startind, endind):
            h.append(np.asarray(self.lines[i].split(), dtype=float))

        count = 0
        for i in range(self.ndim):
            for j in range(i+1):
                zeile  = count//5
                spalte = count%5
                self.hessian[i][j] = h[zeile][spalte]
                self.hessian[j][i] = self.hessian[i][j]
                count += 1

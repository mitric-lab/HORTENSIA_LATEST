#!/usr/bin/env python3

import numpy as np

import hortensia_latest.misc as misc
import hortensia_latest.wigner.HarmonicBase as HarmonicBase

class HarmonicQChem(HarmonicBase.HarmonicBase):
    def __init__(self,linear,filename1,filename2):
        """
        This is a class handling QChem hessian and
        normal mode calculations
        
        The first file needs to be the standard output file, the second file the
        HESS file which is only generated when using the -save flag in the qchem
        calculation
        """

        HarmonicBase.HarmonicBase.__init__(self, linear)
        with open(filename1) as f:
            self.f1 = f.readlines()
        with open(filename2) as f:
            self.f2 = f.readlines()


    def setStructure(self):
        for i, j in enumerate(self.f1[-1::-1]):
            if "Standard Nuclear" in j:
                nat = 0
                while self.f1[-i+2+nat].split()[0][0] != "-":
                    nat += 1

                xyz = np.zeros((nat, 3), dtype=float)
                s   = np.zeros(nat, dtype='U2')
                for n, x in enumerate(self.f1[-i+2:-i+2+nat]):
                    xyz[n] = x.split()[2:]
                    s[n]   = x.split()[1]
                break
        
        self.sym  = s
        self.xyz  = xyz
        self.coordinates = self.xyz.flatten()

        self.nat  = nat
        self.ndim = 3 * nat


    def setMasses(self):
        masses = []
        for i in self.sym:
            masses.append(misc.masses[i])
        self.masses = np.array(masses)


    def setHessian(self):
        """
        The HESS file gives the lower triangle of the Hessian matrix with a 
        maximum of 3 values per line
        """

        self.hessian = np.zeros((self.ndim, self.ndim), dtype=float)

        nline, nrow = 0, 0
        for i in self.f2[2:]:
            if i[0] == "$":
                break
            for j in i.split():
                self.hessian[nline, nrow] = j
                self.hessian[nrow, nline] = j
                nrow += 1
            if nrow > nline:
                nline += 1
                nrow   = 0

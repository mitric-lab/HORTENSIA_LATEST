#!/usr/bin/env python

import numpy as np
import numpy.linalg as la

class HarmonicBase:
    def __init__(self, linear):
        """
        This is a basis class for handling hessian and normal mode calculations
        """

        self.linear   = linear
        self.autoangs = 0.5291772
        self.toaumass = 1822.88853


    def readFile(self, filename):
        with open(filename) as f:
            self.lines = f.readlines()


    def setStructure(self):
        pass


    def setMasses(self):
        pass


    def setHessian(self):
        pass


    def setMassWeightedHessian(self):
        for i in range(self.ndim):
            for j in range(self.ndim):
                self.hessian[i][j] /= np.sqrt(self.masses[i]*self.masses[j])


    def setTransformationMatrix(self):
        """
        get eigenvalues (squares of eigenfrequencies) and eigenvectors of
        mass weighted hessian
        lmat is transformation matrix between mass weighted coordinates and
        mass weighted normal coordinates
        """

        eigv      = la.eig(self.hessian)
        self.ind  = np.argsort(eigv[0])
        self.lmat = np.zeros((self.ndim, self.ndim))

        for i in range(self.ndim):
            for j in range(self.ndim):
                # compared to Numeric, indices of eigenvector are interchanged,
                # hence the transpose
                self.lmat[i,j] = eigv[1].T[self.ind[i]][j]

        # normal frequencies for translation and rotation have to be exactly 0
        self.normfreq = np.sort(eigv[0])
        if self.linear == 1:
            m = 5
        else:
            m = 6

        for i in range(len(self.normfreq)):
            if i < m:
                self.normfreq[i] = 0.0
            else:
                self.normfreq[i] = np.sqrt(self.normfreq[i])


    def getTransformationMatrix(self):
        return self.lmat


    def getNormalFrequencies(self):
        return self.normfreq


    def getHessian(self):
        return self.hessian


    def getCoordinates(self):
        return self.coordinates


    def getSymbols(self):
        return self.sym


    def getMasses(self):
        return self.masses

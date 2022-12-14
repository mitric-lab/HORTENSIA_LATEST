#!/usr/bin/env python

import numpy as np

import hortensia_latest.wigner.HarmonicBase as HarmonicBase

class HarmonicTurbomole(HarmonicBase.HarmonicBase):
    def __init__(self,linear,filename1,filename2):
        """
        This is a class handling turbomole hessian and normal mode calculations
        """

        HarmonicBase.HarmonicBase.__init__(self,linear)
        self.readFile(filename1)
        self.filename2 = filename2


    def setStructure(self):
        with open(self.filename2) as f:
            lines = f.readlines()

        self.coordinates = []
        self.sym = []
        for i in range(1, len(lines)):
            if "$" in lines[i]:
                break
            h = lines[i].split()
            self.sym.append(h[3])
            for j in range(3):
                self.coordinates.append(float(h[j]))

        self.coordinates = np.array(self.coordinates)
        self.ndim = len(self.coordinates)
        self.nat  = self.ndim/3


    def setMasses(self):
        self.masses=[]
        for i in range(len(self.lines)):
            if "ATOMIC WEIGHTS" in self.lines[i]:
                    ind1=i
            if "CARTESIAN FORCE" in self.lines[i]:
                    ind2=i
        for i in range(ind1+3,ind2-2):
            self.masses.append(self.toaumass * float(self.lines[i].split()[2]))
            self.masses.append(self.toaumass * float(self.lines[i].split()[2]))
            self.masses.append(self.toaumass * float(self.lines[i].split()[2]))

        self.masses = np.array(self.masses)


    def setHessian(self):

        def parseHessian(element1,element2,irow):
            if (element1 in self.lines[i]) and (element2 not in self.lines[i]):
                tmp = self.lines[i].split()
                ind = tmp.index(element1)
                if element1 == "dx":
                    irow = int(3*(int(tmp[0])-1))

                temphess = []
                for element in tmp[ind+1:]:
                    if "-" in element:
                        if element.index("-") > 0:
                            tmp1=element.split("-")
                            temphess.append(float(tmp1[0]))
                            for el in tmp1[1:]:
                                    temphess.append(-1.0*float(el))

                        if element.index("-") == 0:
                            tmp1=element.split("-")
                            for el in tmp1[1:]:
                                temphess.append(-1.0*float(el))
                    else:
                        temphess.append(float(element))

                for ii in range(len(temphess)):
                    self.hessian[irow][icolumn+ii] = \
                        temphess[ii]
                    self.hessian[icolumn+ii][irow] = \
                        self.hessian[irow][icolumn+ii]

        self.hessian = np.zeros((self.ndim, self.ndim))

        # reading hessian from aoforce output
        for i in range(len(self.lines)):
            if "CARTESIAN FORCE" in self.lines[i]:
                ind1 = i
            if "projected hessian written" in self.lines[i]:
                ind2 = i

        flag = 0
        for i in range(ind1+4, ind2-2):
            if "ATOM" in self.lines[i]:
                icolumn = 3*(int(self.lines[i].split()[1])-1)

            if ("dx" in self.lines[i]) and ("dy" not in self.lines[i]):
                tmp  = self.lines[i].split()
                irow = int(3*(int(tmp[0])-1))
                flag = 1

            if flag == 1:
                parseHessian("dx","dy",irow)
                parseHessian("dy","dx",irow+1)
                parseHessian("dz","dy",irow+2)


    def setDipoleDerivatives(self):

        # cartesian dipole derivatives in a.u. are in ddip file
        # in the order
        #  d mu_x/dx  d mu_x/dy  d mu_x/dz
        #  d mu_y/dx  d mu_y/dy  d mu_y/dz
        #  d mu_z/dx  d mu_z/dy  d mu_z/dz
        # for each atom

        with open("ddip") as f:
            lines1 = f.readlines()

        dipolex = []
        dipoley = []
        dipolez = []

        # reading and mass weighting cartesian dipole derivatives
        for i in range(1, len(lines1)-3, 3):
            h = lines1[i].split()
            for j in range(len(h)):
                dipolex.append(float(h[j]) / np.sqrt(self.masses[(i-1)+j]))
            h = lines1[i+1].split()
            for j in range(len(h)):
                dipoley.append(float(h[j]) / np.sqrt(self.masses[i+j-1]))
            h = lines1[i+2].split()
            for j in range(len(h)):
                dipolez.append(float(h[j]) / np.sqrt(self.masses[i+j-1]))

        # transformation to normal coordinate derivatives
        dipx = np.dot(self.lmat,np.array(dipolex))
        dipy = np.dot(self.lmat,np.array(dipoley))
        dipz = np.dot(self.lmat,np.array(dipolez))

        self.dipgrad = np.array([np.zeros(3) for i in range(self.ndim)])
        for i in range(self.ndim):
            self.dipgrad[i,0] = dipx[i]
            self.dipgrad[i,1] = dipy[i]
            self.dipgrad[i,2] = dipz[i]


    def getDipoleDerivatives(self):
        return self.dipgrad

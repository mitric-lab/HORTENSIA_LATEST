#!/usr/bin/env python3

import numpy as np
import random

class WignerEnsemble:
    def __init__(self, interface, temp, ncond, fname, dist_type, excmodes):
        self.interface = interface
        self.temp  = temp
        self.ncond = ncond

        self.interface.setStructure()
        self.interface.setMasses()
        self.interface.setHessian()
        self.interface.setMassWeightedHessian()
        self.interface.setTransformationMatrix()

        self.ndim = len(self.interface.getCoordinates())
        self.nat  = self.ndim // 3

        if self.interface.linear == 1:
            self.nzerofreq = 5
        else:
            self.nzerofreq = 6

        self.autoangs  = 0.5291772
        self.fname     = fname
        self.dist_type = dist_type # Wigner or Husimi

        # list of excited vibrational modes,
        # numbering should start at 1 for lowest frequency mode
        # (only for Husimi, currently only v=1 possible)
        self.excmodes = excmodes

        # list of normal modes with quantum numbers 0 or 1
        self.nmode = np.zeros(self.ndim, dtype="int")

        for i in range(len(self.excmodes)):
            self.nmode[self.excmodes[i]+self.nzerofreq-1] = 1


    def initialConditions(self, i):
        A = np.tanh(1.0 * self.interface.getNormalFrequencies()[i] /
                    (2.0 * 3.1671e-6 * self.temp))
        dq = np.sqrt(1.0 / (2*self.interface.getNormalFrequencies()[i]*A))
        dp = np.sqrt(1.0 * self.interface.getNormalFrequencies()[i]/(2*A))
        q0 = random.gauss(0.0, dq)
        p0 = random.gauss(0.0, dp)

        return q0,p0


    def initialConditionsHusimi(self, i):
        def f1(r):
            # arbitrary comparison function,
            # must be always larger than prob distrib
            return 3 * np.exp(-r)

        def f0(r):
            # arbitrary comparison function,
            # must be always larger than prob distrib
            return 2 * np.exp(-r)

        def P1(r):
            # this form comes from int_0^infty r^2 exp(-r^2) r dr *2pi=1
            # for the 2D distribution
            return 2*r**3 * np.exp(-r**2)

        def P0(r):
            # this form comes from int_0^infty exp(-r^2) r dr *2pi=1
            # for the 2D distribution
            return 2*r * np.exp(-r**2)

        flag = True
        while flag:
            r = -np.log(1 - random.random())
            if self.nmode[i] == 1:
                y = f1(r) * random.random()
                p = P1(r)
            elif self.nmode[i] == 0:
                y = f0(r) * random.random()
                p = P0(r)
            if y <= p:
                phi = 2*np.pi * random.random()
                q0 = np.sqrt(2.0/self.interface.getNormalFrequencies()[i]) \
                     * r * np.cos(phi)
                p0 = np.sqrt(2.0*self.interface.getNormalFrequencies()[i]) \
                     * r * np.sin(phi)
                flag = False

        return q0, p0


    def initialConditionsV1(self, i):
        # calculates the probability distribution |psi(q)|^2*|phi(p)|^2 for
        # v=1 of harmonic oscillator

        def f0(r, c):
            # arbitrary comparison function,
            # must be always larger than prob distrib
            return 2*c * np.exp(-r)

        def P0(r, c):
            # probability for v=1
            return c*r**2 * np.exp(-r**2)

        # !! the r values are between 0 and infinity.
        # Negative signs are determined stochastically

        # normalization factor for range 0 --> infinity
        C1 = 4*np.sqrt(self.interface.getNormalFrequencies()[i] / np.pi)
        C2 = 4*np.sqrt(1./self.interface.getNormalFrequencies()[i] / np.pi)

        flag = True
        while flag:
            r = -np.log(1 - random.random())
            y = f0(r, C1) * random.random()
            p = P0(r, C1)
            if y <= p:
                s = random.random()
                if s > 0.5:
                    fac = 1.0
                else:
                    fac = -1.0
                q0 = fac * r * \
                     np.sqrt(1.0/self.interface.getNormalFrequencies()[i])
                flag = False

        flag = True
        while flag:
            r = -np.log(1 - random.random())
            y = f0(r, C2) * random.random()
            p = P0(r, C2)
            if y <= p:
                s = random.random()
                if s > 0.5:
                    fac = 1.0
                else:
                    fac = -1.0
                p0 = fac * r * np.sqrt(self.interface.getNormalFrequencies()[i])
                flag = False

        return q0, p0


    def getWignerEnsemble(self):
        self.invmassmatrix = np.array(
                             [np.zeros(self.ndim) for j in range(self.ndim)])
        for i in range(self.ndim):
            self.invmassmatrix[i,i] = 1/np.sqrt(self.interface.getMasses()[i])

        c0   = self.interface.getCoordinates()
        lmat = np.transpose(self.interface.getTransformationMatrix())

        allout = open("allstruct.xyz","w")
        allout.write(str(self.ncond*self.nat)+"\n")
        allout.write("\n")

        for ii in range(self.ncond):
            q = []
            p = []
            for j in range(self.ndim):
                if j < self.nzerofreq:
                    q0 = 0.0
                    p0 = 0.0
                else:
                    if self.dist_type == "W":
                        q0,p0 = self.initialConditions(j)
                    elif self.dist_type == "H":
                        q0,p0 = self.initialConditionsHusimi(j)
                    elif self.dist_type == "V1":
                        if self.nmode[j] == 0:
                            q0,p0 = self.initialConditions(j)
                        elif self.nmode[j] == 1:
                            q0,p0 = self.initialConditionsV1(j)

                q.append(q0)
                p.append(p0)

            q   = np.array(q)
            p   = np.array(p)
            c   = np.dot(self.invmassmatrix, np.dot(lmat, q))

            vel = np.dot(self.invmassmatrix, np.dot(lmat, p))
            coord = c + c0

            file = open(self.fname+"_"+str(ii).zfill(1)+".in","w")
            file.write(str(self.nat)+"\n")
            for i in range(self.nat):
                file.write("%s %20.14f %20.14f %20.14f\n" %(
                    self.interface.getSymbols()[i],
                    coord[3*i], coord[3*i+1], coord[3*i+2]))
                allout.write("%s %20.14f %20.14f %20.14f\n" %(
                    self.interface.getSymbols()[i],
                    coord[3*i]*self.autoangs, coord[3*i+1]*self.autoangs,
                    coord[3*i+2]*self.autoangs))

            for i in range(self.nat):
                file.write("%20.14f %20.14f %20.14f\n" %(
                    vel[3*i], vel[3*i+1], vel[3*i+2]))
            file.close()
        allout.close()

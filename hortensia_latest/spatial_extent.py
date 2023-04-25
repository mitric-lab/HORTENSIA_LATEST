#!/usr/bin/env python

import os
import subprocess as sp

import numpy as np
from scipy.special import erf
from scipy.optimize import fsolve

import hortensia_latest.calcR2 as calcR2
import hortensia_latest.misc as misc

class readGaussian:
    """
    Gaussian calculation must have been performed with the modified version 
    loaded by module load chem/g09-source
    """
    def __init__(self, qcspecs, fo, doCalc=True):
        if doCalc:
            self.externalQC(qcspecs[0], qcspecs[1], qcspecs[2], fo)
        else:
            fo.write("\nVDE of state %i negative, "%qcspecs[0][1])
            fo.write("starting resonance decay\n")

        with open("anDecay.log", "r") as f:
            self.l = f.readlines()

        state = qcspecs[0][2]
        moflag = 'homo' if state == 0 else 'lumo'

        self.getStructure()
        self.readBasis()
        self.getMOGaussian(moflag)


    def externalQC(self, stateSet, calcSet, molSet, fo):
        """
        Generates input and carries out quantum-chemical calculations with the
        Gaussian program package
        """

        def makeInput(rerun=False):
            """
            Generates input for Gaussian
            """
            # Head
            inp  = "%Chk=anDecay.chk \n%Nproc="+"%i\n"%Nproc
            inp += "%Mem="+"%iMb\n"%(2500*Nproc)
            inp += "#p TD(nstates=%i) %s/GEN "%(max(states), func)
            inp += len(cartBFs)*"%s"%tuple([i for i in cartBFs])
            inp += "Iop(6/7=3) Iop(9/40=14) "
            if rerun:
                inp += "scf(maxcyc=%i,qc,conver=%i) nosym "%(maxcyc, conv)
            else:
                inp += "scf(maxcyc=%i,conver=%i) "%(maxcyc, conv)
                inp += "nosym Guess(Read,TCheck) "
            inp += "GFInput Pop=(minimal,nto,savento) "
            if statenr == 0:
                inp += "Density=(transition=1)"
            else:
                inp += "Density=(transition=%i)"%statenr
            inp += "\n\nComment\n\n%i %i\n"%(charge, mult)

            # Structure
            struct     = ""
            coord_angs = coord * misc.bohr_to_angs
            for i,j in enumerate(coord_angs):
                struct += "%s "%atS[i]
                struct += "%.9f %.9f %.9f\n"%(j[0], j[1], j[2])
            struct += "\n"

            # Tail
            with open("basis","r") as f:
                lines = f.readlines()
            tail = ""
            for i in lines:
                tail += i

            return inp+struct+tail

        states , state   = stateSet[0], stateSet[1]
        statenr, excited = stateSet[2], stateSet[3]

        Nproc  , func    = calcSet[0] , calcSet[1]
        cartBFs, conv    = calcSet[2] , calcSet[3]
        maxcyc , method  = calcSet[4] , calcSet[5].lower()

        coord  , atS     = molSet[0], molSet[1]
        charge , mult    = molSet[2], molSet[3]

        if not excited:
            os.system("cp an.log1 anDecay.log")
            return


        aninp = makeInput()

        with open("anDecay.in", "w") as f:
            f.write(aninp)

        # setting enviromental variable to ensure we calculate in '.'
        cwd = os.getcwd()
        os.environ['GAUSS_SCRDIR'] = cwd

        if os.path.exists("anDecay.log"):
            os.system("rm anDecay.log")

        # actual qc calculation
        fo.write("\nVDE of state %i negative, "%state)
        fo.write("starting resonance decay\n")
        fo.write("Calculating anionic system with NTOs...\n")
        sp.call('%s anDecay.in anDecay.log'%method, shell=True)

        # check for errors in calculation, namely scf convergence and
        # (if excited states are needed) negative excitation energies
        anfail  = False
        cmd     = "grep 'Normal termination' anDecay.log"
        grepout = sp.Popen([cmd], shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
        out,err = grepout.communicate()
        if out == b"":
            anfail = True
            fo.write("anion NTO calculation not converged, " +
                        "restarting calculation\n")
        with open("anDecay.log","r") as f:
            for i in f:
                if "Excited State" in i:
                    exEn1 = float(i.split()[4])
                    if exEn1 < 0:
                        anfail = True
                        fo.write("Negative excitation energy in anion " +
                                "calculation, restarting calculation\n")
                        break

        # if qc calculation failed, calculation is restarted once with the
        # qc scf option and without .chk file of previous calculation
        if anfail:
            fo.write("an.chk removed, switched to QC SCF procedure\n")
            fo.write("Retrying calculation of anion...\n")
            os.system("rm an.chk")

            aninp = makeInput(True)
            with open("anDecay.in","w") as f:
                f.write(aninp)

            sp.call('%s anDecay.in anDecay.log'%method, shell=True)

            with open("anDecay.log","r") as f:
                for i in f:
                    if "Convergence criterion not met" in i:
                        fo.write("anion calculation still not converged, " +
                                    "aborting program\n")
                        quit("SCF for anion not converged!")
                    if "Excited State" in i:
                        exEn1 = float(i.split()[4])
                        if exEn1 < 0:
                            fo.write("Still negative excitation energy in "+
                                        "anion calculation, aborting "+
                                        "program\n")
                            quit("Negative excitation energy in anion!")


    def getMOGaussian(self, moflag):
        for i in range(len(self.l)):
            if "Molecular Orbital Coefficients:" in self.l[i]:
                start = i
                break

        for i in range(start, start+(self.nbas+4)*(self.nbas//5), self.nbas+3):
            if moflag == "homo":
                if "O" in self.l[i+2]:
                    start2 = i+2
            elif moflag == "lumo":
                if "V" in self.l[i+2]:
                    start2 = i+2
                    break

        h = self.l[start2].split()
        for i in range(len(h)):
            if h[i] == "O": 
                ind = i
            if h[i] == "V": 
                if moflag == "lumo": 
                    ind = i
                break

        homo=[]
        count = 0
        for i in range(start2+2, start2+2+self.nbas):
            h = float(self.l[i][21+ind*21:21+ind*21+20])
            orbtype = self.l[i][13:15]
            homo.append(h)
            if orbtype in ["XX", "YY", "ZZ"]: 
                # double the coefficient for quadratic d functions 
                # (necessary for s component in time prop.)
                homo.append(h)
            count += 1

        # return array of AO coefficients for the chosen MO (HOMO or LUMO)
        self.C = np.array(homo)


    def getStructure(self):
        nat = 1
        for i in range(len(self.l)):
            if "NAtoms" in self.l[i]:
                nat = int(self.l[i].split()[1])
            if "Input orientation" in self.l[i]:
                ind = i
                break

        tmpstruct  = ""
        self.atind = []
        self.coord = []
        for i in range(ind+5, ind+5+nat):
            h = self.l[i].split()
            tmpstruct += misc.NrToS(int(h[1])) + " "
            self.atind.append(misc.NrToS(int(h[1])))
            tmp_coord = np.zeros(3)
            for j in range(3, 6):
                tmpstruct += str(float(h[j]) / misc.bohr_to_angs) + " "
                tmp_coord[j-3] = float(h[j]) / misc.bohr_to_angs
            tmpstruct += "; "
            self.coord.append(tmp_coord)
        self.struct = tmpstruct[:-2]
        self.coord  = np.array(self.coord)


    def readBasis(self):
        # basis functions should be provided to Gaussian as sorted by type and 
        # exponent, starting from s orbitals and the largest exponent
        block      = False
        self.basis = []
        self.sym   = []
        for i in range(len(self.l)):
            if "basis functions" in self.l[i]:
                self.nbas = int(self.l[i].split()[0])

            if block and ("****"in self.l[i]): 
                block = False

            if "Centers" in self.l[i]: 
                block = True
                symb  = self.atind[int(self.l[i].split()[1])-1]
                self.sym.append(symb)
                first = True

            if block:
                if "Exponent" in self.l[i]:
                    h = self.l[i].split()
                    if first:
                        orbtype = self.l[i-1].split()[0]
                        first   = False
                        self.basis.append([symb, orbtype])
                        self.basis[-1].append(
                            [[float(h[1].replace("D", "E")), 
                              float(h[3].replace("D", "E"))]])
                    else:
                        self.basis[-1][-1].append(
                            [float(h[1].replace("D", "E")),
                             float(h[3].replace("D", "E"))])
                else:
                    first = True


    def prepR2calc(self):
        nrm        = [] # norms
        al         = [] # exponents
        typel      = [] # orbital types (s,p,d)
        pl         = [] # contraction coefs for each primitive Gaussian
        contl      = [] # contracted Gauss indices for each prim. Gaussian
        Rl         = [] # list of nuclear coordinates
        Ll         = [] # list of l-quantum numbers
        ncont      = 0  # number of contracted basis functions
        icont      = 0  # counting index for contracted Gaussians

        for i in range(len(self.atind)):
            for j in range(len(self.basis)):
                if self.atind[i] != self.basis[j][0]:
                    continue

                if self.basis[j][1] == "S": 
                    lq = [[0,0,0]]
                    ncont = 1
                elif self.basis[j][1] == "P": 
                    lq = [[1,0,0], [0,1,0], [0,0,1]]
                    ncont = 3
                elif self.basis[j][1] == "D": 
                    # Gaussian ordering
                    lq = [[2,0,0], [0,2,0], [0,0,2], 
                        [1,1,0], [1,0,1], [0,1,1]]
                    # 6 d functions + s components for x2 y2 z2
                    ncont = 9

                for k in range(len(self.basis[j][2])):
                    a = self.basis[j][2][k][0] # exponent

                    tmpi = 0

                    for ii in range(len(lq)):
                        nrm.append(misc.normGauss(a, lq[ii][0], lq[ii][1], 
                                                  lq[ii][2]))
                        Rl.append(self.coord[i])
                        Ll.append(lq[ii])
                        al.append(a)
                        typel.append(self.basis[j][1])
                        # contraction coefficient
                        pl.append(self.basis[j][2][k][1])
                        contl.append(icont + tmpi)
                        tmpi += 1

                        if 2 in lq[ii]:
                            # s-type part of the time evolution of 
                            # x**2-type d functions
                            nrm.append(misc.normGauss(a, lq[ii][0], lq[ii][1], 
                                                      lq[ii][2]))
                            Rl.append(self.coord[i])
                            Ll.append([0,0,0])
                            al.append(a)
                            typel.append("DS")
                            # contraction coefficient
                            pl.append(self.basis[j][2][k][1])
                            contl.append(icont + tmpi)
                            tmpi += 1
                icont += ncont

        self.nrm        = np.array(nrm)
        self.al    = np.array(al)
        self.Rl         = np.array(Rl)
        self.Ll         = np.array(Ll)
        self.pl         = np.array(pl)
        self.typel = np.array(typel)
        self.contl      = np.array(contl)

        self.nbas  = len(al)


    def stepR2calc(self, t, prob=0.99, fo=False):
        # list of time-dependent expansion coefs
        bl = np.ones(self.nbas, dtype=complex)

        # update bl and Gaussian exponents at
        at = np.zeros(self.nbas, dtype=complex)
        for i in range(self.nbas):
            f = 1. + 2.j * self.al[i] * t
            at[i] = self.al[i] / f
            bl[i] = 1 / f**1.5
            if self.typel[i] == "P":
                bl[i] /= f
            elif self.typel[i]=="D":
                bl[i] /= f**2 
            elif self.typel[i]=="DS":
                bl[i] *= -1.0j*t / f
        
        r2 = calcR2.getR2(self.nbas, self.C, self.contl, self.pl, self.nrm, 
                          bl, self.Ll, at, self.Rl)

        if t == 0.0:
            self.R0, fvec = self.calcR0(r2, prob)
            fo.write("")
            fo.write("Evaluated R0 value %.3f for probability of %.2f %%\n"%(
                self.R0, prob*100))
            fo.write("Absolute error in P(R0) evaluation: %.5e\n"%fvec)

        # Gaussian width of an 1s function having the calculated <r2> value
        ati = 3./(4.*r2)

        pval = erf(np.sqrt(2*ati) * self.R0) - \
            np.sqrt(8*ati / np.pi) * self.R0 * np.exp(-2*ati * self.R0**2)

        return pval


    def calcR0(self, r2, prob):
        """
        Calculates the distance value R0 for which a 1s Gaussian orbital with 
        <r**2> of HONTO/LUNTO of input molecule has 99 % probability
        """

        def optR0(x, a, prob):
            return erf(np.sqrt(2*a) * x) - \
                np.sqrt(8*a / np.pi) * x * np.exp(-2*a * x**2) - prob
        
        # Alpha value of a 1s Gaussian orbital with same <r**2> value
        a = 3/(4*r2)

        # Calls the fsolve function of scipy.optimize to find the roots of 
        # P(r <= R0) - prob 
        # x0 = 2.0 Angstrom is an arbitrary starting point for R0
        R0, i = fsolve(optR0, x0=2.0, args=(a, prob), full_output=True)[:2]

        return R0[-1], i['fvec'][-1]

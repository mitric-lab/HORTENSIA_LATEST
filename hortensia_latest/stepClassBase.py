#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from configparser import ConfigParser

import numpy as np
import numpy.polynomial.hermite as hp

import hortensia_latest.misc as misc
import hortensia_latest.dysonAux as dyA

config = ConfigParser()
config.read('config.ini')

qcmethod = config['QC']['qcmethod'].lower()
if qcmethod in ["g09", "g16"]:
    import hortensia_latest.dysonG09 as dy
elif qcmethod == "qchem":
    import hortensia_latest.dysonQChem as dy

################################################################################
#                                                                              #
#                 Class for properties of a single time step                   #
#                                                                              #
################################################################################


class Step:
    def __init__(self, states, state, Eshift, grid, an, neu, EA=True):
        """
        parses .log files of quantum chemical calculations
        """

        self.excited = True
        if (len(states) == 1) and (state == 0):
            self.excited = False

        self.anion    = dy.ResultsTDDFT("an" +an, self.excited)
        self.ionized  = dy.ResultsTDDFT("neu"+neu)
        self.ionized.total_energy += Eshift/misc.hartree_to_eV

        if isinstance(EA, bool):
            self.EA = self.ionized.total_energy - self.anion.total_energy
        else:
            self.EA = EA

        self.coord    = self.anion.coordinates.T
        self.energies = self.anion.energies
        self.atnr     = self.anion.atomic_numbers
        self.nat      = self.anion.nat
        self.state    = state
        self.states   = states
        self.nBound   = len(self.states)
        self.grid     = grid
        self.nstates  = self.nBound + self.grid.points
        self.charge   = int(config['QC']['charge'])
        self.mult     = int(config['QC']['mult'])
        self.spin     = self.mult - 1


    def ftDysonPre(self, fo):
        """
        Calculation of Fourier transformed Dyson orbitals, time-independent part
        """

        precision = float(config['States']['precision'])
        max_micro = int(config['States']['max_micro'])
        self.dorbs = dyA.dyson_orbitals(self.ionized, self.anion, fo,
                                        precision, max_micro)
        exp        = self.anion.exponents
        dyson_to0  = self.dorbs[:, 0, :]
        powers     = self.anion.basis.powersD
        newpow, newcoef = self.powerCoefConvert(powers, dyson_to0[:, 0], exp)
        powers     = newpow

        kx  = self.grid.kbasis[:, 0]
        ky  = self.grid.kbasis[:, 1]
        kz  = self.grid.kbasis[:, 2]
        orb = np.zeros((self.anion.basis.nbfsD, self.grid.points),
                       dtype=np.complex)

        ni  = 0
        for i, A in enumerate(self.coord):
            for k, alpha in enumerate(exp[i]):
                c = (1/(2*alpha))**(3/2.)
                if len(powers[i][k][0]) != 1:
                    for l in powers[i][k]:
                        n        = len(l[0]) + len(l[1]) + len(l[2]) - 3
                        orb[ni] += (c * (1j/(2*np.sqrt(alpha)))**n *
                                    self.herm_calc(alpha, kx, ky, kz, l))
                        ni      += 1
                else:
                    l        = powers[i][k]
                    orb[ni] += c * self.herm_calc(alpha, kx, ky, kz, l)
                    ni      += 1

        return orb


    def herm_calc(self, alpha, kx, ky, kz, hermpol):
        """
        evaluates the Hermite polynomials as well as the exp(-k^2/4a)
        term of the Fourier transformed Dyson orbital
        """

        orb = (np.exp(-(kx**2 + ky**2 + kz**2) / (4*alpha)) *
               hp.hermval( ( kx/(2*np.sqrt(alpha)) ), hermpol[0]) *
               hp.hermval( ( ky/(2*np.sqrt(alpha)) ), hermpol[1]) *
               hp.hermval( ( kz/(2*np.sqrt(alpha)) ), hermpol[2]))

        return orb


    def herm_wfn(self, fo, preFT, firstStep=False, refEn=False):
        """
        Calculation of Dyson orbitals and FT of Dyson orbitals from anionic
        states to the neutral ground state
        """

        precision = float(config['States']['precision'])
        max_micro = int(config['States']['max_micro'])
        self.dorbs  = dyA.dyson_orbitals(self.ionized, self.anion, fo,
                                         precision, max_micro)
        self.dorbE  = dyA.dyson_energies(self.ionized, self.anion)
        self.dnorms = dyA.dyson_norms(self.dorbs, self.ionized, self.anion, fo)

        exp         = self.anion.exponents
        powers      = self.anion.basis.powersD
        dyson_to0   = self.dorbs[:, 0, :]

        return_wfn = []
        if firstStep:
            if not config['Dynamics']['restart'] == 'True':
                refEn = 1.0 * self.anion.total_energy

        for i in self.states:
            newpow, newcoef = self.powerCoefConvert(powers, dyson_to0[:, i],
                                                    exp)
            wfn = self.center_orb(exp, newcoef, newpow, preFT)/self.anion.nelec
            return_wfn.append(wfn)

        self.ft_dyson = np.array(return_wfn)
        self.overlap()

        if firstStep:
            return refEn


    def center_orb(self, exp, coef, powers, preFT):
        """
        calculation of Fourier transformed Dyson orbitals
        """

        kx  = self.grid.kbasis[:, 0]
        ky  = self.grid.kbasis[:, 1]
        kz  = self.grid.kbasis[:, 2]
        orb = np.zeros(self.grid.points, dtype=np.complex)

        ni  = 0
        for i, A in enumerate(self.coord):
            for k, alpha in enumerate(exp[i]):
                c = np.exp(1j * (kx*A[0] + ky*A[1] + kz*A[2]))

                if len(powers[i][k][0]) != 1:
                    for m, l in enumerate(powers[i][k]):
                        orb += c * coef[i][k][m] * preFT[ni]
                        ni  += 1
                else:
                    orb += c * coef[i][k] * preFT[ni]
                    ni  += 1

        return orb


    def overlap(self):
        """
        calculates the (relevant portion of the) overlap matrix between the
        bound and free-electron states of the system

        <Psi_(N-1)*k_i|Psi_(N-1)*k_j> = 0 and therefore not calculated
        """

        overlap = np.zeros((self.nBound, self.nstates), dtype=np.complex)

        for i in range(self.nBound):
            overlap[i, self.nBound:] = self.ft_dyson[i] * self.grid.kfact

        self.Smatrix = overlap


    def calcHmatrix(self):
        """
        calculates (parts of) the H matrix of the system in the basis of the
        electronic states of the system, namely the diagonal elements and the
        diabatic couplings between bound and free-electron states
        """

        pass


    def powerCoefConvert(self, powers, coef, exp):
        """
        Convertes the powers and coefficients of the basis functions of a
        molecular orbital from Alexander's formatting to the one needed for the
        calculation of the FT of the MOs
        """

        powDict = {
            (0, 0, 0): [[1], [1], [1]],
            (1, 0, 0): [[0, 1], [1], [1]],
            (0, 1, 0): [[1], [0, 1], [1]],
            (0, 0, 1): [[1], [1], [0, 1]],
            (2, 0, 0): [[0, 0, 1], [1], [1]],
            (0, 2, 0): [[1], [0, 0, 1], [1]],
            (0, 0, 2): [[1], [1], [0, 0, 1]],
            (1, 1, 0): [[0, 1], [0, 1], [1]],
            (1, 0, 1): [[0, 1], [1], [0, 1]],
            (0, 1, 1): [[1], [0, 1], [0, 1]],
            (3, 0, 0): [[0, 0, 0, 1], [1], [1]],
            (0, 3, 0): [[1], [0, 0, 0, 1], [1]],
            (0, 0, 3): [[1], [1], [0, 0, 0, 1]],
            (2, 1, 0): [[0, 0, 1], [0, 1], [1]],
            (2, 0, 1): [[0, 0, 1], [1], [0, 1]],
            (1, 2, 0): [[0, 1], [0, 0, 1], [1]],
            (1, 0, 2): [[0, 1], [1], [0, 0, 1]],
            (1, 1, 1): [[0, 1], [0, 1], [0, 1]],
            (0, 2, 1): [[1], [0, 0, 1], [0, 1]],
            (0, 1, 2): [[1], [0, 1], [0, 0, 1]],
            (4, 0, 0): [[0, 0, 0, 0, 1], [1], [1]],
            (0, 4, 0): [[1], [0, 0, 0, 0, 1], [1]],
            (0, 0, 4): [[1], [1], [0, 0, 0, 0, 1]],
            (3, 1, 0): [[0, 0, 0, 1], [0, 1], [1]],
            (3, 0, 1): [[0, 0, 0, 1], [1], [0, 1]],
            (2, 2, 0): [[0, 0, 1], [0, 0, 1], [1]],
            (2, 0, 2): [[0, 0, 1], [1], [0, 0, 1]],
            (2, 1, 1): [[0, 0, 1], [0, 1], [0, 1]],
            (1, 3, 0): [[0, 1], [0, 0, 0, 1], [1]],
            (1, 0, 3): [[0, 1], [1], [0, 0, 0, 1]],
            (1, 2, 1): [[0, 1], [0, 0, 1], [0, 1]],
            (1, 1, 2): [[0, 1], [0, 1], [0, 0, 1]],
            (0, 3, 1): [[1], [0, 0, 0, 1], [0, 1]],
            (0, 1, 3): [[1], [0, 1], [0, 0, 0, 1]],
            (0, 2, 2): [[1], [0, 0, 1], [0, 0, 1]]}

        newpow, newcoef = [], []
        newpowtemp, newcoeftemp  = [], []
        p_orbs, d_orbs, f_orbs, g_orbs = [], [], [], []
        p_coef, d_coef, f_coef, g_coef = [], [], [], []

        for i in range(len(powers[0])):
            if sum(powers[:, i]) == 0:                          # s
                newpowtemp.append(powDict[tuple(powers[:, i])])
                newcoeftemp.append(coef[i])
            elif sum(powers[:, i]) == 1:                        # p
                p_orbs.append(powDict[tuple(powers[:, i])])
                p_coef.append(coef[i])
                if len(p_orbs) == 3:
                    newpowtemp.append(p_orbs)
                    newcoeftemp.append(p_coef)
                    p_orbs, p_coef = [], []
            elif sum(powers[:, i]) == 2:                        # d
                d_orbs.append(powDict[tuple(powers[:, i])])
                d_coef.append(coef[i])
                if len(d_orbs) == 6:
                    newpowtemp.append(d_orbs)
                    newcoeftemp.append(d_coef)
                    d_orbs, d_coef = [], []
            elif sum(powers[:, i]) == 3:                        # f
                f_orbs.append(powDict[tuple(powers[:, i])])
                f_coef.append(coef[i])
                if len(f_orbs) == 10:
                    newpowtemp.append(f_orbs)
                    newcoeftemp.append(f_coef)
                    f_orbs, f_coef = [], []
            elif sum(powers[:, i]) == 4:                        # g
                g_orbs.append(powDict[tuple(powers[:, i])])
                g_coef.append(coef[i])
                if len(g_orbs) == 15:
                    newpowtemp.append(g_orbs)
                    newcoeftemp.append(g_coef)
                    g_orbs, g_coef = [], []

        n = 0
        for i in exp:
            temp, coeftemp = [], []
            for j in range(len(i)):
                temp.append(newpowtemp[n])
                coeftemp.append(newcoeftemp[n])
                n += 1
            newpow.append(temp)
            newcoef.append(coeftemp)

        return newpow, newcoef

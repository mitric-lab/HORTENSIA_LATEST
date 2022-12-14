#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from configparser import ConfigParser

import numpy as np

import hortensia_latest.dysonAux as dyA
import hortensia_latest.stepClassBase as scBase

config = ConfigParser()
config.read('config.ini')

################################################################################
#                                                                              #
#                 Class for properties of a single time step                   #
#                                                                              #
################################################################################


class Step(scBase.Step):
    def __init__(self, states, state, Eshift, grid, an, neu, EA=True):
        super().__init__(states, state, Eshift, grid, an, neu, EA)


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

        self.orthoRenorm()

        if firstStep:
            return refEn


    def orthoRenorm(self):
        """
        Calculates the normalization constant that arises due to the different
        orthogonalization procedures
        """

        self.renormD = np.ones(self.nstates - self.nBound, dtype=complex)

        ftDsq = np.zeros(self.grid.points)
        for i in range(self.nBound):
            ftDsq += np.abs(self.Smatrix[i, self.nBound:])**2

        Sdy = dyA.dyson_overlap_K(self, self, self.nBound, self.states)

        # Dyson orbitals not orthogonal, so normalization constant is
        # complex...
        ftD2  = np.zeros(self.grid.points, dtype=complex)
        for i in range(self.nBound):
            for j in range(self.nBound):
                ftD2 += np.conjugate(self.Smatrix[i, self.nBound:]) * \
                        Sdy[i, j] * self.Smatrix[j, self.nBound:]

        self.renormD /= np.sqrt(1 - 2*ftDsq + ftD2)


    def calcHmatrix(self, refEn, exInt=False):
        """
        calculates (parts of) the H matrix of the system in the basis of the
        electronic states of the system, namely the diagonal elements and the
        diabatic couplings between bound and free-electron states
        """

        en  = np.zeros(self.nstates)
        dia = np.zeros((self.nBound, self.nstates), dtype=np.complex)
        diaTemp = np.zeros((self.nBound, self.nstates), dtype=np.complex)

        ### diagonal bound state part of H-Matrix
        for i, j in enumerate(self.states):
            en[i] = self.energies[j] - refEn

        # Phi_ion(k_i) -> E_ion + k_i^2/2
        # NOT FINISHED, JUST FOR TESTING !!!!!
        en[self.nBound:] = self.ionized.total_energy - refEn + \
                            self.grid.kinEn
        self.Hdiag = en

        self.diaTemp = diaTemp
        self.diacoup = dia

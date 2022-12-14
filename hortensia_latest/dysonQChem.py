#!/usr/bin/env python3

"""
compute Dyson orbitals from two TD-DFT calculations, one for a molecule with
N electrons the other for a molecule with N-1 electrons.
"""

import numpy as np
import numpy.linalg as la
import logging

import hortensia_latest.dysonAux as dyA

########################################################################
#
# PARSING GAUSSIAN 09 OUTPUT
#
########################################################################

class UncontractedBasisSet:
    def __init__(self, nbfs):
        """
        initialize basis for `nbfs` primitive Gaussian functions
        """

        self.nbfs = nbfs
        self.exponents = np.zeros(nbfs,     dtype=float)
        self.powers    = np.zeros((3,nbfs), dtype=int)
        self.centers   = np.zeros((3,nbfs), dtype=float)
        # index of atomic center to which the basis function belongs
        self.atom_ids  = np.zeros(nbfs,     dtype=int)
        # index of angular momentum shell
        self.shell_index = np.zeros(nbfs,   dtype=int)


class ParserError(Exception):
    pass


class ResultsTDDFT:
    """
    Class for parsing QC calculations from QChem output files
    """

    def __init__(self, log_file, excited=False, check_results=True):
        """parse QChem output file for DFT calculation"""
        logging.info("reading data from QChem output '%s'" % log_file)
        self.excited  = excited
        self.log_file = "%s.log%s" % (log_file[:-1], log_file[-1])
        self.fchk     = "%s.fchk%s" % (log_file[:-1], log_file[-1])
        with open(self.log_file) as f:
            llog  = f.readlines()
            if not "cleanup process" in llog[-1]:
                raise ParserError("No geometry found in '%s'" % log_file)
        with open(self.fchk) as f:
            lfchk = f.readlines()

        self.readGeometry(lfchk)
        self.readExponents(lfchk)
        self.setupBasis(lfchk)
        self.readTotalEnergy(llog)
        self.readExcitedStates(lfchk)
        self.readMolecularOrbitals(lfchk)

        if check_results == True:
            self.computeCIScoefficients()
            self.printResults()


    def readGeometry(self, lines):
        """
        read molecular geometry in standard or input orientation
        """

        nl = len(lines)
        # find number of atoms
        for i in range(nl):
            if "Number of atoms" in lines[i]:
                self.nat = int(lines[i].split()[4])

        self.atomic_numbers = np.zeros(self.nat,  dtype=int)
        self.coordinates = np.zeros((3,self.nat), dtype=float)

        for i in range(nl):
            if "Atomic numbers" in lines[i]:
                for j in range(self.nat):
                    self.atomic_numbers[j] = lines[i+1].split()[j]
            if "Current cartesian" in lines[i]:
                j = 0
                u,v = 0,0
                while j < self.nat*3:
                    i += 1
                    temp = lines[i].split()
                    for k in temp:
                        self.coordinates[u,v] = float(k)
                        u = (u+1)%3
                        if u == 0:
                            v += 1
                        j += 1

        logging.info("number of atoms: %s" % self.nat)


    def readExponents(self, lines):
        """
        read exponents for Gaussian basis functions on each atom
        """

        nl = len(lines)

        # find basis section
        for i in range(nl):
            if "Shell to atom map" in lines[i]:
                nshells     = int(lines[i].split()[6])
                shellToAtom = np.zeros(nshells, dtype=int)

                j  = 0
                while j < nshells:
                    i += 1
                    temp = lines[i].split()
                    for k, tk in enumerate(temp):
                        shellToAtom[j+k] = int(tk)-1
                    j += len(temp)

                i += 1
                primExp = np.zeros(nshells)

                j  = 0
                while j < nshells:
                    i += 1
                    temp = lines[i].split()
                    for k, tk in enumerate(temp):
                        primExp[j+k] = float(tk)
                    j += len(temp)

                i += 1
                j  = 0
                while j < nshells:
                    i += 1
                    temp = lines[i].split()
                    for k, tk in enumerate(temp):
                        if float(tk)-1.0 > 1e-8:
                            text  = "Only fully uncontracted basis functions "
                            text += "(with all coefficients = 1.0) are "
                            text += "supported. Please decontract the basis!"
                            raise ParserError(text)
                    j += len(temp)

                break

        self.atom_ids = shellToAtom

        self.exponents = []
        for j in range(self.nat):
            self.exponents.append([])

            for k in range(nshells):
                if shellToAtom[k] == j:
                    self.exponents[j].append(primExp[k])


    def setupBasis(self, lines):
        """
        Generation of UncontractedBasisSet instance for storage of basis sets
        """

        nl = len(lines)

        # read number of electrons
        for i in range(nl):
            if "Number of electrons" in lines[i]:
                self.nelec       = int(lines[i].split()[4])
                self.nelec_alpha = int(lines[i+1].split()[5])
                self.nelec_beta  = int(lines[i+2].split()[5])
                break

        # read number of MOs
        for j in range(i, nl):
            if "Number of independent functions" in lines[j]:
                self.nmo = int(lines[j].split()[5])
                break

        self.basis = UncontractedBasisSet(self.nmo)


    def readTotalEnergy(self, lines):
        """
        Reading of total energy of molecule
        """

        nl = len(lines)
        for i in range(nl):
            if "Total energy in" in lines[i]:
                self.total_energy = float(lines[i].split()[8])
                break
        else:
            raise ParserError("No final SCF energy could be found!")


    def readExcitedStates(self, lines):
        """
        Reading of excitation coefficients from output
        """

        if not self.excited:
            # If only ground state is used in dynamics, skips excited state
            # evaluation and create single microstate, the ground state
            # determinant
            na = self.nelec_alpha
            nb = self.nelec_beta
            norb = self.nmo
            self.norb = norb
            self.microstates = list(dyA.enumerate_microstates(na, nb, norb))
            nstates = 1
            # TD-DFT excitation coefficients
            X = np.array([[1.0]])
            # TD-DFT deexcitation coefficients
            Y = np.array([[0.0]])
            self.excitation_energies = np.array([0.0])
            self.oscillator_strengths = np.array([0.0])
            state = 0
            nstates = state + 1
            self.nstates = nstates
            self.XpY = X + Y
            self.XmY = X - Y

            # compute total energies of states
            self.energies = np.zeros(nstates, dtype=float)
            self.energies = self.total_energy + \
                            self.excitation_energies[:nstates]

            self.nmicro = 1

            return lines

        nl = len(lines)

        # Number of excited states
        for i in range(nl):
            if "Number of Excited States" in lines[i]:
                nstates = int(lines[i].split()[5]) + 1
                break
        else:
            raise ParserError("Could not find 'roots to seek' in logfile")

        na = self.nelec_alpha
        nb = self.nelec_beta
        norb = self.nmo
        self.norb = norb
        self.microstates = list(dyA.enumerate_microstates(na,nb,norb))
        # Number of microstates
        nmicro = 1 + na*(norb-na) + nb*(norb-nb)
        assert nmicro == len(self.microstates)
        self.nmicro = nmicro
        logging.info("number of micro states: %s" % nmicro)
        logging.info("number of excited states: %s" % (nstates-1))

        # TD-DFT excitation coefficients
        X = np.zeros((nmicro, nstates), dtype=float)
        # TD-DFT deexcitation coefficients
        Y = np.zeros((nmicro, nstates), dtype=float)
        # Extended versions of X and Y with more info for easier recognition
        # of microstates
        Xk = np.zeros((nmicro, nstates, 4))
        Xk[0,0,0] = 1.0
        Yk = np.zeros((nmicro, nstates, 4))
        Yk[0,0,0] = 0.0
        XAB = np.zeros(nmicro, dtype=str)
        YAB = np.zeros(nmicro, dtype=str)

        self.excitation_energies = np.zeros(nstates, dtype=float)
        self.oscillator_strengths = np.zeros(nstates, dtype=float)
        # Ground state
        X[0,0] = 1.0
        Y[0,0] = 0.0

        for i in range(i, nl):
            if "Excitation Energies" in lines[i]:
                for j in range(1, nstates):
                    self.excitation_energies[j]  = lines[i+1].split()[j-1]
                    self.oscillator_strengths[j] = lines[i+3].split()[j-1]
                break

        # Alpha excitations
        for i in range(i, nl):
            if "Alpha X Amplitudes" in lines[i]:
                nocc, nvirt = na, norb - na
                if not nocc * nvirt * (nstates-1) == int(lines[i].split()[5]):
                    raise ParserError("Something wrong with nocc/nvirt " +
                                      "in Alpha Amplitudes")
                nampl = int(lines[i].split()[5])
                state = 1
                occ   = 0
                virt  = 0

                n = 0
                while True:
                    i += 1
                    if n >= nampl:
                        break
                    temp = lines[i].split()
                    for j, tj in enumerate(temp):
                        coef  = tj
                        index = 1 + (norb-na)*occ + virt
                        X[index, state]     = coef
                        Xk[index, state, 0] = coef
                        Xk[index, state, 1] = occ
                        Xk[index, state, 2] = virt + na
                        Xk[index, state, 3] = index
                        XAB[index] = "A"
                        virt += 1
                        if virt%(norb - na) == 0:
                            occ += 1
                            virt = 0
                            if occ%nocc == 0:
                                state += 1
                                occ    = 0
                    n += len(temp)
                    if n >= nampl:
                        break
                i += 1
                break

        # Beta excitations
        for i in range(i, nl):
            if "Beta X Amplitudes" in lines[i]:
                nocc, nvirt = nb, norb - nb
                if not nocc * nvirt * (nstates - 1) == int(lines[i].split()[5]):
                    raise ParserError("Something wrong with nocc/nvirt " +
                                      "in Beta Amplitudes")
                nampl = int(lines[i].split()[5])
                state = 1
                occ   = 0
                virt  = 0

                n = 0
                while True:
                    i += 1
                    if n >= nampl:
                        break
                    temp = lines[i].split()
                    for j, tj in enumerate(temp):
                        coef  = tj
                        index = 1 + (norb - nb) * occ + virt + na * (norb - na)
                        X[index, state]     = coef
                        Xk[index, state, 0] = coef
                        Xk[index, state, 1] = occ
                        Xk[index, state, 2] = virt + nb
                        Xk[index, state, 3] = index
                        XAB[index] = "B"
                        virt += 1
                        if virt%(norb - nb) == 0:
                            occ += 1
                            virt = 0
                            if occ%nocc == 0:
                                state += 1
                                occ    = 0
                    n += len(temp)
                    if n >= nampl:
                        break
                i += 1
                break

        # Alpha deexcitations
        for i in range(i, nl):
            if "Alpha Y Amplitudes" in lines[i]:
                nocc, nvirt = na, norb - na
                if not nocc * nvirt * (nstates-1) == int(lines[i].split()[5]):
                    raise ParserError("Something wrong with nocc/nvirt " +
                                      "in Alpha Amplitudes")
                nampl = int(lines[i].split()[5])
                state = 1
                occ   = 0
                virt  = 0

                n = 0
                while True:
                    i += 1
                    if n >= nampl:
                        break
                    temp = lines[i].split()
                    for j, tj in enumerate(temp):
                        coef  = tj
                        index = 1 + (norb-na)*occ + virt
                        Y[index, state]     = coef
                        Yk[index, state, 0] = coef
                        Yk[index, state, 1] = occ
                        Yk[index, state, 2] = virt + na
                        Yk[index, state, 3] = index
                        YAB[index] = "A"
                        virt += 1
                        if virt%(norb - na) == 0:
                            occ += 1
                            virt = 0
                            if occ%nocc == 0:
                                state += 1
                                occ    = 0
                    n += len(temp)
                    if n >= nampl:
                        break
                i += 1
                break

        # Beta deexcitations
        for i in range(i, nl):
            if "Beta Y Amplitudes" in lines[i]:
                nocc, nvirt = nb, norb - nb
                if not nocc * nvirt * (nstates - 1) == int(lines[i].split()[5]):
                    raise ParserError("Something wrong with nocc/nvirt " +
                                      "in Beta Amplitudes")
                nampl = int(lines[i].split()[5])
                state = 1
                occ   = 0
                virt  = 0

                n = 0
                while True:
                    i += 1
                    if n >= nampl:
                        break
                    temp = lines[i].split()
                    for j, tj in enumerate(temp):
                        coef  = tj
                        index = 1 + (norb - nb) * occ + virt + na * (norb - na)
                        Y[index, state]     = coef
                        Yk[index, state, 0] = coef
                        Yk[index, state, 1] = occ
                        Yk[index, state, 2] = virt + nb
                        Yk[index, state, 3] = index
                        YAB[index] = "B"
                        virt += 1
                        if virt%(norb - nb) == 0:
                            occ += 1
                            virt = 0
                            if occ%nocc == 0:
                                state += 1
                                occ    = 0
                    n += len(temp)
                    if n >= nampl:
                        break

        self.Xk   = Xk[:, :, :]
        self.Yk   = Yk[:, :, :]
        self.XpYk = Xk[:, :, :] * 1.0
        self.XpYk[:, :, 0] += self.Yk[:, :, 0]
        self.XAB  = XAB
        self.YAB  = YAB

        self.nstates = nstates
        self.XpY = X+Y
        self.XmY = X-Y

        err = la.norm( np.dot(self.XpY.transpose(), self.XmY) - np.eye(nstates) )
        # check normalization
        logging.info("|<X+Y|X-Y> - Id|= %e" % err)

        # compute total energies of states
        self.energies = np.zeros(nstates, dtype=float)
        self.energies = self.total_energy + self.excitation_energies[:nstates]


    def readMolecularOrbitals(self, lines):
        """
        Parsing of AO coefficients of MOs
        """

        nl   = len(lines)
        nbfs = self.basis.nbfs
        nmo  = self.nmo

        for i in range(nl):
            if "Shell types" in lines[i]:
                nshell = int(lines[i].split()[4])
                shelltype = np.zeros(nbfs, dtype=int)

                nid = 0
                n   = 0
                m   = 1
                while True:
                    i += 1
                    parts = lines[i].split()
                    for k,j in enumerate(parts):
                        if self.atom_ids[k] != nid:
                            m = 1
                        if int(j) == 0:
                            shelltype[n]              = 0
                            self.basis.shell_index[n] = m
                            self.basis.powers[:, n]   = (0, 0, 0)
                            n += 1
                        elif int(j) == 1:
                            shelltype[n:n+3]              = 1
                            self.basis.shell_index[n:n+3] = m
                            self.basis.powers[:, n]       = (1, 0, 0)
                            self.basis.powers[:, n + 1]   = (0, 1, 0)
                            self.basis.powers[:, n + 2]   = (0, 0, 1)
                            n += 3
                        elif int(j) == 2:
                            shelltype[n:n+6]              = 2
                            self.basis.shell_index[n:n+6] = m
                            self.basis.powers[:, n]     = (2, 0, 0)
                            self.basis.powers[:, n + 1] = (0, 2, 0)
                            self.basis.powers[:, n + 2] = (0, 0, 2)
                            self.basis.powers[:, n + 3] = (1, 1, 0)
                            self.basis.powers[:, n + 4] = (1, 0, 1)
                            self.basis.powers[:, n + 5] = (0, 1, 1)
                            n += 6
                        elif int(j) == 3:
                            shelltype[n:n+10]              = 3
                            self.basis.shell_index[n:n+10] = m
                            self.basis.powers[:, n]     = (3, 0, 0)
                            self.basis.powers[:, n + 1] = (0, 3, 0)
                            self.basis.powers[:, n + 2] = (0, 0, 3)
                            self.basis.powers[:, n + 3] = (1, 2, 0)
                            self.basis.powers[:, n + 4] = (2, 1, 0)
                            self.basis.powers[:, n + 5] = (2, 0, 1)
                            self.basis.powers[:, n + 6] = (1, 0, 2)
                            self.basis.powers[:, n + 7] = (0, 1, 2)
                            self.basis.powers[:, n + 8] = (0, 2, 1)
                            self.basis.powers[:, n + 9] = (1, 1, 1)
                            n += 10
                        else:
                            raise ParserError("Atomic orbitals above F " +
                                              "not supported!")
                        m += 1

                    if n >= nshell:
                        break
        self.basis.shelltype = shelltype

        # Alpha MOs
        orbs = np.zeros((nbfs, nmo), dtype=float)
        for i in range(nl):
            if "Alpha MO coefficients" in lines[i]:
                n,m = 0,0
                while True:
                    i += 1
                    parts = lines[i].split()
                    for j, pj in enumerate(parts):
                        orbs[m,n] = pj
                        if (n + 1) * (m + 1) >= nmo * nbfs:
                            break
                        m += 1
                        if m%nmo == 0:
                            m  = 0
                            n += 1
                    if (n+1)*(m+1) >= nmo*nbfs:
                        break

        self.orbs_alpha = orbs

        # Beta MOs
        orbs = np.zeros((nbfs, nmo), dtype=float)
        for i in range(nl):
            if "Beta MO coefficients" in lines[i]:
                n,m = 0,0
                while True:
                    i += 1
                    parts = lines[i].split()
                    for j, pj in enumerate(parts):
                        orbs[m,n] = pj
                        if (n + 1) * (m + 1) >= nmo * nbfs:
                            break
                        m += 1
                        if m%nmo == 0:
                            m  = 0
                            n += 1
                    if (n+1)*(m+1) >= nmo*nbfs:
                        break
        self.orbs_beta = orbs

        orbe = np.zeros(nmo, dtype=float)
        for i in range(nl):
            if "Alpha Orbital Energies" in lines[i]:
                n = 0
                while True:
                    i += 1
                    parts = lines[i].split()
                    for j, pj in enumerate(parts):
                        orbe[n+j] = pj
                    n += len(parts)
                    if n >= nmo:
                        break
        self.orbe_alpha = orbe

        orbe = np.zeros(nmo, dtype=float)
        for i in range(nl):
            if "Beta Orbital Energies" in lines[i]:
                n = 0
                while True:
                    i += 1
                    parts = lines[i].split()
                    for j, pj in enumerate(parts):
                        orbe[n + j] = pj
                    n += len(parts)
                    if n >= nmo:
                        break
        self.orbe_beta = orbe

        n = 0
        for i, ids in enumerate(self.atom_ids):
            if shelltype[n] == 0:
                self.basis.atom_ids[n] = ids
                n += 1
            else:
                for j in range((shelltype[i] + 1) * (shelltype[i] + 2) // 2):
                    self.basis.atom_ids[n] = ids
                    n += 1

        for k in range(self.basis.nbfs):
            self.basis.centers[:, k] = self.coordinates[:,
                                                        self.basis.atom_ids[k]]

        center = 0
        nAO = 0
        n = 0
        while n < self.basis.nbfs:
            if center != self.basis.atom_ids[n]:
                nAO = 0
                center = self.basis.atom_ids[n]

            if shelltype[n] == 0:
                self.basis.exponents[n] = self.exponents[center][nAO]
                nAO += 1
                n += 1
            else:
                sindex = shelltype[n]
                for i in range((sindex+1)*(sindex+2)//2):
                    self.basis.exponents[n] = self.exponents[center][nAO]
                    n += 1
                nAO += 1

        return


    def computeCIScoefficients(self):
        Sao = dyA.basis_overlap(self.basis, self.basis)
        self.Sao = Sao.copy()

        self.cis = np.zeros(self.XpY.shape)
        # Ground state
        self.cis[0, 0] = 1.0

        if not self.excited:
            return

        # A-B
        AmB = dyA.AminusB_diagonal(self.orbe_alpha, self.orbe_beta,
                                   self.nelec_alpha, self.nelec_beta, self.norb)

        ### Kevin 07.02.20
        self.cisk = np.zeros(self.XpYk.shape)
        self.cisk[0, 0, 0] = 1.0
        ###

        # excited states
        for i in range(1, self.nstates):
            omega = self.excitation_energies[i]
            self.cis[:, i] = np.sqrt(omega) * self.XpY[:, i] / np.sqrt(AmB)
            ### Kevin 07.02.20
            self.cisk[:, i, 0]  = np.sqrt(omega) * self.XpYk[:, i, 0] / \
                                  np.sqrt(AmB)
            self.cisk[:, i, 1:] = self.XpYk[:, i, 1:].copy()

        logging.info(np.dot(self.cis.transpose(), self.cis))


    def printResults(self):
        logging.info("%d alpha electrons  %d beta electrons" %
                        (self.nelec_alpha, self.nelec_beta))
        logging.info("%d basis functions" % self.basis.nbfs)
        logging.info("Total SCF energy: %s Hartree\n" % self.total_energy)
        logging.info("number of electronic states: %d" % self.nstates)

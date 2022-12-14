#!/usr/bin/env python3

"""
compute Dyson orbitals from two TD-DFT calculations, one for a molecule with
N electrons the other for a molecule with N-1 electrons.

Requirements for Gaussian jobs:

- The source code of Gaussian 09 has to be modified to write out M.O.
  coefficients, TD-DFT eigenvectors and energies and oscillator strengths with
  higher precision. The modified Fortran file can be found in 'g09-source/g09'.
"""

import numpy as np
import numpy.linalg as la
import logging

import hortensia_latest.misc as misc
import hortensia_latest.dysonAux as dyA

################################################################################
#                                                                              #
#                        Parsing Gaussian09 Output                             #
#                                                                              #
################################################################################


class UncontractedBasisSet:
    def __init__(self, nbfsC, nbfsD):
        """
        initialize basis for `nbfs` primitive Gaussian functions
        C indicates contracted basis functions, D decontracted
        """

        self.nbfsC        = nbfsC
        self.exponentsC   = np.zeros(nbfsC,      dtype=float)
        self.exponentsC   = self.exponentsC.tolist()
        self.powersC      = np.zeros((3, nbfsC), dtype=int  )
        self.centersC     = np.zeros((3, nbfsC), dtype=float)
        self.nbfsD        = nbfsD
        self.exponentsD   = np.zeros(nbfsD,      dtype=float)
        self.powersD      = np.zeros((3, nbfsD), dtype=int  )
        self.centersD     = np.zeros((3, nbfsD), dtype=float)
        # index of atomic center to which the basis function belongs
        self.atom_idsC    = np.zeros(nbfsC, dtype=int)
        self.atom_idsD    = np.zeros(nbfsD, dtype=int)
        # index of angular momentum shell
        self.shell_indexC = np.zeros(nbfsC, dtype=int)
        self.shell_indexD = np.zeros(nbfsD, dtype=int)


    def decontractBasis(self):
        """
        Takes the contracted basis set and decontracts it
        """

        self.Omatrix = np.zeros((self.nbfsC, self.nbfsD))
        self.Pmatrix = np.zeros((self.nbfsC, self.nbfsD), dtype=int)

        nc, nd = 0, 0
        while nc < self.nbfsC:
            l = len(self.exponentsC[nc])
            psum  = sum(self.powersC[:, nc])
            nfunc = (psum + 1) * (psum + 2) // 2

            tmp = 0
            for j in range(l):
                for k in range(nfunc):
                    ni, nj = nc+tmp%nfunc, nd+tmp
                    self.exponentsD[nj]  = self.exponentsC[ni][j][0]
                    self.Omatrix[ni, nj] = self.exponentsC[ni][j][1]
                    self.Pmatrix[ni, nj] = 1.0
                    tmp += 1
            nc += nfunc
            nd += nfunc * l

        self.powersD   = np.dot(self.powersC, self.Pmatrix)
        self.centersD  = np.dot(self.centersC, self.Pmatrix)
        self.atom_idsD = np.dot(self.atom_idsC, self.Pmatrix)

        na, nsh = 0, 1
        i = 0
        while i < self.nbfsD:
            if na != self.atom_idsD[i]:
                nsh = 1
            na = self.atom_idsD[i]

            l = sum(self.powersD[:, i])
            for _ in range((l+1)*(l+2)//2):
                self.shell_indexD[i] = nsh
                i += 1
            nsh += 1

        self.S   = dyA.basis_overlap(self, self)
        l, m, n  = self.powersD[0], self.powersD[1], self.powersD[2]
        alpha    = self.exponentsD
        self.N   = dyA.norm(alpha, l, m, n)

        return self.Omatrix


class ParserError(Exception):
    pass


class ResultsTDDFT:
    """
    Class for parsing QC calculations from Gaussian output files
    """
    def __init__(self, log_file, excited=False, check_results=True):
        """
        parse Gaussian log file for TD-DFT calculation
        """

        logging.info("reading data from Gaussian log-file '%s'\n" % log_file)
        log_file = "%s.log%s" % (log_file[:-1], log_file[-1])
        self.log_file = log_file
        self.excited  = excited
        fh = open(log_file)
        lines = fh.readlines()
        lines = self.readGeometry(lines)
        lines = self.readExponents(lines)
        lines = self.setupBasis(lines)
        lines = self.readTotalEnergy(lines)
        lines = self.readExcitedStates(lines)
        if self.nelec_alpha != self.nelec_beta:
            lines = self.readMolecularOrbitals(lines)
        lines = self.readMolecularOrbitals(lines, True)
        fh.close()

        if check_results == True:
            self.computeCIScoefficients()
            self.printResults()


    def readGeometry(self, lines):
        """
        read molecular geometry in standard or input orientation
        """

        nl = len(lines)
        # find number of atoms
        for i in range(0, nl):
            if "NAtoms=" in lines[i]:
                parts = lines[i].split()
                self.nat = int(parts[1])
        # try to find standard orientation
        for i in range(0, nl):
            if "Standard orientation:" in lines[i]:
                break
        else:
            # try to find input orientation (in case 'NoSym' was specified)
            for i in range(0, nl):
                if "Input orientation:" in lines[i]:
                    break
            else:
                raise ParserError("No geometry found in '%s'" % self.log_file)
        istart = i + 5
        self.atomic_numbers = np.zeros(self.nat, dtype=int)
        self.coordinates = np.zeros((3, self.nat), dtype=float)
        logging.info("number of atoms: %s" % self.nat)
        for j in range(0, self.nat):
            i = istart+j
            if "----" in lines[i]:
                break
            parts = lines[i].split()
            self.atomic_numbers[j] = int(parts[1])
            self.coordinates[:, j] = parts[3:]
        # convert from Angstroms to bohr
        self.coordinates /= misc.bohr_to_angs

        return lines[i:]


    def readExponents(self, lines):
        """
        read exponents for Gaussian basis functions on each atom
        """

        nl = len(lines)
        # check that cartesian basis functions are used
        for i in range(0, nl):
            l = lines[i]
            if "General basis read from cards:" in l:
                if not "(6D, 10F)" in l:
                    raise ParserError(
                        "Only cartesian basis functions supported. Please add"+
                        " the keywords '6D' and '10F' to the route section!")
            break
        # find basis section
        for i in range(i, nl):
            if "AO basis set in the form of general basis input" in lines[i]:
                break

        self.exponents = []
        for j in range(self.nat):
            self.exponents.append([])
            i += 2
            nb = 0
            while not "****" in lines[i]:
                self.exponents[j].append([])
                ncont = int(lines[i].split()[1])
                for k in range(ncont):
                    i += 1
                    tmpexp = float(lines[i].split()[0].replace(
                             "D", "E").replace("d", "E"))
                    tmpfac = float(lines[i].split()[1].replace(
                             "D", "E").replace("d", "E"))
                    self.exponents[j][nb].append([tmpexp, tmpfac])
                i  += 1
                nb += 1

        return lines[i:]


    def setupBasis(self, lines):
        """
        Generation of UncontractedBasisSet instance for storage of basis sets
        """

        nl = len(lines)
        for i in range(0, nl):
            l = lines[i]
            if ("basis functions" in l) and ("primitive gaussians" in l):
                nbfs_contr = int(l.split()[0])
                nbfs_decon = int(l.split()[3])
                break
        i += 1
        # read number of alpha and beta electrons
        parts = lines[i].split()
        self.nelec_alpha = int(parts[0])
        self.nelec_beta  = int(parts[3])
        self.nelec = self.nelec_alpha + self.nelec_beta

        # read number of M.O.'s
        for i in range(i, nl):
            l = lines[i]
            if "NBsUse=" in l:
                parts = l.split()
                self.nmo = int(parts[1])
                break

        self.basis = UncontractedBasisSet(nbfs_contr, nbfs_decon)

        return lines[i:]


    def readTotalEnergy(self, lines):
        """
        Reading of total energy of molecule
        """

        nl = len(lines)
        for i in range(0, nl):
            if "SCF Done:" in lines[i]:
                parts = lines[i].split()
                self.total_energy = float(parts[4])
                break
        else:
            raise ParserError("No final SCF energy could be found!")
        return lines[i:]


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
        # number of correlated orbitals
        for i in range(0, nl):
            l = lines[i]
            if "Range of M.O.s used for correlation:" in l:
                parts = l.split()
                orb_corr_first, orb_corr_last = int(parts[-2])-1, int(parts[-1])
                if (orb_corr_first != 0):
                    raise ParserError("All electrons have to be correlated. "+
                        "Add the option 'Full' to the keyword 'TD'!")
                break
        # number of excited states
        for i in range(i, nl):
            l = lines[i]
            if "roots to seek" in l:
                parts = l.split()
                nstates = int(parts[6]) + 1
                break
        else:
            raise ParserError("Could not find 'roots to seek' in logfile")

        na = self.nelec_alpha
        nb = self.nelec_beta
        norb = orb_corr_last - orb_corr_first
        self.norb = norb
        self.microstates = list(dyA.enumerate_microstates(na, nb, norb))
        # Number of microstates
        nmicro = 1 + na*(norb-na) + nb*(norb-nb)
        assert nmicro == len(self.microstates)
        self.nmicro = nmicro
        logging.info("number of micro states: %s" % nmicro)
        logging.info("number of excited states: %s" % (nstates - 1))

        # TD-DFT excitation coefficients
        X = np.zeros((nmicro, nstates), dtype=float)
        # TD-DFT deexcitation coefficients
        Y = np.zeros((nmicro, nstates), dtype=float)
        # Extended versions of X and Y with more info for easier recognition
        # of microstates
        Xk = np.zeros((nmicro, nstates, 4))
        Xk[0, 0, 0] = 1.0
        Yk = np.zeros((nmicro, nstates, 4))
        Yk[0, 0, 0] = 0.0
        XAB = np.zeros(nmicro, dtype=str)
        YAB = np.zeros(nmicro, dtype=str)

        self.excitation_energies = np.zeros(nstates, dtype=float)
        self.oscillator_strengths = np.zeros(nstates, dtype=float)
        # Ground state
        X[0, 0] = 1.0
        Y[0, 0] = 0.0

        state = 0
        for i in range(i, nl):
            l = lines[i]
            if "Excited State" in l:
                state += 1
                parts = l.split()
                en = float(parts[4]) / misc.hartree_to_eV
                f = float(parts[8].replace("f=", ""))
                self.excitation_energies[state] = en
                self.oscillator_strengths[state] = f
            elif "->" in l:
                # Excitation coefficients
                parts = l.split()
                if len(parts) == 4:
                    occ, virt, coef = parts[0], parts[2], float(parts[3])
                elif len(parts) == 3:
                    occ, virt, coef = parts[0], parts[1][2:], float(parts[2])
                if "A" in occ:
                    # Alpha excitation
                    occ = int(occ.replace("A", "")) - 1
                    virt = int(virt.replace("A", "")) - na - 1
                    index = 1 + (norb-na)*occ + virt
                    X[index, state] = coef
                    Xk[index, state, 0] = coef
                    Xk[index, state, 1] = occ
                    Xk[index, state, 2] = virt + na
                    Xk[index, state, 3] = index
                    XAB[index] = "A"
                elif "B" in occ:
                    # Beta excitation
                    occ = int(occ.replace("B", "")) - 1
                    virt = int(virt.replace("B", "")) - nb - 1
                    index = 1 + (norb-nb)*occ + virt + na*(norb-na)
                    X[index, state] = coef
                    Xk[index, state, 0] = coef
                    Xk[index, state, 1] = occ
                    Xk[index, state, 2] = virt + nb
                    Xk[index, state, 3] = index
                    XAB[index] = "B"
                else:
                    # Excitation in closed-shell system
                    occ = int(occ) - 1
                    virt = int(virt) - na - 1
                    assert na == nb
                    indexA = 1 + (norb-na)*occ + virt
                    indexB = indexA + na*(norb-na)
                    X[indexA, state] = coef
                    X[indexB, state] = coef
                    Xk[indexA, state, 0] = coef
                    Xk[indexA, state, 1] = occ
                    Xk[indexA, state, 2] = virt + na
                    Xk[indexA, state, 3] = indexA
                    Xk[indexB, state, 0] = coef
                    Xk[indexB, state, 1] = occ
                    Xk[indexB, state, 2] = virt + nb
                    Xk[indexB, state, 3] = indexB
                    XAB[indexA] = "A"
                    XAB[indexB] = "B"
            elif "<-" in l:
                # Deexcitation coefficients
                parts = l.split()
                if len(parts) == 4:
                    occ, virt, coef = parts[0], parts[2], float(parts[3])
                elif len(parts) == 3:
                    occ, virt, coef = parts[0], parts[1][2:], float(parts[2])
                if "A" in occ:
                    # Alpha deexcitation
                    occ = int(occ.replace("A", "")) - 1
                    virt = int(virt.replace("A", "")) - na - 1
                    index = 1 + (norb-na)*occ + virt
                    Y[index, state] = coef
                    Yk[index, state, 0] = coef
                    Yk[index, state, 1] = occ
                    Yk[index, state, 2] = virt + na
                    Yk[index, state, 3] = index
                    YAB[index] = "A"
                elif "B" in occ:
                    # Beta deexcitation
                    occ = int(occ.replace("B", "")) - 1
                    virt = int(virt.replace("B", "")) - nb - 1
                    index = 1 + (norb-nb)*occ + virt + na*(norb-na)
                    Y[index, state] = coef
                    Yk[index, state, 0] = coef
                    Yk[index, state, 1] = occ
                    Yk[index, state, 2] = virt + nb
                    Yk[index, state, 3] = index
                    YAB[index] = "B"
                else:
                    # Deexcitation in closed-shell system
                    occ = int(occ) - 1
                    virt = int(virt) - na - 1
                    assert na == nb
                    indexA = 1 + (norb-na)*occ + virt
                    indexB = indexA + na*(norb-na)
                    Y[indexA, state] = coef
                    Y[indexB, state] = coef
                    Yk[indexA, state, 0] = coef
                    Yk[indexA, state, 1] = occ
                    Yk[indexA, state, 2] = virt + na
                    Yk[indexA, state, 3] = indexA
                    Yk[indexB, state, 0] = coef
                    Yk[indexB, state, 1] = occ
                    Yk[indexB, state, 2] = virt + nb
                    Yk[indexB, state, 3] = indexB
                    YAB[indexA] = "A"
                    YAB[indexB] = "B"
            elif "SavETr" in l:
                break

        # number of states that were actually read
        logging.info("read %d excited states" % state)
        nstates = state+1
        X = X[:, :nstates]
        Y = Y[:, :nstates]

        self.Xk   = Xk[:, :nstates, :]
        self.Yk   = Yk[:, :nstates, :]
        self.XpYk = Xk[:, :nstates, :] * 1.0
        self.XpYk[:, :, 0] += self.Yk[:, :, 0]
        self.XAB  = XAB
        self.YAB  = YAB

        self.nstates = nstates
        self.XpY = X+Y
        self.XmY = X-Y

        err = la.norm(np.dot(self.XpY.transpose(), self.XmY) - np.eye(nstates))
        # check normalization
        logging.info("|<X+Y|X-Y> - Id|= %e" % err)

        # compute total energies of states
        self.energies = np.zeros(nstates, dtype=float)
        self.energies = self.total_energy + self.excitation_energies[:nstates]

        return lines[i:]


    def readMolecularOrbitals(self, lines, last=False):
        """
        Parsing of AO coefficients of MOs
        """

        nl    = len(lines)
        nbfsC = self.basis.nbfsC
        nmo   = self.nmo
        orbe  = np.zeros(nmo, dtype=float)
        orbsC = np.zeros((nbfsC, nmo), dtype=float)
        for i in range(0, nl):
            l = lines[i]
            if "Molecular Orbital Coefficients" in l:
                if "Alpha" in l:
                    spin = "alpha"
                elif "Beta" in l:
                    spin = "beta"
                else:
                    spin = "both"
                break

        # read blocks of MO coefficients
        if nmo % 5 == 0:
            nblocks = nmo // 5
        else:
            nblocks = nmo // 5 + 1

        # count orbitals
        iorb = 0
        for n in range(0, nblocks):
            i += 3
            # read eigenvalues
            l = lines[i]
            parts = l.split()
            assert parts[0] == "Eigenvalues"
            norb = len(parts[2:])
            orbe[iorb:iorb+norb] = parts[-norb:]

            for k in range(0, nbfsC):
                i += 1
                l = lines[i]
                parts = l.split()
                orbsC[k, iorb:iorb+norb] = parts[-norb:]

                if n == 0:
                    # read atom centered basis
                    if len(parts) == norb + 4:
                        center = int(parts[1])-1
                    # powers (nx,ny,nz)
                    if len(parts) < norb+4:
                        powstr = parts[1]
                    else:
                        powstr = parts[3]

                    (nx, ny, nz) = (powstr.count('X'), powstr.count('Y'),
                                    powstr.count('Z'))
                    self.basis.powersC[:, k] = (nx, ny, nz)

                    # angular momentum shell (all orbitals within a shell share
                    # the same exponent and l-value
                    self.basis.shell_indexC[k] = int(powstr.replace(
                        "S", "").replace("X", "").replace("Y", "").replace(
                        "Z", "").replace("P", "").replace("D-2", "").replace(
                        "D-1", "").replace("D+1", "").replace(
                        "D+2", "").replace("D", "").replace("F-3", "").replace(
                        "F-2", "").replace("F-1", "").replace(
                        "F+1", "").replace("F+2", "").replace(
                        "F+3", "").replace("F", ""))

                    # centers of basis functions
                    self.basis.atom_idsC[k] = center
                    self.basis.centersC[:, k] = self.coordinates[:, center]

                    # exponents
                    exp_index = ""
                    for s in powstr:
                        if not s.isdigit():
                            break
                        else:
                            exp_index += s
                    exp_index = int(exp_index) - 1
                    self.basis.exponentsC[k] = self.exponents[center][exp_index]

            iorb += norb

        Omatrix = self.basis.decontractBasis()
        orbsD   = np.dot(Omatrix.T, orbsC)

        assert iorb == nmo
        if spin == "alpha":
            self.orbe_alpha = orbe
            self.orbs_alpha = orbsD
        elif spin == "beta":
            self.orbe_beta = orbe
            self.orbs_beta = orbsD
        else:
            self.orbe_alpha = orbe
            self.orbs_alpha = orbsD
            self.orbe_beta = orbe
            self.orbs_beta = orbsD

        if last:
            exptmp = []
            n = -1
            for j in self.exponents:
                exptmp.append([])
                n += 1
                for k in j:
                    for l in k:
                        exptmp[n].append(l[0])
            self.exponents = []
            self.exponents = exptmp

        return lines[i:]


    def computeCIScoefficients(self):
        Sao = dyA.basis_overlap(self.basis, self.basis)
        self.Sao = Sao.copy()

        self.cis = np.zeros(self.XpY.shape)
        # Ground state
        self.cis[0, 0] = 1.0
        self.cisk = np.ones((1, 1, 1))

        if not self.excited:
            self.ucis = self.cis / la.norm(self.cis, axis=0)[None, :]
            return

        # A-B
        AmB = dyA.AminusB_diagonal(self.orbe_alpha, self.orbe_beta,
                                   self.nelec_alpha, self.nelec_beta, self.norb)

        self.cisk = np.zeros(self.XpYk.shape)
        self.cisk[0, 0, 0] = 1.0

        # excited states
        for i in range(1, self.nstates):
            omega = self.excitation_energies[i]
            self.cis[:, i] = np.sqrt(omega) * self.XpY[:, i] / np.sqrt(AmB)
            self.cisk[:, i, 1:] = self.XpYk[:, i, 1:].copy()

        self.ucis = self.cis / la.norm(self.cis, axis=0)[None, :]
        self.cisk[:, :, 0] = self.ucis

        logging.info(np.dot(self.cis.transpose(), self.cis))


    def printResults(self):
        logging.info("%d alpha electrons  %d beta electrons" %
                        (self.nelec_alpha, self.nelec_beta))
        logging.info("%d basis functions" % self.basis.nbfsD)
        logging.info("Total SCF energy: %s Hartree\n" % self.total_energy)
        logging.info("number of electronic states: %d" % self.nstates)

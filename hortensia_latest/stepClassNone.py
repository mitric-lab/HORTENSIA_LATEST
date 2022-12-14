#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from configparser import ConfigParser

import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite as hp

from pyscf import gto

from joblib import Parallel, delayed

import hortensia_latest.misc as misc
import hortensia_latest.aid4c as aid4c
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

        self.exact2elInt = config['2e-Integrals']['exact2elInt'] == 'True'


    def ftMOpre(self):
        """
        calculates the overlap of anionic MOs and plane wave,
        time-independent part
        """

        exN        = self.anion.basis.exponentsD
        pow1       = self.anion.basis.powersD
        powers     = misc.powerToHerm(pow1)

        kx    = self.grid.kbasis[:, 0]
        ky    = self.grid.kbasis[:, 1]
        kz    = self.grid.kbasis[:, 2]

        c2    = 1/(2*exN)**1.5 * \
                dyA.norm(exN, pow1[0], pow1[1], pow1[2]) * \
                (1j / (2*np.sqrt(exN)))**np.sum(pow1, axis=0)
        eTerm = np.exp(-(kx[:, None]**2 + ky[:, None]**2 + kz[:, None]**2) /
                       (4 * exN[None, :]))

        H = np.zeros((len(kz), len(exN)))
        for i, exNi in enumerate(exN):
            H[:, i] = hp.hermval((kx / (2*np.sqrt(exNi))), powers[i, 0]) * \
                      hp.hermval((ky / (2*np.sqrt(exNi))), powers[i, 1]) * \
                      hp.hermval((kz / (2*np.sqrt(exNi))), powers[i, 2])

        SkMOpre  = c2[None, :] * eTerm * H

        return SkMOpre


    def MOoverlap(self, preMO):
        """
        calculates the overlap of anionic MOs and plane wave
        """

        kx    = self.grid.kbasis[:, 0]
        ky    = self.grid.kbasis[:, 1]
        kz    = self.grid.kbasis[:, 2]
        A     = self.anion.basis.centersD
        orbs  = self.anion.orbs_alpha

        c1    = np.exp(1j * (kx[:, None] * A[0][None, :] + \
                             ky[:, None] * A[1][None, :] + \
                             kz[:, None] * A[2][None, :]))

        SkMO  = c1 * preMO

        self.SkMOa  = aid4c.calcSkMOa(SkMO, orbs, self.grid.points,
                                      self.anion.basis.nbfsD, self.anion.nmo)


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
        en[self.nBound:] = self.ionized.total_energy - refEn + \
                            self.grid.kinEn
        self.Hdiag = en

        ### NOT DONE!!!
        self.diaTemp = diaTemp
        self.diacoup = dia

        return
        ###

        ### core-electron and electron-electron contribution to H-Matrix
        nel_a = self.anion.nelec_alpha
        nel_b = self.anion.nelec_beta
        nmo   = self.anion.nmo

        if len(self.states) == 1 and self.state == 0:
            # For if only the anionic ground state is considered for dynamics
            # The ground state only consists of a single Slater determinant,
            # therefore does not need to iterate over microstates
            Atemp = np.arange(nel_a)
            Btemp = np.arange(nel_b)

            Alist = np.zeros(nmo, dtype=bool)
            Blist = np.zeros(nmo, dtype=bool)

            for k in Atemp:
                Alist[k] = True
            for k in Btemp:
                Blist[k] = True

            if self.exact2elInt:
                eeInt = self.twoElInt(Alist, Blist, exInt) * \
                        self.grid.kfact
            else:
                eeIntTemp = self.twoElInt(Alist, Blist)
                eeInt = eeIntTemp * self.grid.kfact
                diaTemp[0, self.nBound:] += eeIntTemp

            dia[0, self.nBound:] += eeInt

        else:
            # For if excited anion states are involved
            # This makes iteration over microstates necessary
            precision = float(config['States']['precision']) / 100
            max_micro = int(config['States']['max_micro'])
            for i, state in enumerate(self.states):
                print("State: ", state)
                SDet   = self.anion.cisk[:, state, :]
                SDet   = SDet[(-abs(SDet[:, 0])).argsort()]

                cutoff = 0.0
                nDet   = 0
                for j, exc in enumerate(SDet):
                    print("  Microstate: ", j)
                    Atemp = np.arange(nel_a)
                    Btemp = np.arange(nel_b)

                    if self.anion.XAB[int(exc[3])] == "A":
                        Atemp[int(exc[1])] = exc[2]
                    else:
                        Btemp[int(exc[1])] = exc[2]

                    Alist = np.zeros(nmo, dtype=bool)
                    Blist = np.zeros(nmo, dtype=bool)

                    for k in Atemp:
                        Alist[k] = True
                    for k in Btemp:
                        Blist[k] = True

                    if self.exact2elInt:
                        eeInt = self.twoElInt(Alist, Blist, exInt) * \
                                self.grid.kfact
                    else:
                        eeIntTemp = self.twoElInt(Alist, Blist)
                        eeInt = eeIntTemp * self.grid.kfact
                        diaTemp[i, self.nBound:] += eeIntTemp * exc[0]

                    dia[i, self.nBound:] += eeInt * exc[0]

                    cutoff += exc[0]**2
                    nDet   += 1
                    if cutoff > precision or nDet == max_micro:
                        break

        self.diaTemp = diaTemp
        self.diacoup = dia


    def twoElInt(self, Alist, Blist, exInt=False):
        """
        Calculation of diabatic couplings

        For better understanding of the formulae in the following code, the
        diabatic couplings are separated into the terms 1, 2, 3 and 4,
        corresponding to
                                1              3
                            --------   ------------------
        V_i0 = N_k sum_lmn [(<kl|mn> - sum_r B_r <rl|mn>) (A_lmn + Ap_lmn) -
                            (<kl|nm> - sum_r B_r <rl|nm>) *A_lmn]
                            --------   ------------------
                                2              4

        Terms 1 and 2 as well as terms 3 and 4 are grouped together and
        calculated simultaneously. The numbers indicate the terms WITHOUT sign,
        therefore the correct total coupling term is finally produced with

        V_i0 = N_k * (term1minus2 - term3minus4)

        Also, although in the paper V_i0 is given, the program actually
        calculates V_0i, which is simply the complex conjugate of V_i0. This is
        due to the co-authors not being able to agreed upon which convention is
        more intuitive

        Annotation:
        orbs[i,j] : i = AOs; j = MOs
        """

        nel_a    = self.anion.nelec_alpha
        nel_b    = self.anion.nelec_beta
        nbfs     = self.anion.basis.nbfsD
        atom_ids = self.anion.basis.atom_idsD

        orbsAnA  = self.anion.orbs_alpha[:, Alist]
        orbsAnB  = self.anion.orbs_beta[:, Blist]
        orbsNeB  = self.ionized.orbs_beta[:, :nel_b]
        orbsNeA  = np.zeros((nbfs, nel_a))
        orbsNeA[:, :-1] = self.ionized.orbs_alpha[:, :(nel_a-1)]

        Sao      = self.anion.Sao
        SmoA     = la.multi_dot([orbsNeA.T, Sao, orbsAnA])
        SmoB     = la.multi_dot([orbsNeB.T, Sao, orbsAnB])
        detSmoB  = la.det(la.multi_dot([orbsNeB.T, Sao, orbsAnB]))

        detSA  = np.zeros((nel_a-1, nel_a, nel_a))
        for iMO in range(nel_a-1):
            for jMO in range(nel_a):
                for kMO in range(jMO+1, nel_a):
                    detSA[iMO, jMO, kMO] = self.SdetA(SmoA[:-1], iMO, jMO, kMO)
        detSA  *= detSmoB
        detSAp  = np.zeros((nel_b, nel_a, nel_b))
        for iMO in range(nel_b):
            for jMO in range(nel_a):
                for kMO in range(nel_b):
                    detSAp[iMO, jMO, kMO] = self.SdetAp(SmoA[:-1], SmoB,
                                                        iMO, jMO, kMO)

        A      = self.A(SmoA[:-1], detSA, Alist)
        Ap     = self.Ap(SmoB, detSAp, Alist, Blist, True)
        B      = np.einsum('jk,ij->ik', self.SkMOa[self.occMOs, :],
                                        self.anion.orbs_alpha[:, self.occMOs])

        term3minus4 = self.PYSCF(A, Ap, B)

        if self.exact2elInt:
            ### !! beta terms not implemented yet
            term1minus2 = self.exactPWterms(exInt, A)
            return self.renormM * (term1minus2 - term3minus4)
        else:
            term1minus2 = self.eConstApprox(nbfs, atom_ids, A, Ap)
            return self.renormM * (term1minus2 - term3minus4)


    def exactPWterms(self, exInt, A):
        Int    = exInt.twoElInt(self.coord)

        eTerm1 = np.einsum('ijk,mjik->m', A, Int)
        eTerm2 = np.einsum('ijk,mkij->m', A, Int)

        eTerm  = exInt.interpolate(self.grid, eTerm1-eTerm2)

        return eTerm


    def eConstApprox(self, nbfs, atom_ids, A, Ap):
        """
        Approximation where
        exp(ik*r) = exp(ik*R) * exp(-ik*R) * exp(ik*r)
                  = exp(ik*R) * exp(ik*(r-R))
                  ~ exp(ik*R) * (1 + ik*(r-R))
        leading to
        <k(1)l(2)|m(1)n(2)> = exp(ik*Rm)*(    <l(2)|m (1)n(2)>
                                          +ik*<l(2)|m'(1)n(2)>)
        So in fact, the approximation is NOT "eConst"
        """

        expkR = np.zeros((nbfs, self.grid.points), dtype=complex)
        for i in range(nbfs):
            expkR[i] = (2*np.pi)**(-1.5) * np.exp(1.0j *
                        np.dot(self.grid.kbasis, self.coord[atom_ids[i]]))

        # list of k components involved in first-order corrections to
        # constant PW approximation
        kcomp1 = 1j * self.grid.kbasis[:, 0]
        kcomp2 = 1j * self.grid.kbasis[:, 1]
        kcomp3 = 1j * self.grid.kbasis[:, 2]

        An = A + Ap

        # Zeroth order terms
        eTerms1a   = np.einsum('ijk,ikj->j', An   , self.int3c  )
        eTerms1C   = np.einsum('jm,j->m'   , expkR, eTerms1a    )
        eTerms2a   = np.einsum('ijk,ijk->k', A    , self.int3c  )
        eTerms2C   = np.einsum('km,k->m'   , expkR, eTerms2a    )

        eTerms11xa = np.einsum('ijk,ikj->j', An   , self.int3ckx)
        eTerms11xC = np.einsum('jm,m,j->m' , expkR, kcomp1      , eTerms11xa)
        eTerms21xa = np.einsum('ijk,ijk->k', A    , self.int3ckx)
        eTerms21xC = np.einsum('km,m,k->m' , expkR, kcomp1      , eTerms21xa)

        eTerms11ya = np.einsum('ijk,ikj->j', An   , self.int3cky)
        eTerms11yC = np.einsum('jm,m,j->m' , expkR, kcomp2      , eTerms11ya)
        eTerms21ya = np.einsum('ijk,ijk->k', A    , self.int3cky)
        eTerms21yC = np.einsum('km,m,k->m' , expkR, kcomp2      , eTerms21ya)

        eTerms11za = np.einsum('ijk,ikj->j', An   , self.int3ckz)
        eTerms11zC = np.einsum('jm,m,j->m' , expkR, kcomp3      , eTerms11za)
        eTerms21za = np.einsum('ijk,ijk->k', A    , self.int3ckz)
        eTerms21zC = np.einsum('km,m,k->m' , expkR, kcomp3      , eTerms21za)

        eTerms1m2  = (eTerms1C + eTerms11xC + eTerms11yC + eTerms11zC) - \
                     (eTerms2C + eTerms21xC + eTerms21yC + eTerms21zC)

        return eTerms1m2


    def A(self, S, detS, Alist):
        """
        Aijk = sum_ijk (-1)**P * det(tilde(S)) * (Ci - sum_l Sil Cl) * Cj * Ck
        """

        nel_a   = self.anion.nelec_alpha
        orbsAnA = self.anion.orbs_alpha[:, Alist]
        orbsNeA = self.ionized.orbs_alpha[:, :(nel_a-1)]

        P = np.ones((nel_a-1, nel_a, nel_a))
        P[::2, ::2, ::2]   = -1
        P[1::2, 1::2, ::2] = -1
        P[1::2, ::2, 1::2] = -1
        P[::2, 1::2, 1::2] = -1

        overtemp = np.einsum('il,al->ai', S, orbsAnA)
        Icoefs   = orbsNeA - overtemp

        temp1    = np.einsum('ijk,ijk,ai->ajk', P, detS, Icoefs)
        temp2    = np.einsum('ajk,bj->abk', temp1, orbsAnA)
        A        = np.einsum('abk,ck->abc', temp2, orbsAnA)

        return A


    def Ap(self, S, detS, Alist, Blist, beta=False):
        """
        Aijk = sum_ijk (-1)**P * det(tilde(S)) * (Ci - sum_l Sil Cl) * Cj * Ck
        """

        nel_a   = self.anion.nelec_alpha
        nel_b   = self.anion.nelec_beta
        orbsAnA = self.anion.orbs_alpha[:, Alist]
        orbsAnB = self.anion.orbs_beta[:, Blist]
        orbsNeB = self.ionized.orbs_beta[:, :(nel_b)]

        P = np.ones((nel_b, nel_a, nel_b))
        P[::2, ::2, ::2]   = -1
        P[1::2, 1::2, ::2] = -1
        P[1::2, ::2, 1::2] = -1
        P[::2, 1::2, 1::2] = -1

        overtemp = np.einsum('il,al->ai', S, orbsAnB)
        Icoefs   = orbsNeB - overtemp

        temp1    = np.einsum('ijk,ijk,ai->ajk', P, detS, Icoefs)
        temp2    = np.einsum('ajk,bj->abk', temp1, orbsAnA)
        A        = np.einsum('abk,ck->abc', temp2, orbsAnB)

        return A


    def SdetA(self, S, i, j, k):
        """
        calculation of det(Salpha) with a modified overlap matrix S where
        the ith/last line and jth/kth column are omitted
        """

        nel_a      = self.anion.nelec_alpha
        temp       = np.ones((nel_a-1, nel_a), dtype=bool)
        temp[i, :] = 0
        temp[:, j] = 0
        temp[:, k] = 0

        Sdet = S[temp]

        return la.det(Sdet.reshape(nel_a-2, nel_a-2))


    def SdetAp(self, Sa, Sb, i, j, k):
        """
        calculation of det(Stotal) with a modified overlap matrix S where
        the ith beta and last alpha line and jth alpha and kth beta column
        are omitted
        """

        nel_a       = self.anion.nelec_alpha
        nel_b       = self.anion.nelec_beta
        tempA       = np.ones((nel_a-1, nel_a), dtype=bool)
        tempA[:, j] = 0
        tempB       = np.ones((nel_b, nel_b), dtype=bool)
        tempB[i, :] = 0
        tempB[:, k] = 0

        SdetA      = Sa[tempA]
        SdetB      = Sb[tempB]

        Sdet       = la.det(SdetA.reshape(nel_a-1, nel_a-1)) * \
                     la.det(SdetB.reshape(nel_b-1, nel_b-1))

        return Sdet


    def PYSCF(self, A, Ap, B, Int2=False):
        """
        calculates two electron integrals using the PYSCF program
        evaluating either two, three or four center integrals
        """

        def createCoord(coord, atnr):
            coordS  = ""
            for i, R in enumerate(coord*misc.bohr_to_angs):
                coordS += "%s %.9f %.9f %.9f;"%(misc.NrToS(atnr[i]),
                                                 R[0], R[1], R[2])
                atseq.append(misc.NrToS(atnr[i]))
            coordS  = coordS[:-1]

            return coordS

        def createMol(coord, powVal, atUsed, auxVal=0):
            """
            creation of molecule object in PYSCF
            """

            spin1  = (np.sum(self.atnr[atUsed])-self.charge)%2

            mol    = gto.M(atom=coord, charge=self.charge, spin=spin1)

            basis  = []
            atsym  = []
            atnbfs = []

            n = 0
            for i, j in enumerate(self.anion.exponents):
                if not atUsed[i]:
                    n += nPowPerC[i]
                    continue

                if misc.NrToS(self.atnr[i]) in atsym:
                    n += atnbfs[atsym.index(misc.NrToS(self.atnr[i]))]
                else:
                    basistemp = ""
                    nbfs      = 0
                    for k in j:
                        if powers[n] in powVal:
                            basistemp += "%s %s\n"%(misc.NrToS(self.atnr[i]),
                                                    ldict[powers[n]+auxVal])
                            basistemp += "%.9f 1.0\n"%k
                        n    += 2*powers[n] + 1 + powers[n] * (powers[n]-1) // 2
                        nbfs += 2*powers[n] + 1 + powers[n] * (powers[n]-1) // 2

                    basis.append(basistemp)
                    atsym.append(misc.NrToS(self.atnr[i]))
                    atnbfs.append(nbfs)

            mol.basis = {}
            for i, j in enumerate(atsym):
                mol.basis[j] = gto.basis.parse(basis[i])

            mol.build()

            return mol

        def renorm(mol, atNr, powVal):
            exp = []
            for i in atNr:
                a = mol.basis[atseq[i]]
                for j in a:
                    exp.append(j[1][0])

            exp  = np.asarray(exp)
            ni   = 2*powVal + 1 + powVal * (powVal-1) // 2
            nj   = 2*powVal + 3 + powVal * (powVal+1) // 2
            nrms = np.zeros((len(exp), ni, nj))

            li, lj = lorder[powVal], lorder[powVal+1]
            for a, expa in enumerate(exp):
                for i in range(ni):
                    for j in range(nj):
                        nrms[a, i, j] = (misc.normGauss(expa, li[i][0],
                                         li[i][1], li[i][2]) /
                                         misc.normGauss(expa, lj[j][0],
                                         lj[j][1], lj[j][2]))

            return nrms


        ldict   = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G"}
        inddict = {0: [0], 1: [1, 2, 3], 2: [4, 5, 6, 7, 8, 9]}
        # order in PYSCF
        lorder  = {0: [[0, 0, 0]],
                   1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                   2: [[2, 0, 0], [1, 1, 0], [1, 0, 1],
                       [0, 2, 0], [0, 1, 1], [0, 0, 2]],
                   3: [[3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0], [1, 1, 1],
                       [1, 0, 2], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]}

        atseq   = []
        coord   = ""
        for i, R in enumerate(self.coord*misc.bohr_to_angs):
            coord += "%s %.9f %.9f %.9f;"%(misc.NrToS(self.atnr[i]),
                                            R[0], R[1], R[2])
            atseq.append(misc.NrToS(self.atnr[i]))
        coord  = coord[:-1]

        powers = sum(self.anion.basis.powersD)

        atNr   = np.arange(self.nat)
        nPowPerC = np.zeros(len(self.coord), dtype=int)
        for i in self.anion.basis.atom_idsD:
            nPowPerC[i] += 1
        maxPow = np.zeros(np.shape(nPowPerC), dtype=int)

        n = 0
        for i, powi in enumerate(nPowPerC):
            ptemp = max(powers[n:n+powi])
            maxPow[i] = ptemp
            n += powi

        # generation of molecule in PYSCF
        atUsed   = np.ones(len(self.coord), dtype=bool)
        mol      = createMol(coord, range(np.max(powers)+1), atUsed)
        self.mol = mol

        # auxiliary molecule with p functions bearing the s exponents
        smol   = createMol(coord, [0], atUsed, 1)
        molps  = mol + smol
        # renormalization of p-corrections for s functions ( *s-norm/p-norm)
        nrms   = renorm(smol, atNr, 0)

        # auxiliary molecule with d functions bearing the p exponents
        atUsed = maxPow >= 1
        pAvail = atUsed.any()
        if pAvail:
            pcoord = createCoord(self.coord[atUsed], self.atnr[atUsed])
            pmol   = createMol(pcoord, [1], atUsed, 1)
            moldp  = mol + pmol
            # renormalization of d-corrections for p functions ( *p-norm/d-norm)
            nrmp   = renorm(pmol, atNr[atUsed], 1)

        # auxiliary molecule with f functions bearing the d exponents
        atUsed = maxPow >= 2
        dAvail = atUsed.any()
        if dAvail:
            dcoord = createCoord(self.coord[atUsed], self.atnr[atUsed])
            dmol   = createMol(dcoord, [2], atUsed, 1)
            molfd  = mol + dmol
            # renormalization of d-corrections for p functions ( *d-norm/f-norm)
            nrmd   = renorm(dmol, atNr[atUsed], 2)

        # indices indicating the concrete form of gto:
        # 0:s, 1:px, 2:py, 3:pz, 4:dx2, 5:dxy, 6:dxz, 7:dy2, 8:dyz, 9: dz2
        self.indexlist  = []
        # indices indicating the sequence of l quantum numbers in the basis set:
        # 0:s, 1:p, 2:d
        self.indexlist2 = []

        n = 0
        for i, j in enumerate(self.anion.exponents):
            for k in j:
                for l in inddict[powers[n]]:
                    self.indexlist.append(l)

                self.indexlist2.append(powers[n])
                n += 2*powers[n]+1 + powers[n]*(powers[n]-1)//2

        # PYSCF to G09 conversion
        # PYSCF: x2, xy, xz, y2, yz, z2
        # G09  : x2, y2, z2, xy, xz, yz
        # PYSCF to G09 index change: [ 0, 2, 2,-2, 1,-3]
        # G09 to PYSCF index change: [ 0, 2, 3,-2,-2,-1]
        # since we are calculating in G09 convention, the former applies
        pyToG09 = np.array([0, 2, 2, -2, 1, -3])
        i = 0
        sortAO = np.arange(len(powers))
        while i < len(powers):
            if powers[i] == 2:
                sortAO[i:i+6] += pyToG09
                i += 5
            i += 1

        g09ToPy = np.array([0, 2, 3, -2, -2, -1])
        i = 0
        sortAOb = np.arange(len(powers))
        while i < len(powers):
            if powers[i] == 2:
                sortAOb[i:i+6] += g09ToPy
                i += 5
            i += 1

        ###############################

        # renormalization from overlap due to wrong normalization of d functions
        overlap  = self.mol.intor('int1e_ovlp_cart')
        nrm1     = 1/np.sqrt(np.diagonal(overlap))

        if pAvail:
            overlap2 = moldp.intor('int1e_ovlp_cart', shls_slice = (
                                    mol.nbas, mol.nbas+pmol.nbas,
                                    mol.nbas, mol.nbas+pmol.nbas))
            nrm2     = 1/np.sqrt(np.diagonal(overlap2))

        if dAvail:
            overlap3 = molfd.intor('int1e_ovlp_cart', shls_slice = (
                                    mol.nbas, mol.nbas+dmol.nbas,
                                    mol.nbas, mol.nbas+dmol.nbas))
            nrm3     = 1/np.sqrt(np.diagonal(overlap3))

        ###############################

        ### two-center integrals (currently not used)
        if Int2:
            self.int2c = mol.intor('int2c2e_cart')
            self.int2c = self.int2c[sortAO[:, None], sortAO[None, :]]

        ### three-center integrals used for constant-wave approximation
        self.int3c  = mol.intor('int3c2e_cart')
        # renormalization of integrals due to non-normalized d functions
        self.int3c *= nrm1[:, None, None] * nrm1[None, :, None] * \
                      nrm1[None, None, :]
        # reordering of terms due to different sorting of d-orbitals in PYSCF
        self.int3c  = self.int3c[sortAO[:, None, None],
                                 sortAO[None, :, None],
                                 sortAO[None, None, :]]

        # correction to the three-center integrals for s functions
        int3cs  = molps.intor('int3c2e_cart', shls_slice = (
                                    0       , mol.nbas,
                                    0       , mol.nbas,
                                    mol.nbas, mol.nbas+smol.nbas))
        int3cs *= nrm1[:, None, None] * nrm1[None, :, None]
        int3cs  = int3cs[sortAO[:, None], sortAO[None, :], :].copy()

        # correction to the three-center integrals for p functions
        if pAvail:
            int3cp  = moldp.intor('int3c2e_cart', shls_slice = (
                                        0       , mol.nbas,
                                        0       , mol.nbas,
                                        mol.nbas, mol.nbas+pmol.nbas))
            int3cp *= nrm1[:, None, None] * nrm1[None, :, None] * \
                      nrm2[None, None, :]
            int3cp  = int3cp[sortAO[:, None], sortAO[None, :], :].copy()

        # correction to the three-center integrals for d functions
        if dAvail:
            int3cd  = molfd.intor('int3c2e_cart', shls_slice = (
                                        0       , mol.nbas,
                                        0       , mol.nbas,
                                        mol.nbas, mol.nbas+dmol.nbas))
            int3cd *= nrm1[:, None, None] * nrm1[None, :, None] * \
                      nrm3[None, None, :]
            int3cd  = int3cd[sortAO[:, None], sortAO[None, :], :].copy()

        ### four-center integrals for terms arising from the
        ### orthogonalization of the plane wave

        # Transform into PYSCF sorting
        A  = A[sortAOb[:, None, None],
               sortAOb[None, :, None],
               sortAOb[None, None, :]]
        Ap = Ap[sortAOb[:, None, None],
                sortAOb[None, :, None],
                sortAOb[None, None, :]]
        B  = B[sortAOb]

        lnCdict = {'s': 1, 'p': 3, 'd': 6}
        lnSdict = {'s': 1, 'p': 3, 'd': 5}
        b = []
        n = 0
        for i in mol.ao_labels():
            if n > 0:
                n -= 1
                continue
            for j in i.split()[2]:
                if j in ['s', 'p', 'd']:
                    b.append(lnCdict[j])
                    n = lnSdict[j]-1
                    break

        # division of basis function iteration into Nproc parallel processes
        Nproc = int(config['QC']['Nproc'])
        cpu1, cpu2 = [], []
        for i in range(Nproc):
            cpu1.append([round(mol.nbas / Nproc * i),
                         round(mol.nbas / Nproc * (i+1))])
            cpu2.append(sum(b[:cpu1[-1][0]]))
        cpu1 = np.asarray(cpu1, dtype=np.int32)

        q = Parallel(n_jobs=Nproc)(
            delayed(aid4c.int2eParallel)(
            A, Ap, B, nrm1, self.anion.basis.nbfsD, mol, cpu1[i], cpu2[i])
            for i in range(Nproc))

        term3m4 = np.zeros(len(B[0]), dtype=np.complex)
        for i in q:
            term3m4 += i
        q = None

        ##############################
        # sorting of int3cs and int3cp correction terms into arrays of the same
        # shape as self.int3c, i.e., the corrections are put at the same index
        # values as the terms in self.int3c
        shape3c      = np.shape(self.int3c)
        self.int3ckx = np.zeros(shape3c)
        self.int3cky = np.zeros(shape3c)
        self.int3ckz = np.zeros(shape3c)

        # sequence of d functions in pyscf: x2, xy, xz, y2, yz, z2
        # f functions: x3, x2y, x2z, xy2, xyz, xz2, y3, y2z, yz2, z3
        conv = np.array([[0, 1, 2], [1, 3, 4], [2, 4, 5],
                         [3, 6, 7], [4, 7, 8], [5, 8, 9]])

        scnt, pcnt, dcnt = 0, 0, 0
        nbas = len(nrm1)
        for i in range(nbas):
            if self.indexlist[i] == 0:   #s
                self.int3ckx[:, :, i] = int3cs[:, :, 3*scnt] \
                                        * nrms[scnt, 0, 0]
                self.int3cky[:, :, i] = int3cs[:, :, 3*scnt+1] \
                                        * nrms[scnt, 0, 0]
                self.int3ckz[:, :, i] = int3cs[:, :, 3*scnt+2] \
                                        * nrms[scnt, 0, 0]
                scnt += 1

            elif self.indexlist[i] == 1: #p
                for j in range(3):
                    self.int3ckx[:, :, i+j] = int3cp[:, :, 6*pcnt+conv[j, 0]] \
                                              * nrmp[pcnt, j, conv[j, 0]]
                    self.int3cky[:, :, i+j] = int3cp[:, :, 6*pcnt+conv[j, 1]] \
                                              * nrmp[pcnt, j, conv[j, 1]]
                    self.int3ckz[:, :, i+j] = int3cp[:, :, 6*pcnt+conv[j, 2]] \
                                              * nrmp[pcnt, j, conv[j, 2]]
                pcnt += 1

            elif self.indexlist[i] == 4: #d
                for j in range(6):
                    self.int3ckx[:, :, i+j] = int3cd[:, :, 10*dcnt+conv[j, 0]] \
                                              * nrmd[dcnt, j, conv[j, 0]]
                    self.int3cky[:, :, i+j] = int3cd[:, :, 10*dcnt+conv[j, 1]] \
                                              * nrmd[dcnt, j, conv[j, 1]]
                    self.int3ckz[:, :, i+j] = int3cd[:, :, 10*dcnt+conv[j, 2]] \
                                              * nrmd[dcnt, j, conv[j, 2]]
                dcnt += 1

        return term3m4

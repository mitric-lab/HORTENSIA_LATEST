#!/usr/bin/env python3

import logging
from configparser import ConfigParser

import numpy as np
import numpy.linalg as la
from scipy.integrate import ode
from numpy.random import rand

import hortensia_latest.misc as misc
import hortensia_latest.dysonAux as dyA
import hortensia_latest.exactIntegral as exI

config = ConfigParser()
config.read('config.ini')

orthotype = config['States']['orthotype'].lower()

if orthotype == 'mo':
    import hortensia_latest.stepClassMO as stepC
elif orthotype == 'state':
    import hortensia_latest.stepClassState as stepC
elif orthotype == 'dyson':
    import hortensia_latest.stepClassDyson as stepC
elif orthotype == 'none':
    import hortensia_latest.stepClassNone as stepC

################################################################################
#                                                                              #
#                    Population Dynamics Calculation Class                     #
#                                                                              #
################################################################################


class PopulationDynamics:
    def __init__(self, states, state, Eshift, grid, coord, atS, fo, index=""):
        """
        initializes population dynamics
        """

        self.states  = states
        if config['Dynamics']['restart'] == 'True':
            self.statenr = int(np.loadtxt("boundPop.dat")[-1, -1])
            self.state   = self.states.index(self.statenr)
        else:
            self.state = state
            self.statenr = states[state]

        self.orthotype = orthotype

        self.Eshift  = Eshift
        self.grid    = grid
        self.nBound  = len(self.states)
        self.nstates = self.nBound + self.grid.points
        self.dt      = float(config['Dynamics']['timestep']) / misc.autime2fs
        self.dtfs    = float(config['Dynamics']['timestep'])
        # Creates an instance of the Step class to parse quantum-chemistry
        # output files
        self.Step1   = stepC.Step(self.states, self.statenr, self.Eshift,
                                  self.grid, "1", "1")
        self.Step1.coord = coord
        # Precalculates time-independent parts of the terms in the underlying
        # dynamics equations
        self.preFT   = self.Step1.ftDysonPre(fo)
        if self.orthotype in ["mo", "none"]:
            self.preMO   = self.Step1.ftMOpre()
            # Calculates the overlap of all plane waves and all anion MOs at
            # time step 1
            self.Step1.MOoverlap(self.preMO)

        if config['2e-Integrals']['exact2elInt'] == 'True':
            self.exInt = exI.ExactInt(atS, self.grid)

        self.index   = str(index)
        self.prob    = np.zeros(self.nstates+1, dtype=np.complex)
        # oldsign accounts for sign changes in wavefunctions of different steps
        # oldsign[0]    -> signs for anionic wavefunctions
        # oldsign[1, 0] -> sign for ionized groundstate,
        #                  other arguments on axis 1 are not used
        self.oldsign = np.ones((2, self.nBound))

        if config['Dynamics']['restart'] == 'True':
            self.trajpop = [int(np.loadtxt("trajpop%s.dat"%self.index)[-1, 1])]
            self.trajpop.append(
                int(np.loadtxt("trajpopWP%s.dat"%self.index)[-1, 1]))
            tempcoef     = np.loadtxt("tempcoef%s.dat"%self.index)
            self.c       = tempcoef[0, :] + 1j * tempcoef[1, :]
            if self.orthotype != "none":
                self.Step1.orthoRenorm(fo)
        else:
            self.c             = np.zeros(self.nstates, dtype=np.complex)
            self.c[self.state] = 1.0
            self.trajpop       = [int(config['Dynamics']['nrpertraj']),
                                  int(config['Dynamics']['nrpertraj'])]

            with open("trajpop%s.dat"%self.index, "w") as out:
                out.write("%.2f %i\n"%(0.0, self.trajpop[0]))

            with open("trajpopWP%s.dat"%self.index, "w") as out:
                out.write("%.2f %i\n"%(0.0, self.trajpop[1]))

            with open("boundPop%s.dat"%self.index, "w") as out:
                out.write("%.2f %i\n"%(0.0, self.statenr))

            if int(config['Dynamics']['printlevel']) >= 0:
                with open("norm%s.dat"%self.index, "w") as out:
                    out.write("%7.2f %5.9f\n"%(0.0, 1.0))

            if config['Dynamics']['coupprint'] == 'full':
                with open("couplings%s.dat"%self.index, "w") as out:
                    out.write("%7.2f "%0.0 + len(self.c) *
                              " %9.5e "%tuple(0.0 for _ in self.c) + "\n")
                for j in range(self.nBound):
                    if self.nBound == 1 and self.statenr == 0:
                        break
                    coupName = "couplings%s_to_%i.dat"%(self.index,
                                                        self.states[j])
                    with open(coupName, "w") as out:
                        out.write("%7.2f "%0.0 + len(self.c) *
                                  " %9.5e "%tuple(0.0 for _ in self.c) + "\n")
            elif config['Dynamics']['coupprint'] == 'summed':
                with open("couplings%s.dat"%self.index, "w") as out:
                    out.write("%7.2f "%0.0)
                    for i in range(self.nBound):
                        out.write("%9.5e "%0.0)
                    out.write("%9.5e\n"%0.0)
                for j in range(self.nBound):
                    if self.nBound == 1 and self.statenr == 0:
                        break
                    coupName = "couplings%s_to_%i.dat"%(self.index,
                                                        self.states[j])
                    with open(coupName, "w") as out:
                        out.write("%7.2f "%0.0)
                        for i in range(self.nBound):
                            out.write("%9.5e "%0.0)
                        out.write("%9.5e\n"%0.0)

            if config['Dynamics']['coefprint'] == 'full':
                with open("coefficients%s.dat"%self.index, "w") as out:
                    out.write("%7.2f " % 0.0 + len(self.c) *
                              " %9.5e " % tuple(abs(i)**2 for i in self.c) +
                              "\n")
            elif config['Dynamics']['coefprint'] == 'summed':
                with open("coefficients%s.dat"%self.index, "w") as out:
                    out.write("%7.2f "%0.0)
                    for i in range(self.nBound):
                        out.write("%9.5e "%abs(self.c[i]))
                    out.write("%9.5e\n"%0.0)

            if config['Dynamics']['probprint'] == 'full':
                with open("probabilities%s.dat"%self.index, "w") as out:
                    out.write("%7.2f " % 0.0 + len(self.c) *
                              " %9.5e " % tuple(abs(i)**2 for i in self.c) +
                              "\n")
            elif config['Dynamics']['probprint'] == 'summed':
                with open("probabilities%s.dat"%self.index, "w") as out:
                    out.write("%7.2f "%0.0)
                    for i in range(self.nBound):
                        out.write("%9.5e "%abs(self.c[i]))
                    out.write("%9.5e %9.5e\n"%(0.0, 0.0))

            fo.write("Calculating Dyson orbitals\n")
            self.refEn = self.Step1.herm_wfn(fo, self.preFT, firstStep=True)

            if config['2e-Integrals']['exact2elInt'] == 'True':
                self.Step1.calcHmatrix(self.refEn, self.exInt)
            else:
                self.Step1.calcHmatrix(self.refEn)

            if config['Dynamics']['enerprint'] == 'full':
                with open("energies%s.dat"%self.index, "w") as out:
                    out.write("%7.2f %12.9e "%(0.0,
                              self.Step1.Hdiag[self.state]))
                    out.write(len(self.Step1.Hdiag) *
                              "%12.9e "%tuple(i for i in self.Step1.Hdiag) +
                              "\n")
            elif config['Dynamics']['enerprint'] == 'none':
                with open("energies%s.dat"%self.index, "w") as out:
                    neuEn = self.Step1.ionized.total_energy - self.refEn
                    out.write("%7.2f %12.9e "%(0.0,
                              self.Step1.Hdiag[self.state]))
                    for i in range(self.nBound):
                        out.write("%12.9e "%self.Step1.Hdiag[i])
                    out.write("%12.9e\n"%neuEn)


    def integrate(self, step, coordold, coord, kinEn, fo):
        """
        Integration of electronic Schrödinger equation and evaluation of
        hopping instances through hopping probabilities
        """

        def func(t, a, arg):
            """
            Arguments in arg:
            arg[0] = transformation matrix to interaction picture without time
            arg[1] = transpose of arg[0]
            arg[2] = first rows of the total coupling matrix
            arg[3] = first columns of the total coupling matrix
            arg[4] = electronic resonance decay constant
            """
            T1 = arg[0]**t
            T2 = arg[1]**t
            apBound = np.einsum("ij,ij,j->i", T1, arg[2], a)
            apFree  = np.einsum("ji,ij,i->j", T2, arg[3], a[:self.nBound])
            apFree[:self.nBound] += apBound - arg[4] * a[:self.nBound]

            return apFree

        # Creates an instance of the Step class to parse quantum-chemistry
        # output files
        self.Step2 = stepC.Step(self.states, self.statenr, self.Eshift,
                                self.grid, "2", "2")
        self.Step2.coord = coord
        if self.orthotype in ["mo", "none"]:
            # Calculates the overlap of all plane waves and all anion MOs at
            # time step 2
            self.Step2.MOoverlap(self.preMO)
        # Calculates Dyson orbitals and their Fourier transforms
        self.Step2.herm_wfn(fo, self.preFT, refEn=self.refEn)

        # Calculates the Hamiltonian matrix (diagonal + diabatic couplings)
        if config['2e-Integrals']['exact2elInt'] == 'True':
            # Not tested, therefore not working
            self.Step2.calcHmatrix(self.refEn, self.exInt)
        else:
            self.Step2.calcHmatrix(self.refEn)

        dt = self.dt / int(config['Dynamics']['popsteps'])
        # Calculates non-adiabatic couplings
        # Hod are the off-diagonal elements of the Hamiltonian matrix,
        # Hd the diagonal terms and D the non-adiabatic coupling terms
        Hod, Hd, D = self.couplings(coordold, coord, fo)
        arg = []

        # Constructs necessary terms for integration of el. Schrödinger equation
        # using the interaction picture
        Tinter1 = np.zeros((self.nBound, self.nstates), dtype=complex)
        Tinter2 = np.zeros((self.nstates, self.nBound), dtype=complex)
        for i in range(self.nBound):
            Tinter1[i, self.nBound:] = np.exp(-1.0j * \
                                              (Hd[self.nBound:] - Hd[i]))
            Tinter2[self.nBound:, i] = np.exp(-1.0j * \
                                              (Hd[i] - Hd[self.nBound:]))
            for j in range(i+1, self.nBound):
                Tinter1[i, j] = np.exp(-1.0j * (Hd[j] - Hd[i]))
                Tinter1[j, i] = np.exp(-1.0j * (Hd[i] - Hd[j]))
                Tinter2[i, j] = np.exp(-1.0j * (Hd[j] - Hd[i]))
                Tinter2[j, i] = np.exp(-1.0j * (Hd[i] - Hd[j]))

        arg.append(Tinter1)
        arg.append(Tinter2)
        arg.append(-1.0j*Hod - 1/(2*self.dt) * D)
        arg.append(-1.0j*np.conjugate(Hod) + 1/(2*self.dt) * np.conjugate(D))

        neuEn = self.Step2.ionized.total_energy
        # If IE is negative, turns on electronic resonance decay
        elDecay = np.zeros(self.nBound)
        for i, en in Hd[:self.nBound]:
            if self.refEn + en <= neuEn:
                elDecay[i] = misc.autime2fs / \
                             (2 * float(config['States']['tau']))
        arg.append(elDecay)

        r = ode(func).set_integrator("zvode", method="adams",
                                     with_jacobian=False)
        r.set_initial_value(self.oldc, 0.0)
        r.set_f_params(arg)

        # Actual integration
        while r.successful() and r.t < self.dt:
            r.integrate(r.t + dt)

        # Retransformation from interaction picture
        self.c    = r.y * np.exp(-1.0j * Hd * r.t)
        # Calls the method evaluating the hopping probabilities
        self.prob = self.calcProb(self.c, self.oldc, fo)

        # Evaluation of autoionization hops
        fo.write("\nCURRENT ELECTRONIC STATE: %i\n"%(self.states[self.state]))
        self.trajpop[0] = self.singleWaveHop(rand(self.trajpop[0]),
                                             self.trajpop[0], step, kinEn, fo)
        self.trajpop[1] = self.wavepacketHop(rand(self.trajpop[1]),
                                             self.trajpop[1], step, kinEn, fo)
        # Evaluation of bound state (surface hopping) hop
        if self.nBound > 1:
            kinEn = self.boundHop(rand(1), kinEn, fo)

        # Writes all electronic properties to output files
        self.writeDat(step)
        # Shifts the Step2 instance on Step1
        self.shiftSteps() #!!!

        return [self.state, self.statenr, kinEn]


    def couplings(self, coord1, coord2, fo):
        """
        Calculation of non-adiabatic couplings
        """

        self.oldc   = self.c.copy()

        self.Cross1 = stepC.Step(self.states, self.statenr, self.Eshift,
                                 self.grid, "1", "2", self.Step2.EA)
        self.Cross2 = stepC.Step(self.states, self.statenr, self.Eshift,
                                 self.grid, "2", "1", self.Step1.EA)
        self.Cross1.coord = coord1
        self.Cross2.coord = coord2
        if self.orthotype in ["mo", "none"]:
            self.Cross1.MOoverlap(self.preMO)
            self.Cross2.MOoverlap(self.preMO)
        fo.write("\ncalculating cross term <phi_N(t)|phi_N-1(t+dt)psi^k>\n")
        self.Cross1.herm_wfn(fo, self.preFT, refEn=self.refEn)
        fo.write("\ncalculating cross term <phi_N(t+dt)|phi_N-1(t)psi^k>\n")
        self.Cross2.herm_wfn(fo, self.preFT, refEn=self.refEn)

        # D_ij = 1/(2*dt) * [<i(t)|j(t+dt)>-<i(t+dt)|j(t)>] = 1/(2*dt) * self.D
        self.D       = np.zeros((self.nBound, self.nstates), dtype=np.complex)
        self.Dtemp   = np.zeros((self.nBound, self.nstates), dtype=np.complex)
        self.diaTemp = np.zeros((self.nBound, self.nstates), dtype=np.complex)
        self.diaCoup = np.zeros((self.nBound, self.nstates), dtype=np.complex)

        # Non-adiabatic coupling between bound states
        Dbound, newsign = self.nonAdBound(self.oldsign)
        self.D[:self.nBound, :self.nBound]     = Dbound
        self.Dtemp[:self.nBound, :self.nBound] = Dbound

        if self.orthotype == "mo":
            self.nonAdOrthoMO(newsign)
        elif self.orthotype == "state":
            self.nonAdOrthoState(newsign)
        elif self.orthotype == "dyson":
            self.nonAdOrthoDyson(newsign)
        elif self.orthotype == "none":
            ### NOT DONE!!!
            self.nonAdOrthoNone(newsign)

        self.oldsign = newsign.copy()

        return self.diaCoup, self.Step2.Hdiag, self.D


    def nonAdOrthoMO(self, newsign):
        """
        Non-adiabatic couplings between bound anion and free electron states
        Plane waves orthogonalized with respect to the anion MOs
        """
        oldsign = self.oldsign.copy()

        self.diaTemp = self.Step1.diaTemp * oldsign[0][:, None]
        self.diaCoup = self.Step1.diacoup * oldsign[0][:, None]

        # Anion AO basis overlap <AO(t+dt)|AO(t)>
        Sao1 = dyA.basis_overlap(self.Step2.anion.basis,
                                 self.Cross1.anion.basis)
        # Anion AO basis overlap <AO(t)|AO(t+dt)>
        Sao2 = dyA.basis_overlap(self.Step1.anion.basis,
                                 self.Cross2.anion.basis)
        for i in range(self.nBound):
            orbAnA1 = self.Step1.anion.orbs_alpha[:, self.Step1.occMOs]
            orbAnA2 = self.Step2.anion.orbs_alpha[:, self.Step2.occMOs]
            # Overlap of Dyson orbital and anion MOs
            SmoDy1  = la.multi_dot([orbAnA2.T, Sao1,
                                    self.Cross1.dorbs[:, 0, i]])
            SmoDy2  = la.multi_dot([orbAnA1.T, Sao2,
                                    self.Cross2.dorbs[:, 0, i]])

            # non-adiabatic couplings between bound and free-electron state
            self.D[i, self.nBound:] = self.grid.kfact * (
                self.Step2.renormM * oldsign[0, i] * newsign[1, 0] * (
                self.Cross1.Smatrix[i, self.nBound:] -
                np.einsum('ik,i->k', self.Step2.SkMOa[self.Step2.occMOs, :],
                                     SmoDy1))
                -
                self.Step1.renormM * oldsign[1, 0] * newsign[0, i] * (
                self.Cross2.Smatrix[i, self.nBound:] -
                np.einsum('ik,i->k', self.Step1.SkMOa[self.Step1.occMOs, :],
                                     SmoDy2)))
            self.Dtemp[i, self.nBound:] = self.D[i, self.nBound:] / \
                                          self.grid.kfact


    def nonAdOrthoState(self, newsign):
        """
        Non-adiabatic couplings between bound anion and free electron states
        Free electron states orthogonalized with respect to anion states
        """

        oldsign = self.oldsign.copy()
        overt = self.stateOverlapTime()

        for i in range(self.nBound):
            j = 1 * self.nBound

            self.D[i, j:] = self.grid.kfact * (
                self.Step2.renormS * oldsign[0, i] * newsign[1, 0] * (
                self.Cross1.Smatrix[i, j:] -
                np.einsum("j,jk->k", overt[i], self.Step2.Smatrix[:, j:]))
                -
                self.Step1.renormS * oldsign[1, 0] * newsign[0, i] * (
                self.Cross2.Smatrix[i, j:] -
                np.einsum("j,jk->k", overt[:, i], self.Step1.Smatrix[:, j:])))

            self.Dtemp[i, j:] = self.D[i, j:] / self.grid.kfact


    def nonAdOrthoDyson(self, newsign):
        """
        Non-adiabatic couplings between bound anion and free electron states
        Plane waves orthogonalized with respect to the Dyson orbitals
        """

        oldsign = self.oldsign.copy()

        Sdy = dyA.dyson_overlap_K(self.Step1, self.Step2,
                                  self.nBound, self.states)

        for i in range(self.nBound):
            j = 1 * self.nBound

            self.D[i, j:] = self.grid.kfact * (
                self.Step2.renormD * oldsign[0, i] * newsign[1, 0] * (
                self.Cross1.Smatrix[i, j:] -
                np.einsum("n,nk->k", Sdy[i], self.Step2.Smatrix[:, j:]))
                -
                self.Step1.renormD * oldsign[1, 0] * newsign[0, i] * (
                self.Cross2.Smatrix[i, j:] -
                np.einsum("n,nk->k", Sdy[:, i], self.Step1.Smatrix[:, j:])))

            self.Dtemp[i, j:] = self.D[i, j:] / self.grid.kfact


    def nonAdOrthoNone(self, newsign):
        """
        Non-adiabatic couplings between bound anion and free electron states
        Plane waves orthogonalized with respect to the anion MOs
        """
        oldsign = self.oldsign.copy()

        self.diaTemp = self.Step1.diaTemp * oldsign[0][:, None]
        self.diaCoup = self.Step1.diacoup * oldsign[0][:, None]

        for i in range(self.nBound):
            # non-adiabatic couplings between bound and free-electron state
            self.D[i, self.nBound:] = self.grid.kfact * (
                oldsign[0, i] * newsign[1, 0] *
                self.Cross1.Smatrix[i, self.nBound:] -
                oldsign[1, 0] * newsign[0, i] *
                self.Cross2.Smatrix[i, self.nBound:])
            self.Dtemp[i, self.nBound:] = self.D[i, self.nBound:] / \
                                          self.grid.kfact


    def stateOverlapTime(self):
        """
        """

        # Anion properties
        nel_a = self.Step1.anion.nelec_alpha
        nel_b = self.Step1.anion.nelec_beta
        nmo   = self.Step1.anion.nmo

        MOa1 = self.Step1.anion.orbs_alpha
        MOa2 = self.Step2.anion.orbs_alpha
        MOb1 = self.Step1.anion.orbs_beta
        MOb2 = self.Step2.anion.orbs_beta

        # overlap between AO bases of anion at t and t+dt
        Sao = dyA.basis_overlap(self.Step1.anion.basis, self.Step2.anion.basis)
        # overlap between MOs of anion at t and t+dt
        SmoA = la.multi_dot([MOa1.T, Sao, MOa2])
        SmoB = la.multi_dot([MOb1.T, Sao, MOb2])

        over = np.zeros((self.nBound, self.nBound))

        # If only bound ground state is considered
        if self.nBound == 1 and self.state == 0:
            Atemp = np.arange(nel_a)
            Btemp = np.arange(nel_b)

            Alist = np.zeros(nmo, dtype=bool)
            Blist = np.zeros(nmo, dtype=bool)

            for k in Atemp:
                Alist[k] = True
            for k in Btemp:
                Blist[k] = True

            over[0,0] = la.det(SmoA[Alist][:, Alist]) * \
                        la.det(SmoB[Blist][:, Blist])

            return over

        precision = float(config['States']['precision']) / 100
        max_micro = int(config['States']['max_micro'])

        # If excited states are considered
        for i, istate in enumerate(self.states):
            # <i(t)|
            iSDet = self.Step1.anion.cisk[:, istate, :]
            iSDet = iSDet[(-abs(iSDet[:, 0])).argsort()]

            for j, jstate in enumerate(self.states):
                # |j(t+dt)>
                jSDet = self.Step2.anion.cisk[:, jstate, :]
                jSDet = jSDet[(-abs(jSDet[:, 0])).argsort()]

                icutoff = 0.0
                inDet   = 0

                for ii, iexc in enumerate(iSDet):
                    # Construction for microstate of state i at t
                    iAtemp = np.arange(nel_a)
                    iBtemp = np.arange(nel_b)
                    if self.Step1.anion.XAB[int(iexc[3])] == "A":
                        iAtemp[int(iexc[1])] = iexc[2]
                    else:
                        iBtemp[int(iexc[1])] = iexc[2]

                    iAlist = np.zeros(nmo, dtype=bool)
                    iBlist = np.zeros(nmo, dtype=bool)
                    for k in iAtemp:
                        iAlist[k] = True
                    for k in iBtemp:
                        iBlist[k] = True

                    # if True, spin-adaptation is considered for i
                    if (iAlist.astype(int) + iBlist == 1).sum() == 2:
                        iSpin = True
                    else:
                        iSpin = False

                    jcutoff = 0.0
                    jnDet   = 0

                    for jj, jexc in enumerate(jSDet):
                        # Construction for microstate of state i at t+dt
                        jAtemp = np.arange(nel_a)
                        jBtemp = np.arange(nel_b)
                        if self.Step2.anion.XAB[int(jexc[3])] == "A":
                            jAtemp[int(jexc[1])] = jexc[2]
                        else:
                            jBtemp[int(jexc[1])] = jexc[2]

                        jAlist = np.zeros(nmo, dtype=bool)
                        jBlist = np.zeros(nmo, dtype=bool)
                        for k in jAtemp:
                            jAlist[k] = True
                        for k in jBtemp:
                            jBlist[k] = True

                        # if True, spin-adaptation is considered for j
                        if (jAlist.astype(int) + jBlist == 1).sum() == 2:
                            jSpin = True
                        else:
                            jSpin = False

                        # Overlap <i(t)|j(t+dt)>
                        if iSpin and jSpin:
                            # if i and j are open-shell singlets
                            over[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist][:, jAlist]) *
                                la.det(SmoB[iBlist][:, jBlist])
                                +
                                la.det(SmoA[iBlist][:, jAlist]) *
                                la.det(SmoB[iAlist][:, jBlist])
                                +
                                la.det(SmoA[iAlist][:, jBlist]) *
                                la.det(SmoB[iBlist][:, jAlist])
                                +
                                la.det(SmoA[iBlist][:, jBlist]) *
                                la.det(SmoB[iAlist][:, jAlist])) / 2
                        elif iSpin:
                            # if i is an open-shell singlet
                            over[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist][:, jAlist]) *
                                la.det(SmoB[iBlist][:, jBlist])
                                +
                                la.det(SmoA[iBlist][:, jAlist]) *
                                la.det(SmoB[iAlist][:, jBlist])) / np.sqrt(2)
                        elif jSpin:
                            # if j is an open-shell singlet
                            over[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist][:, jAlist]) *
                                la.det(SmoB[iBlist][:, jBlist])
                                +
                                la.det(SmoA[iAlist][:, jBlist]) *
                                la.det(SmoB[iBlist][:, jAlist])) / np.sqrt(2)
                        else:
                            # This is the overlap for all states that do not
                            # involve open-shell singlets, so the "default"
                            over[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist][:, jAlist]) *
                                la.det(SmoB[iBlist][:, jBlist]))

                        jcutoff += jexc[0]**2
                        jnDet   += 1
                        if jcutoff > precision or jnDet == max_micro:
                            break

                    icutoff += iexc[0]**2
                    inDet   += 1
                    if icutoff > precision or inDet == max_micro:
                        break

        return over


    def nonAdBound(self, oldsign):
        """
        calculation of non-adiabatic couplings between bound states

        annotation:
        orbs[i,j] : i = AOs; j = MOs
        """

        newsign = oldsign.copy()

        # Standard terms of non-adiabatic coupling with
        # D_ij = 1/(2*dt) * [<i(t)|j(t+dt)> - <i(t+dt)|j(t)>]
        #      = 1/(2*dt) * [D1[i,j]-D2[i,j]] = 1/(2*dt) * naCoup[i,j]
        D1     = np.zeros((self.nBound, self.nBound))
        D2     = np.zeros((self.nBound, self.nBound))

        # Anion properties
        nel_a = self.Step1.anion.nelec_alpha
        nel_b = self.Step1.anion.nelec_beta
        nmo   = self.Step1.anion.nmo

        MOa1 = self.Step1.anion.orbs_alpha
        MOa2 = self.Step2.anion.orbs_alpha
        MOb1 = self.Step1.anion.orbs_beta
        MOb2 = self.Step2.anion.orbs_beta

        # overlap between AO bases of anion at t and t+dt
        Sao = dyA.basis_overlap(self.Step1.anion.basis, self.Step2.anion.basis)
        # overlap between MOs of anion at t and t+dt
        SmoA = la.multi_dot([MOa1.T, Sao, MOa2])
        SmoB = la.multi_dot([MOb1.T, Sao, MOb2])

        ########################################################################
        # Properties of ionized molecule
        nel_ai = self.Step1.ionized.nelec_alpha
        nel_bi = self.Step1.ionized.nelec_beta
        nmoi   = self.Step1.ionized.nmo

        MOa1i = self.Step1.ionized.orbs_alpha
        MOa2i = self.Step2.ionized.orbs_alpha
        MOb1i = self.Step1.ionized.orbs_beta
        MOb2i = self.Step2.ionized.orbs_beta

        # overlap between AO bases of ionized molecule at t and t+dt
        Saoi = dyA.basis_overlap(self.Step1.ionized.basis,
                                 self.Step2.ionized.basis)
        # overlap between MOs of ionized molecule at t and t+dt
        SmoAi = la.multi_dot([MOa1i.T, Saoi, MOa2i])
        SmoBi = la.multi_dot([MOb1i.T, Saoi, MOb2i])

        Atempi = np.arange(nel_ai)
        Btempi = np.arange(nel_bi)

        Alisti = np.zeros(nmoi, dtype=bool)
        Blisti = np.zeros(nmoi, dtype=bool)

        for k in Atempi:
            Alisti[k] = True
        for k in Btempi:
            Blisti[k] = True

        overi = la.det(SmoAi[Alisti][:, Alisti]) * \
                la.det(SmoBi[Blisti][:, Blisti])

        newsign[1,0] *= np.sign(overi)

        ########################################################################
        # Non-adiabatic couplings
        if self.nBound == 1 and self.state == 0:
            Atemp = np.arange(nel_a)
            Btemp = np.arange(nel_b)

            Alist = np.zeros(nmo, dtype=bool)
            Blist = np.zeros(nmo, dtype=bool)

            for k in Atemp:
                Alist[k] = True
            for k in Btemp:
                Blist[k] = True

            over = la.det(SmoA[Alist][:, Alist]) * \
                   la.det(SmoB[Blist][:, Blist])

            newsign[0,0] *= np.sign(over)

            naCoup  = np.zeros((1, 1))

            return naCoup, newsign

        precision = float(config['States']['precision']) / 100
        max_micro = int(config['States']['max_micro'])

        over = np.zeros((self.nBound, self.nBound))
        for i, istate in enumerate(self.states):
            # <i(t)|
            iSDet1 = self.Step1.anion.cisk[:, istate, :]
            iSDet1 = iSDet1[(-abs(iSDet1[:, 0])).argsort()]
            # <i(t+dt)|
            iSDet2 = self.Step2.anion.cisk[:, istate, :]
            iSDet2 = iSDet2[(-abs(iSDet2[:, 0])).argsort()]

            for j, jstate in enumerate(self.states):
                # |j(t+dt)>
                jSDet2 = self.Step2.anion.cisk[:, jstate, :]
                jSDet2 = jSDet2[(-abs(jSDet2[:, 0])).argsort()]
                # |j(t)>
                jSDet1 = self.Step1.anion.cisk[:, jstate, :]
                jSDet1 = jSDet1[(-abs(jSDet1[:, 0])).argsort()]

                icutoff1 = 0.0
                icutoff2 = 0.0
                inDet    = 0

                for ii, iexc in enumerate(iSDet1):
                    # Construction for microstate of state i at t
                    iAtemp1 = np.arange(nel_a)
                    iBtemp1 = np.arange(nel_b)

                    if self.Step1.anion.XAB[int(iexc[3])] == "A":
                        iAtemp1[int(iexc[1])] = iexc[2]
                    else:
                        iBtemp1[int(iexc[1])] = iexc[2]

                    iAlist1 = np.zeros(nmo, dtype=bool)
                    iBlist1 = np.zeros(nmo, dtype=bool)

                    for k in iAtemp1:
                        iAlist1[k] = True
                    for k in iBtemp1:
                        iBlist1[k] = True

                    # Construction for microstate of state i at t+dt
                    iAtemp2 = np.arange(nel_a)
                    iBtemp2 = np.arange(nel_b)

                    if self.Step2.anion.XAB[int(iSDet2[ii][3])] == "A":
                        iAtemp2[int(iSDet2[ii][1])] = iSDet2[ii][2]
                    else:
                        iBtemp2[int(iSDet2[ii][1])] = iSDet2[ii][2]

                    iAlist2 = np.zeros(nmo, dtype=bool)
                    iBlist2 = np.zeros(nmo, dtype=bool)

                    for k in iAtemp2:
                        iAlist2[k] = True
                    for k in iBtemp2:
                        iBlist2[k] = True

                    # if True, spin-adaptation is considered for i
                    if (iAlist1.astype(int) + iBlist1 == 1).sum() == 2:
                        iSpin = True
                    else:
                        iSpin = False

                    jcutoff1 = 0.0
                    jcutoff2 = 0.0
                    jnDet   = 0

                    for jj, jexc in enumerate(jSDet2):
                        # Construction for microstate of state i at t+dt
                        jAtemp2 = np.arange(nel_a)
                        jBtemp2 = np.arange(nel_b)

                        if self.Step2.anion.XAB[int(jexc[3])] == "A":
                            jAtemp2[int(jexc[1])] = jexc[2]
                        else:
                            jBtemp2[int(jexc[1])] = jexc[2]

                        jAlist2 = np.zeros(nmo, dtype=bool)
                        jBlist2 = np.zeros(nmo, dtype=bool)

                        for k in jAtemp2:
                            jAlist2[k] = True
                        for k in jBtemp2:
                            jBlist2[k] = True

                        # Construction for microstate of state i at t
                        jAtemp1 = np.arange(nel_a)
                        jBtemp1 = np.arange(nel_b)

                        if self.Step1.anion.XAB[int(jSDet1[jj][3])] == "A":
                            jAtemp1[int(jSDet1[jj][1])] = jSDet1[jj][2]
                        else:
                            jBtemp1[int(jSDet1[jj][1])] = jSDet1[jj][2]

                        jAlist1 = np.zeros(nmo, dtype=bool)
                        jBlist1 = np.zeros(nmo, dtype=bool)

                        for k in jAtemp1:
                            jAlist1[k] = True
                        for k in jBtemp1:
                            jBlist1[k] = True

                        # if True, spin-adaptation is considered for j
                        if (jAlist1.astype(int) + jBlist1 == 1).sum() == 2:
                            jSpin = True
                        else:
                            jSpin = False

                        # Overlap <i(t)|j(t+dt)> for signs of wavefunction
                        # spin-adaption not necessary for determination of signs
                        over[i, j] += iexc[0] * jexc[0] * (
                            la.det(SmoA[iAlist1][:, jAlist2]) *
                            la.det(SmoB[iBlist1][:, jBlist2]))

                        # Standard terms of non-adiabatic coupling with
                        # D_ij = 1/(2*dt) * [<i(t)|j(t+dt)> - <i(t+dt)|j(t)>]
                        #      = 1/(2*dt) * [D1 - D2] = 1/(2*dt) * naCoup[i, j]
                        if iSpin and jSpin:
                            # if i and j are open-shell singlets
                            D1[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist1][:, jAlist2]) *
                                la.det(SmoB[iBlist1][:, jBlist2])
                                +
                                la.det(SmoA[iBlist1][:, jAlist2]) *
                                la.det(SmoB[iAlist1][:, jBlist2])
                                +
                                la.det(SmoA[iAlist1][:, jBlist2]) *
                                la.det(SmoB[iBlist1][:, jAlist2])
                                +
                                la.det(SmoA[iBlist1][:, jBlist2]) *
                                la.det(SmoB[iAlist1][:, jAlist2])) / 2

                            D2[i, j] += iSDet2[ii][0] * jSDet1[jj][0] * (
                                la.det(SmoA[jAlist1][:, iAlist2]) *
                                la.det(SmoB[jBlist1][:, iBlist2])
                                +
                                la.det(SmoA[jAlist1][:, iBlist2]) *
                                la.det(SmoB[jBlist1][:, iAlist2])
                                +
                                la.det(SmoA[jBlist1][:, iAlist2]) *
                                la.det(SmoB[jAlist1][:, iBlist2])
                                +
                                la.det(SmoA[jBlist1][:, iBlist2]) *
                                la.det(SmoB[jAlist1][:, iAlist2])) / 2

                        elif iSpin:
                            # if i is an open-shell singlet
                            D1[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist1][:, jAlist2]) * \
                                la.det(SmoB[iBlist1][:, jBlist2])
                                +
                                la.det(SmoA[iBlist1][:, jAlist2]) *
                                la.det(SmoB[iAlist1][:, jBlist2])) / np.sqrt(2)

                            D2[i, j] += iSDet2[ii][0] * jSDet1[jj][0] * (
                                la.det(SmoA[jAlist1][:, iAlist2]) * \
                                la.det(SmoB[jBlist1][:, iBlist2])
                                +
                                la.det(SmoA[jAlist1][:, iBlist2]) * \
                                la.det(SmoB[jBlist1][:, iAlist2])) / np.sqrt(2)

                        elif jSpin:
                            # if j is an open-shell singlet
                            D1[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist1][:, jAlist2]) *
                                la.det(SmoB[iBlist1][:, jBlist2])
                                +
                                la.det(SmoA[iAlist1][:, jBlist2]) *
                                la.det(SmoB[iBlist1][:, jAlist2])) / np.sqrt(2)

                            D2[i, j] += iSDet2[ii][0] * jSDet1[jj][0] * (
                                la.det(SmoA[jAlist1][:, iAlist2]) *
                                la.det(SmoB[jBlist1][:, iBlist2])
                                +
                                la.det(SmoA[jBlist1][:, iAlist2]) *
                                la.det(SmoB[jAlist1][:, iBlist2])) / np.sqrt(2)

                        else:
                            # This is for all states that do not
                            # involve open-shell singlets, so the "default"
                            D1[i, j] += iexc[0] * jexc[0] * (
                                la.det(SmoA[iAlist1][:, jAlist2]) * \
                                la.det(SmoB[iBlist1][:, jBlist2]))

                            D2[i, j] += iSDet2[ii][0] * jSDet1[jj][0] * (
                                la.det(SmoA[jAlist1][:, iAlist2]) * \
                                la.det(SmoB[jBlist1][:, iBlist2]))

                        jcutoff1 += jSDet1[jj][0]**2
                        jcutoff2 += jexc[0]**2
                        jnDet    += 1
                        if (jcutoff1 > precision and jcutoff2 > precision) or \
                            jnDet == max_micro:
                                break

                    icutoff1 += iexc[0]**2
                    icutoff2 += iSDet2[ii][0]**2
                    inDet    += 1
                    if (icutoff1 > precision and icutoff2 > precision) or \
                        inDet == max_micro:
                            break

        newsign[0] *= np.sign(
                      np.array([sorted(i, key=abs) for i in over])[:,-1])
        naCoup  = D1 * newsign[0][None, :] * oldsign[0][:, None] - \
                  D2 * newsign[0][:, None] * oldsign[0][None, :]

        return naCoup, newsign


    def calcProb(self, c, c0, fo):
        """
        Calculation of hopping probabilities of single time step
        """

        prob      = np.zeros(self.nstates+1, dtype=float)
        c0sq      = np.zeros(self.nstates+1, dtype=float)
        csq       = np.zeros(self.nstates+1, dtype=float)
        c0sq[:-1] = abs(c0)**2
        csq[:-1]  = abs(c)**2
        c0sq[-1]  = 1 - np.sum(abs(c0)**2)
        csq[-1]   = 1 - np.sum(abs(c)**2)

        if c0sq[self.state] > csq[self.state]:
            allowed = csq > c0sq
            akk     = np.einsum("i,i->", csq-c0sq, allowed)
            prob    = -1.0 * allowed * (csq[self.state] - c0sq[self.state]) / \
                                       (akk * c0sq[self.state]) * (csq - c0sq)

        prob[self.state] = 1.0 - sum(prob)

        fo.write("\nHopping Probabilities:\n")
        for i in range(self.nBound):
            fo.write("  State %2i      : %12.9f\n"%(self.states[i], prob[i]))
        fo.write("  Ionization    : %12.9f\n"%np.sum(prob[self.nBound:-1]))
        fo.write("  Neg-VDE Decay : %12.9f\n"%prob[-1])

        return prob


    def singleWaveHop(self, rNr, trajpop, step, Ekin, fo):
        """
        Check if hopping criteria are met
        Hopping is evaluated for single plane waves
        """

        fo.write("\nSINGLE WAVE HOPPING\n")

        for j, randNr in enumerate(rNr):
            probsum = 0.0
            for i, prob in enumerate(self.prob):
                if (prob + probsum) > randNr:
                    if self.nBound <= i < self.nstates:
                        kvec    = self.grid.kbasis[i - self.nBound]
                        Eel     = la.norm(kvec)**2 / 2
                        Efin    = self.Step2.Hdiag[i]
                        trajpop = self.checkEnergy(trajpop, step, Ekin,
                                                   Eel, kvec, Efin, "", fo)
                    elif i == self.nstates:
                        trajpop = self.checkEnRes(trajpop, step, Ekin, "", fo)
                    break
                else:
                    probsum += prob

        return trajpop


    def wavepacketHop(self, rNr, trajpop, step, Ekin, fo):
        """
        Check if hopping criteria are met
        Hopping is evaluated for a weighted average of plane waves
        """

        fo.write("\nWAVEPACKET HOPPING\n")

        Efin = np.einsum("k,k->", self.Step2.Hdiag[self.nBound:],
                         abs(self.c[self.nBound:])**2) / \
               sum(abs(self.c[self.nBound:])**2)

        kvec = np.einsum("kx,k->x", self.grid.kbasis,
                         abs(self.c[self.nBound:])**2) / \
               sum(abs(self.c[self.nBound:])**2)

        Eel  = (la.norm(self.grid.kbasis, axis=1) * \
                abs(self.c[self.nBound:]))**2
        Eel  = sum(Eel) / (2 * sum(abs(self.c[self.nBound:])**2))

        prob = np.zeros(self.nBound+2)
        prob[:self.nBound] = self.prob[:self.nBound]
        prob[-2] = sum(self.prob[self.nBound:-1])
        prob[-1] = self.prob[-1]

        for j, randNr in enumerate(rNr):
            probsum = 0.0
            for i in range(self.nBound+2):
                if (prob[i] + probsum) > randNr:
                    # only hopping into free-electron states is allowed here
                    if i == self.nBound:
                        trajpop = self.checkEnergy(trajpop, step, Ekin,
                                                   Eel, kvec, Efin, "WP", fo)
                    elif i == self.nBound+1:
                        trajpop = self.checkEnRes(trajpop, step, Ekin, "WP", fo)
                    break
                else:
                    probsum += prob[i]

        return trajpop


    def boundHop(self, rNr, Ekin, fo):
        """
        Check for bound state hop
        """
        totprob = np.sum(self.prob[:self.nBound])
        prob    = self.prob[:self.nBound]/totprob

        temp  = "\nBOUND STATE HOPPING\n\n"
        temp += "State probabilities for bound state hop: \n"
        for i in range(self.nBound):
            temp += "  State %i: %12.9f\n"%(self.states[i], prob[i])
        temp += "Random number: %12.9f\n\n"%rNr

        probsum = 0.0
        for i in range(self.nBound):
            if (prob[i] + probsum) > rNr:
                if i == self.state:
                    temp += "Staying in current state\n\n"
                else:
                    Estart = self.Step2.Hdiag[self.state]
                    Efin   = self.Step2.Hdiag[i]
                    temp  += "E[kin]          : %12.9f\n"%Ekin
                    temp  += "E[start]        : %12.9f\n"%Estart
                    temp  += "E[start]+E[kin] : %12.9f\n"%(Estart+Ekin)
                    temp  += "E[end]          : %12.9f\n"%Efin
                    if (Estart + Ekin) > Efin:
                        Ekin = self.rescaleKinEn(Ekin, i)
                        self.state   = i
                        self.statenr = self.states[i]
                        temp += "\nHop to bound state %i\n\n"%self.statenr
                    else:
                        temp += "\nNo bound state hop: state %i "%self.statenr
                        temp += "energically not accessable\n\n"
                break
            else:
                probsum += prob[i]

        fo.write(temp)

        return Ekin


    def checkEnergy(self, trajpop, step, Ekin, Eel, kvec, Efin, wp, fo):
        """
        if hop is allowed by probability, checks if energy criteria
        are fulfilled
        """

        Estart = self.Step2.Hdiag[self.state]

        if Estart + Ekin - Efin >= 0:
            trajpop -= 1
            temp     = "\nHop into continuum! (Option %s) %s\n\n"%(
                       self.index, wp)
        else:
            temp     = "\nNo hop due to insufficient energy! %s\n\n"%wp

        temp += "E[kin]           : %12.9f\n"%Ekin
        temp += "E[start]         : %12.9f\n"%Estart
        temp += "E[start]+E[kin]  : %12.9f\n"%(Estart+Ekin)
        temp += "E[free]          : %12.9f\n"%Efin

        if Estart + Ekin - Efin >= 0:
            temp += "\nHopped at step %i with k-vector: "%step
            temp += "%8.5f %8.5f %8.5f  "%(kvec[0], kvec[1], kvec[2])
            temp += "energy gap neutral/anion: "
            temp += "%e eV\n"%(27.211386 * (Efin - Estart))

            with open("freeEl%s%s.dat"%(self.index, wp), "a") as f:
                f.write("%.2f %.5e %.5e %.5e"%(
                        step * self.dtfs, kvec[0], kvec[1], kvec[2]))
                f.write(" %12.9f %12.9f\n"%(Eel, 27.211386*(Efin - Estart)))

        else:
            temp += "k-vector: %2.5f %2.5f %2.5f\n"%(kvec[0], kvec[1],
                                                       kvec[2])

        temp += "E[photoelectron] : %12.9f\n"%Eel
        temp += "E[kin,neutral]   : %12.9f\n\n"%(Estart+Ekin-Efin)
        fo.write(temp)

        return trajpop


    def checkEnRes(self, trajpop, step, Ekin, wp, fo):
        """
        if hop is allowed by probability, checks if energy criteria
        for decay due to negative VDE if fulfilled
        """
        Estart = self.Step2.Hdiag[self.state]
        neuEn  = self.Step2.ionized.total_energy-self.refEn

        if Estart > neuEn:
            trajpop -= 1
            temp     = "\nHop into continuum! (Option %s) %s\n\n"%(self.index,
                                                                    wp)
        else:
            temp     = "\nNo hop due to insufficient energy!\n\n"

        temp += "E[kin]           : %2.9f\n"%Ekin
        temp += "E[start]         : %2.9f\n"%Estart
        temp += "E[start]+E[kin]  : %2.9f\n"%(Estart+Ekin)
        temp += "E[free]          : %2.9f\n"%neuEn

        if Estart > neuEn:
            temp += "\nHopped at step %i "%step
            temp += "due to limited resonance lifetime with "
            temp += "VDE neutral/anion: "
            temp += "%e eV\n"%(27.211386*(neuEn - Estart))

            with open("freeEl%s%s.dat"%(self.index, wp), "a") as f:
                f.write("%.2f %.5e %.5e %.5e"%(
                        step * self.dtfs, 0.0, 0.0, 0.0))
                f.write(" %.9f %.9f\n"%(
                        Estart - neuEn, 27.211386 * (neuEn - Estart)))

            temp += "E[photoelectron] : %2.9f\n"%(Estart-neuEn)
            temp += "E[kin,neutral]   : %2.9f\n\n"%Ekin

        fo.write(temp)

        return trajpop


    def rescaleKinEn(self, kinEold, newstate):
        """
        Rescales kinetic energy after bound state hop
        """

        Enew, Eold = self.Step2.Hdiag[newstate], self.Step2.Hdiag[self.state]

        return kinEold - (Enew - Eold)


    def shiftSteps(self):
        """
        overwrites Step1 instance with Step2 instance and
        "old" (Step1) with "new" (Step2) hopping probabilities
        """

        self.Step1   = self.Step2
        self.oldprob = self.prob


    def writeDat(self, step):
        """
        Population dynamics output for current step is written
        """

        t        = self.dtfs * step
        H, D     = self.Step2.Hdiag, self.D
        totCoup  = -1.0j*self.diaCoup[self.state] - D[self.state]/(2*self.dt)
        totCoup  = -1.0j*self.diaCoup - D/(2*self.dt)
        l        = len(H)
        apl      = int(config['Dynamics']['printlevel'])
        try:
            skipstep = int(config['Dynamics']['skipstep'])
        except KeyError:
            skipstep = 0

        if ((apl <= -1) and (step%skipstep == 0)) or apl >= 0:

            if config['Dynamics']['enerprint'] == 'full':
                with open("energies%s.dat"%self.index, "a") as out:
                    out.write("%7.2f %12.9e "%(t, H[self.state]))
                    out.write(l*"%12.9e "%tuple(i for i in H)+"\n")
            elif config['Dynamics']['enerprint'] == 'none':
                neuEn = H[self.nBound]-self.grid.kinEn[0]
                with open("energies%s.dat"%self.index, "a") as out:
                    out.write("%7.2f %12.9e "%(t, H[self.state]))
                    for i in range(self.nBound):
                        out.write("%12.9e "%H[i])
                    out.write("%12.9e\n"%neuEn)

            if config['Dynamics']['coupprint'] == 'full':
                with open("couplings%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t + l * " %9.5e "%tuple(
                        abs(totCoup[self.state,i]) for i in range(l)) + "\n")
                for j in range(self.nBound):
                    if self.nBound == 1 and self.statenr == 0:
                        break
                    coupName = "couplings%s_to_%i.dat"%(self.index,
                                                        self.states[j])
                    with open(coupName, "a") as out:
                        out.write("%7.2f "%t + l * " %9.5e "%tuple(
                            abs(totCoup[j, i]) for i in range(l)) + "\n")
            elif config['Dynamics']['coupprint'] == 'summed':
                with open("couplings%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t)
                    for i in range(self.nBound):
                        out.write("%9.5e "%abs(totCoup[self.state, i]))
                    out.write("%9.5e\n"%sum(abs(totCoup[self.state,
                                                        self.nBound:])))
                for j in range(self.nBound):
                    if self.nBound == 1 and self.statenr == 0:
                        break
                    coupName = "couplings%s_to_%i.dat"%(self.index,
                                                        self.states[j])
                    with open(coupName, "a") as out:
                        out.write("%7.2f "%t)
                        for i in range(self.nBound):
                            out.write("%9.5e "%abs(totCoup[j, i]))
                        out.write("%9.5e\n"%sum(abs(totCoup[j, self.nBound:])))

            if config['Dynamics']['coefprint'] == 'full':
                with open("coefficients%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t + l * " %9.5e "%tuple(
                        abs(i)**2 for i in self.c) + "\n")
            elif config['Dynamics']['coefprint'] == 'summed':
                with open("coefficients%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t)
                    for i in range(self.nBound):
                        out.write("%9.5e "%(abs(self.c[i])**2))
                    out.write("%9.5e\n"%sum(abs(self.c[self.nBound:])**2))

            if config['Dynamics']['nadiaprint'] == 'full':
                with open("nonAd%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t+l*" %9.5e "%tuple(
                        abs(-self.Dtemp[self.state, i] / (2*self.dt))
                        for i in range(l))+"\n")
                with open("diaC%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t + l*" %9.5e "%tuple(
                        abs(-1.0j*self.diaTemp[self.state, i])
                        for i in range(l)) + "\n")
            elif config['Dynamics']['nadiaprint'] == 'summed':
                with open("nonAd%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t)
                    for i in range(self.nBound):
                        out.write("%9.5e "%
                                  abs(self.Dtemp[self.state, i] / (2*self.dt)))
                    out.write("%9.5e\n"%
                              sum(abs(self.Dtemp[self.state, self.nBound:] / \
                                      (2*self.dt))))
                with open("diaC%s.dat"%self.index, "a") as out:
                    out.write("%7.2f "%t)
                    for i in range(self.nBound):
                        out.write("%9.5e "%abs(1.0j *
                                               self.diaTemp[self.state, i]))
                    out.write("%9.5e\n"%
                              sum(abs(1.0j * \
                                      self.diaTemp[self.state, self.nBound:])))

        with open("trajpop%s.dat"%self.index, "a") as f:
            f.write("%.2f %i\n"%(t, self.trajpop[0]))
        with open("trajpopWP%s.dat"%self.index, "a") as f:
            f.write("%.2f %i\n"%(t, self.trajpop[1]))
        with open("boundPop%s.dat"%self.index, "a") as out:
            out.write("%.2f %i\n"%(t, self.statenr))

        with open("tempcoef%s.dat"%self.index, "w") as out:
            out.write(l * "%9.5e "%(tuple(i.real for i in self.c)) + "\n")
            out.write(l * "%9.5e "%(tuple(i.imag for i in self.c)) + "\n")

        if config['Dynamics']['probprint'] == 'full':
            with open("probabilities%s.dat"%self.index, "a") as out:
                out.write("%7.2f "%t + (l+1) * " %9.5e "%tuple(
                    i for i in self.prob) + "\n")
        elif config['Dynamics']['probprint'] == 'summed':
            with open("probabilities%s.dat"%self.index, "a") as out:
                out.write("%7.2f "%t)
                for i in range(self.nBound):
                    out.write("%9.5e "%(self.prob[i]))
                out.write("%9.5e "%(sum(abs(self.prob[self.nBound:-1]))))
                out.write("%9.5e\n"%(self.prob[-1]))

        if apl >= 0:
            with open("norm%s.dat"%self.index, "a") as out:
                out.write("%7.2f "%t + "%5.9f\n"%sum([abs(i)**2 for i
                                                      in self.c]))
                logging.info("Norm after integration\n%f"%sum([abs(i)**2 for i
                                                               in self.c]))


    def writeCube(self, step):
        """
        NOT IMPLEMENTED YET
        Writes .cube file for ionized electron
        """

        coord = self.Step2.coord

        x = np.arange(-10.0, 11.0, 1.0)

        cube = np.zeros((len(x), len(x), len(x)), dtype=complex)
        for ki in range(self.grid.points):
            k     = self.grid.kbasis[ki]
            cube += abs(self.c[self.nBound+ki])**2 * np.exp(1.0j*(
                        k[0] * x[:, None, None] +
                        k[1] * x[None, :, None] +
                        k[2] * x[None, None, :]))

        f = open("WP_step%i.cube"%step, "w")
        f.write("Step: %i\nTime: %.2f\n"%(step, step * self.dtfs))
        f.write("%i -10.0 -10.0 -10.0\n"%(self.Step2.nat))
        f.write("11 2.0 0.0 0.0\n")
        f.write("11 0.0 2.0 0.0\n")
        f.write("11 0.0 0.0 2.0\n")
        for i in range(self.Step2.nat):
            f.write("%s %.7f %.7f %.7f\n"%(
                    self.Step2.atnr[i],
                    coord[i, 0], coord[i, 1], coord[i, 2]))

        f.close()
        print(cube)
        quit("cube")

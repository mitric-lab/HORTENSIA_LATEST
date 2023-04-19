#!/usr/bin/env python3

import os
import ast
from datetime import datetime
from configparser import ConfigParser

import numpy as np
import numpy.linalg as la

import hortensia_latest.misc as misc
import hortensia_latest.gridClass as grid
import hortensia_latest.populationDynamics as pop

config = ConfigParser()
config.read('config.ini')

if config['2e-Integrals']['exact2elInt'] == 'True':
    import hortensia_latest.exactIntegral as exI

if config['QC']['qcmethod'].lower() in ['g09', 'g16']:
    import hortensia_latest.handlerG09 as qc
elif config['QC']['qcmethod'].lower() == 'qchem':
    import hortensia_latest.handlerQChem as qc
else:
    quit("Quantum chemistry method not understood!")

################################################################################
#                                                                              #
#                     Nuclear Dynamics Calculation Class                       #
#                                                                              #
################################################################################


class NuclearDynamics:
    """
    Class that handles the aspects of the classical nuclear dynamics including
    the Velocity-Verlet algorithm and elimination of translation and rotation,
    as well as output generation of everything related to nuclear dynamics.
    Furthermore, instances of the PopulationDynamics class are created and
    methods for the calculation of population dynamics are called here.
    """

    def __init__(self, infile):
        """
        initialization of nuclear dynamics
        """

        self.restart = config['Dynamics']['restart'] == 'True'

        self.states = ast.literal_eval(
            config['States']['anion_states'])
        if self.restart:
            self.statenr = int(np.loadtxt("boundPop.dat")[-1, -1])
            self.state = self.states.index(self.statenr)
        else:
            self.state = int(config['States']['starting_state'])
            self.statenr = self.states[self.state]

        if len(self.states) > 1 or self.statenr != 0:
            self.excited = True
        else:
            self.excited = False

        self.Eshift   = float(config['States']['Eshift'])

        self.steps    = int(config['Dynamics']['steps'])
        self.dt       = float(config['Dynamics']['timestep']) / misc.autime2fs

        self.charge   = int(config['QC']['charge'])
        self.mult     = int(config['QC']['mult'])
        self.qcmethod = config['QC']['qcmethod']
        self.func     = config['QC']['func']
        self.Nproc    = int(config['QC']['Nproc'])
        self.conv     = int(config['QC']['convergence'])
        self.maxcyc   = int(config['QC']['maxiter'])

        self.cartBFs  = []

        with open("basis", "r") as f:
            basis = f.readlines()

        ltype = {"S":0, "P":1, "D":2, "F":3,
                 0:"S", 1:"P", 2:"D", 3:"F"}

        newbas = []
        atoms  = []

        # Rewrites the basis so that e.g. 'SP' basis functions are separated
        # into S and P. This makes the output of qc programs more consistent
        at = -1
        for i, j in enumerate(basis):
            if len(j.split()) == 2 and j.split()[1] == "0":
                newbas.append([])
                atoms.append(j.split()[0])
                at += 1
                nf  = [0, 0, 0, 0]
                continue

            if len(j.split()) == 3:
                try:
                    float(j.split()[0].replace("D", "E").replace("d", "E"))
                    continue
                except:
                    pass

                for w, q in enumerate(j.split()[0]):
                    if ltype[q.upper()]+1 > len(newbas[at]):
                        newbas[at].append([])
                    nfunc = int(j.split()[1])
                    newbas[at][ltype[q]].append([])
                    for k in range(nfunc):
                        exp = float(basis[i+k+1].split()[0].replace(
                                    "D", "E").replace("d", "E"))
                        fac = float(basis[i+k+1].split()[1+w].replace(
                                    "D", "E").replace("d", "E"))
                        newbas[at][ltype[q]][nf[ltype[q]]].append([exp, fac])
                    nf[ltype[q]] += 1

        os.system("mv basis basis.old")

        with open("basis", "w") as f:
            for i, j in enumerate(newbas):
                f.write("%s  0\n"%atoms[i])
                for ki, k in enumerate(j):
                    for l in k:
                        f.write("%s %i 1.0\n"%(ltype[ki], len(l)))
                        for m in l:
                            f.write("  %.9f %.9f\n"%(m[0], m[1]))
                f.write("****\n")
            f.write("\n")

        # Looks for higher basis functions than P,
        # so cartesian calculation is requested
        with open("basis", "r") as f:
            basis = f.readlines()
        for i in basis:
            if len(i.split()) == 3:
                if i.split()[0] == "d" or i.split()[0] == "D":
                    self.cartBFs.append("6D ")
                    break
        for i in basis:
            if len(i.split()) == 3:
                if i.split()[0] == "f" or i.split()[0] == "F":
                    self.cartBFs.append("10F ")
                    break
        for i in basis:
            if len(i.split()) == 3:
                if i.split()[0] == "g" or i.split()[0] == "G":
                    self.cartBFs.append("15G ")
                    break

        self.dynPrep(infile)

        with open("out.out", "a") as f:
            f.write("\nREFERENCE ENERGY: %4.15f\n"%self.refEn)


    def dynPrep(self, infile):
        """
        Clears the folder of any unnecessary data and
        writes the first step of the dynamics into the output files
        """

        data  = ""
        if not self.restart:
            data += "dynamics.xyz *.dat "
        data += "Gau* "
        data += "neu.* an.* "              # !!!
        os.system("rm -f %s"%data)

        with open(infile, "r") as f:
            lines = f.readlines()

        self.nat          = int(lines[0])
        self.atnr         = np.zeros(self.nat, dtype=int)
        self.atS          = np.zeros(self.nat, dtype=str)
        self.coord        = np.zeros((self.nat, 3))
        self.coord_angs   = np.zeros((self.nat, 3))
        self.velo         = np.zeros((self.nat, 3))

        if self.restart:
            with open("dynamics.xyz", "r") as f:
                xyzdata = f.readlines()
            self.refEn  = float(xyzdata[1].split()[2])

            with open("velocities.dat", "r") as f:
                velodata   = f.readlines()
            self.startstep = int(velodata[len(velodata)-2-self.nat].split()[1])

            fact = misc.bohr_to_angs
            for i in range(self.nat):
                j                  = self.nat - 1 - i
                s, x, y, z         = xyzdata[len(xyzdata)-1-i].split()
                self.atS[j]        = s
                self.atnr[j]       = misc.NrToS(s)
                self.coord[j]      = (float(x) / fact, float(y) / fact,
                                      float(z) / fact)
                self.coord_angs[j] = (float(x), float(y), float(z))

                vx, vy, vz   = velodata[len(velodata)-1-i].split()
                self.velo[j] = (float(vx), float(vy), float(vz))

        else:
            self.startstep = 0

            with open(infile, "r") as f:
                lines = f.readlines()

            fact  = misc.bohr_to_angs
            for i in range(self.nat):
                s, x, y, z         = lines[i+1].split()
                self.atS[i]        = s
                self.atnr[i]       = misc.NrToS(s)
                self.coord_angs[i] = (fact * float(x), fact * float(y),
                                      fact * float(z))
                self.coord[i]      = (float(x), float(y), float(z))

                vx, vy, vz         = lines[i+self.nat+1].split()
                self.velo[i]       = (float(vx), float(vy), float(vz))

        self.masses = np.zeros(self.nat)
        for i, j in enumerate(self.atnr):
            self.masses[i] = misc.masses(misc.NrToS(j))

        self.eliminateTransRot()

        self.kinEn = self.kinEnergy(self.velo)

        fo = open("out.out", "a")

        # QC calculation of initial step
        stateSet = [self.states, self.state, self.statenr, self.excited]
        calcSet  = [self.Nproc , self.func , self.cartBFs, self.conv,
                    self.maxcyc, self.qcmethod]
        molSet   = [self.coord.copy(), self.atS, self.charge, self.mult]

        self.acc = qc.externalQC(stateSet, calcSet, molSet, "1", fo)

        # Generation of PopulationDynamics instance
        try:
            kDistr = ast.literal_eval(config['Continuum']['kDistr'])
        except ValueError:
            kDistr = config['Continuum']['kDistr']

        if isinstance(kDistr, list):
            self.popDyn = self.genPopDyn(fo, kDistr, calcSet)
        else:
            if kDistr == 'snub':
                options = [float(config['Continuum']['maxEk']),
                           int(config['Continuum']['nEk'])]
            elif kDistr == 'fib':
                options = [float(config['Continuum']['maxEk']),
                           int(config['Continuum']['nEk']),
                           int(config['Continuum']['nkfib'])]
            elif kDistr == 'cubic':
                maxEk   = [float(config['Continuum']['maxkx']),
                           float(config['Continuum']['maxky']),
                           float(config['Continuum']['maxkz'])]
                nk      = [int(config['Continuum']['nkx']),
                           int(config['Continuum']['nky']),
                           int(config['Continuum']['nkz'])]
                options = [maxEk, nk]

            self.grid = grid.Grid(kDistr, options)
            self.popDyn = [pop.PopulationDynamics(self.states, self.state,
                           self.Eshift, self.grid, self.coord, self.atS, 
                           calcSet, fo)]

        if self.restart:
            for i, popDyn in enumerate(self.popDyn):
                popDyn.refEn = self.refEn
                popDyn.Step1.herm_wfn(fo, popDyn.preFT, refEn=self.refEn)
                if config['2e-Integrals']['exact2elInt'] == 'True':
                    exInt = exI.ExactInt(self.atS, self.grid)
                    popDyn.Step1.calcHmatrix(self.refEn, popDyn.preMO, exInt)
                else:
                    popDyn.Step1.calcHmatrix(self.refEn, popDyn.preMO)
        else:
            self.refEn = self.popDyn[0].refEn

            with open("dynamics.xyz", "w") as f:
                f.write("%i\nREFERENCE ENERGY: %4.15f\n"%(
                        self.nat, self.refEn))
                for i, coord in enumerate(self.coord_angs):
                    f.write("%s %16.12f %16.12f %16.12f\n"%(
                            self.atS[i], coord[0], coord[1], coord[2]))

            with open("velocities.dat", "w") as f:
                f.write("Step: 0\n\n")
                for i in self.velo:
                    f.write("%16.12f %16.12f %16.12f\n"%tuple(i))

        fo.close()


    def genPopDyn(self, fo, kDistr, calcSet):
        """
        generates instances of PopulationDynamics class,
        used when multiple k-grids are propagated simulaneously
        """

        popDyn = []
        js     = 0
        jc     = 0
        jfib   = 0

        f = open("info.dat", "w")
        for i, ktype in enumerate(kDistr):
            f.write("Specifications for free electron basis, Option %i\n\n"%i)

            if ktype == 'snub':
                nEk   = ast.literal_eval(config['Continuum']['nEk'])
                maxEk = ast.literal_eval(config['Continuum']['maxEk'])

                options = [maxEk[js], nEk[js]]

                f.write("Distribution   : spherical (Snub cube)\n")
                f.write("Plane waves    : 24 per energy\n")
                f.write("               : %i different energies\n"%nEk[js])
                f.write("               : %i total\n"%(nEk[js]*24))
                f.write("Maximum energy : %f eV\n"%maxEk[js])

                js     += 1

            elif ktype == 'fib':
                nEk   = ast.literal_eval(config['Continuum']['nEk'])
                maxEk = ast.literal_eval(config['Continuum']['maxEk'])
                nkfib = ast.literal_eval(config['Continuum']['nkfib'])
                nEk   = nEk[js]
                maxEk = maxEk[js]
                nkfib = nkfib[jfib]

                options = [maxEk, nEk, nkfib]

                f.write("Distribution   : spherical (Fibonacci)\n")
                f.write("Plane waves    : %i per energy\n"%nkfib)
                f.write("               : %i different energies\n"%nEk)
                f.write("               : %i total\n"%(nEk*nkfib))
                f.write("Maximum energy : %f eV\n"%maxEk)

                js     += 1
                jfib   += 1

            elif ktype == 'cubic':
                maxkx = ast.literal_eval(config['Continuum']['maxkx'])
                maxky = ast.literal_eval(config['Continuum']['maxky'])
                maxkz = ast.literal_eval(config['Continuum']['maxkz'])
                nkx = ast.literal_eval(config['Continuum']['nkx'])
                nky = ast.literal_eval(config['Continuum']['nky'])
                nkz = ast.literal_eval(config['Continuum']['nkz'])

                maxEk   = [maxkx[jc], maxky[jc], maxkz[jc]]
                nk      = [nkx[jc], nky[jc], nkz[jc]]
                options = [maxEk, nk]

                Emax = la.norm(maxEk)**2 / 2 * misc.hartree_to_eV

                f.write("Distribution   : cubic\n")
                f.write("Plane waves    : %i in kx\n"%nk[0])
                f.write("               : %i in ky\n"%nk[1])
                f.write("               : %i in kz\n"%nk[2])
                f.write("               : %i total\n"%(nk[0]*nk[1]*nk[2]))
                f.write("Maximum energy : %f eV\n"%Emax)

                jc      += 1
            f.write("\n\n")

            tempgrid   = grid.Grid(ktype, options)
            temppopdyn = pop.PopulationDynamics(self.states, self.state,
                                                self.Eshift, tempgrid,
                                                self.coord, self.atS, calcSet,
                                                fo, index=i)
            popDyn.append(temppopdyn)
        f.close()

        return popDyn


    def simulation(self):
        """
        calculates the whole dynamics
        """

        fo = open("out.out", "a")
        fo.write("\n\n---------------\nDYNAMICS STARTS\n---------------\n")

        for i in range(self.startstep+1, self.steps+1):
            # Print of time step and actual time
            print("\nStep: %s / %s\n" %
                  (str(i).zfill(len(str(self.steps))),
                   str(self.steps).zfill(len(str(self.steps)))))
            fo.write("\nStep: %i\n" % i)
            fo.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S')+"\n\n")

            # Call the verlet method to propagate the nuclear dynamics
            self.verlet(fo)

            # With the new nuclear dynamics step, propagate the electronic
            # state dynamics and evaluate (surface) hopping
            trajcheck = 0
            for j, popDyn in enumerate(self.popDyn):
                intRet = popDyn.integrate(i, self.coordold, self.coord,
                                          self.kinEn, fo)
                self.state   = intRet[0]
                self.statenr = intRet[1]
                self.velo   *= intRet[2]/self.kinEn
                self.kinEn   = 1.0*intRet[2]

                trajcheck += self.popDyn[j].trajpop[0]

            self.writeXyzVelo(i)

            fo.flush()

            # if all trajs have bound population of 0, finishes calculation
            if trajcheck == 0:
                break

        self.finishDynamics(fo)


    def verlet(self, fo):
        """
        propagates nuclear coordinates and starts qc calculation
        """

        self.accold     = self.acc.copy()
        self.veloold    = self.velo.copy()
        self.kinEnold   = self.kinEn.copy()
        self.coordold   = self.coord.copy()
        # Velocity-Verlet algorithm for nuclear coordinates
        self.coord      = self.coord + (self.dt * self.velo +
                                        0.5 * self.dt**2 * self.acc)
        self.coord_angs = self.coord * misc.bohr_to_angs

        stateSet = [self.states, self.state, self.statenr, self.excited]
        calcSet  = [self.Nproc, self.func, self.cartBFs, self.conv,
                    self.maxcyc, self.qcmethod]
        molSet   = [self.coord.copy(), self.atS, self.charge, self.mult]

        # Call external qc program for calculation of electronic properties
        self.acc   = qc.externalQC(stateSet, calcSet, molSet, "2", fo)

        # Calculation of new velocities with forces from qc output
        self.velo  = self.veloold + 0.5 * self.dt * (self.acc + self.accold)
        # Elimination of translation and rotation
        self.eliminateTransRot()
        # Calculation of new kinetic energy
        self.kinEn = self.kinEnergy(self.velo)

        stateEn    = qc.stateEnergy(self.statenr, max(self.states))

        totOld     = self.popDyn[0].Step1.Hdiag[self.state] + self.kinEnold
        totNew     = stateEn - self.refEn + self.kinEn
        correction = 1 + (totOld - totNew) / self.kinEn
        self.kinEn = self.kinEn * correction

        fo.write("energy conservation error in kinEn: %.9f\n"%correction)
        fo.write("kinetic energy rescaled\n")


    def eliminateTransRot(self):
        """
        Eliminates translation and rotation, therefore correcting velocities
        """

        # translation
        pTotal = np.sum(self.velo * self.masses[:, None], axis=0)
        mTotal = np.sum(self.masses)
        self.velo -= pTotal[None, :]/mTotal

        # angular momentum
        L = np.sum(self.masses[:, None] *
                   np.cross(self.coord, self.velo, axis=1), axis=0)
        # moment of inertia
        I = np.zeros((3, 3))
        Emat = np.eye(3)
        for i in range(self.nat):
            I += self.masses[i]*(np.inner(self.coord[i], self.coord[i])*Emat -
                                 np.outer(self.coord[i], self.coord[i]))
        # angular velocity
        angVel = np.dot(la.inv(I), (L))

        # rotation
        self.velo -= np.cross(angVel[None, :], self.coord, axis=1)


    def kinEnergy(self, velo):
        """
        Calculates the kinetic energy of a molecule
        """

        kinE = 0
        for i in range(self.nat):
            kinE += 0.5 * self.masses[i] * la.norm(velo[i])**2

        return kinE


    def writeXyzVelo(self, nucl_step):
        """
        Writes nuclear coordinates and nuclear velocities to output files
        """

        with open("dynamics.xyz", "a") as f:
            f.write(str(self.nat)+"\n\n")
            for i, coord in enumerate(self.coord_angs):
                f.write("%s "%self.atS[i]
                        + 3 * "%16.12f "%tuple(coord) + "\n")

        with open("velocities.dat", "a") as f:
            f.write("Step: %i\n\n"%nucl_step)
            for i in self.velo:
                f.write(3*"%16.12f "%tuple(i)+"\n")


    def finishDynamics(self, fo):
        fo.write("\n" + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                 + "\nDynamics finished\n")
        fo.close()

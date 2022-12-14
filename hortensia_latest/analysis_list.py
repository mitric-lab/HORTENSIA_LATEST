#!/usr/bin/env python3

import os
import ast
import subprocess as sp
from configparser import ConfigParser

import numpy as np
import numpy.linalg as la

config = ConfigParser()
config.read('config.ini')
kDistr = ast.literal_eval(config['Continuum']['kDistr'])

################################################################################
#                                                                              #
#                       pre- and post-dynamics features                        #
#                                                                              #
################################################################################


def trajCheck(trajs, wavepacket):
    """
    Checks for erroneous trajectories and moves them to .error,
    evaluates the total hopping statistics and writes pop.dat files with the
    total population of every time step
    !!! This function was not properly tested !!!
    """

    out = "out.out"

    states = len(ast.literal_eval(config['States']['anion_states']))
    dynl   = len(kDistr)
    nrpertraj = int(config['Dynamics']['nrpertraj'])
    ntrajs = 0
    nerror = 0
    nhop   = np.zeros(dynl, dtype=int)
    nnohop = np.zeros(dynl, dtype=int)

    pop    = np.zeros((dynl, int(config['Dynamics']['steps']) + 1, states + 1))

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    for i in range(trajs):
        if os.path.exists("TRAJ_%i.out"%i):
            ntrajs       += 1
            cmd           = "grep 'Dynamics finished' TRAJ_%i.out/%s"%(i,out)
            grepout       = sp.Popen([cmd], shell=True, stdout=sp.PIPE,
                            stderr=sp.STDOUT)
            stdout,stderr = grepout.communicate()

            if stdout == b"":
                print("\nError: Trajectory %i did not finish normally"%i)
                print("TRAJ_%i.out moved to TRAJ_%i.error"%(i,i))
                os.system("mv TRAJ_%i.out TRAJ_%i.error"%(i,i))

                nerror += 1
                continue

            print("\nTraj %i ended normally"%i)

            if os.path.exists("TRAJ_%i"%i):
                os.system("rm -r TRAJ_%i"%i)

            for j in range(dynl):
                bound = np.loadtxt("TRAJ_%i.out/boundPop%i.dat"%(i,j))
                temp  = np.loadtxt("TRAJ_%i.out/trajpop%s%i.dat"%(i,wp,j))
                bound = np.asarray(bound[:,1], dtype=int)
                temp  = np.asarray(temp[:,1] , dtype=int)
                print("Option %i: %5i bound, %5i hopped"%(j, temp[-1,1],
                                                          nrpertraj-temp[-1,1]))
                for l,k in enumerate(temp):
                    pop[j,l,bound[l]] += k
                    pop[j,l,-1]       += nrpertraj - k
                pop[j,len(temp):,-1] += nrpertraj

                nnohop[j] += int(temp[-1,1])
                nhop[j]   += nrpertraj - int(temp[-1,1])

    nhop   = pop[:,-1,-1]
    nnohop = np.sum(pop[:,-1,:-1], axis=1)

    for j in range(dynl):
        if nnohop[j]+nhop[j] == 0:
            quit("\nNo trajectories evaluated, quitting program")
        else:
            pop[j] /= (nhop[j] + nnohop[j])
            with open("pop%s%i.dat"%(wp,j),"w") as f:
                for i,k in enumerate(pop[j]):
                    tmpt = i * float(config['Dynamics']['timestep'])
                    f.write("%5.2f "%(tmpt))
                    f.write('%5.5f\n'%k)


    print("\n DYNAMICS OVERALL \n------------------")
    print("Of %i nuclear trajectories, %i returned, %i errors"%(
           trajs,ntrajs,nerror))
    for j in range(dynl):
        print("\nOption %i"%j)
        print("   Trajectories that hopped : %7i"%(nhop[j]))
        print("   Trajectories still bound : %7i"%(nnohop[j]))
        print("\n%5.1f %% OF RETURNED TRAJECTORIES HOPPED"%(
               100 * nhop[j] / (nnohop[j]+nhop[j])))
    print("\n")


def freeElectron(trajs, wavepacket):
    """
    combines all TRAJ_X.out/freeEl.dat files into one
    !!! This function was not properly tested !!!
    """

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    for j, jDistr in enumerate(kDistr):
        with open("freeelectrons%i%s.dat"%(j,wp),"w") as f:
            out  = "#   t/fs      "
            out += "kx           "
            out += "ky           "
            out += "kz          "
            out += "Eel/a.u.       "
            out += "VDE/eV\n"
            f.write(out)

    for i in range(trajs):
        if os.path.exists("TRAJ_%i.out"%i):
            print("Trajectory %i"%i)

            for j, jDistr in enumerate(kDistr):
                if os.path.exists("TRAJ_%i.out/freeEl.dat%i%s"%(i,j,wp)):
                    temp = np.loadtxt("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp))
                else:
                    continue

                if np.shape(temp) == (6,):
                    temp = np.asarray([temp])

                with open("freeelectrons%i%s.dat"%(j,wp),"a") as f:
                    for k in temp:
                        f.write("%8.2f %12.5e "%(k[0], k[1]))
                        f.write("%12.5e %12.5e "%(k[2], k[3]))
                        f.write("%14.7e %14.7e\n"%(k[4], k[5]))


def ionDistance(trajs, atomlist, wavepacket):
    """
    calculates atomic distances for given centers at all hopping instances and
    writes it to HoppingDistances/OptionX/
    !!! This function was not properly tested !!!
    """

    atomlist = ast.literal_eval(atomlist)
    if isinstance(atomlist[0], int):
        atomlist = [atomlist]

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    if not os.path.isdir("HoppingDistances%s"%wp):
        os.system("mkdir HoppingDistances%s"%wp)
    if os.listdir("HoppingDistances%s"%wp):
        os.system("rm -r HoppingDistances%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            xyz = f.readlines()
        nat = int(xyz[0])

        for j, jDistr in enumerate(kDistr):
            path = "HoppingDistances%s/Option%i"%(wp,j)
            if not os.path.isdir(path):
                os.system("mkdir %s"%path)

            if os.path.exists("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp)):
                freeEl = np.loadtxt("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp))
            else:
                print("   No hops in trajectory (option %i)"%j)
                continue

            if np.shape(freeEl) == (6,):
                freeEl = freeEl.reshape((1,6))

            steps = np.asarray(freeEl[:,0]*5, dtype=int)

            for k in atomlist:
                with open("%s/%i-%i%s.dat"%(path, k[0], k[1],wp), "a") as f:
                    for l in steps:
                        a1  = np.asarray(xyz[(nat+2)*l+1+k[0]].split()[1:],
                                        dtype=float)
                        a2  = np.asarray(xyz[(nat+2)*l+1+k[1]].split()[1:],
                                        dtype=float)

                        l12 = la.norm(a1 - a2)

                        f.write("%5i %12.9f\n"%(l,l12))


def hoppingAngle(trajs, inp, wavepacket):
    """
    calculates given angles at all hopping instances and writes it to
    HoppingAngles/OptionX/
    !!! This function was not properly tested !!!
    """

    inp = ast.literal_eval(inp)
    if isinstance(inp[0], int):
        inp = [inp]

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    if not os.path.isdir("HoppingAngles%s"%wp):
        os.system("mkdir HoppingAngles%s"%wp)
    if os.listdir("HoppingAngles%s"%wp):
        os.system("rm -r HoppingAngles%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            xyz = f.readlines()
        nat = int(xyz[0])

        for j, jDistr in enumerate(kDistr):
            path = "HoppingAngles%s/Option%i"%(wp,j)
            if not os.path.isdir(path):
                os.system("mkdir %s"%path)

            if os.path.exists("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp)):
                freeEl = np.loadtxt("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp))
            else:
                print("   No hops in trajectory (option %i)"%j)
                continue

            if np.shape(freeEl) == (6,):
                freeEl = freeEl.reshape((1,6))

            steps = np.asarray(freeEl[:,0]*5, dtype=int)

            for k in inp:
                fn = "%i-%i-%i.dat"%(k[0], k[1], k[2])

                with open("%s/%s"%(path,fn),"a") as f:
                    for l in steps:
                        a1  = np.asarray(xyz[(nat+2)*l+1+k[0]].split()[1:],
                                            dtype=float)
                        a2  = np.asarray(xyz[(nat+2)*l+1+k[1]].split()[1:],
                                            dtype=float)
                        a3  = np.asarray(xyz[(nat+2)*l+1+k[2]].split()[1:],
                                            dtype=float)

                        d1  = a1 - a2
                        d2  = a3 - a2
                        angle = 180 / np.pi * (np.arccos(np.dot(d1, d2) /
                                              (la.norm(d1) * la.norm(d2))))
                        if angle > 180:
                            angle = 360 - angle

                        f.write("%5i %6.3f\n"%(l,angle))


def hoppingDihedrals(trajs, inp, wavepacket):
    """
    calculates dihedral angles at all hopping instances and writes it to
    HoppingDihedrals/OptionX/
    !!! This function was not properly tested !!!
    """

    inp = ast.literal_eval(inp)
    if isinstance(inp[0], int):
        inp = [inp]

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    if not os.path.isdir("HoppingDihedrals%s"%wp):
        os.system("mkdir HoppingDihedrals%s"%wp)
    os.system("rm -rf HoppingDihedrals%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            xyz = f.readlines()
        nat = int(xyz[0])

        for j, jDistr in enumerate(kDistr):
            path = "HoppingAngles%s/Option%i"%(wp,j)
            if not os.path.isdir(path):
                os.system("mkdir %s"%path)

            if os.path.exists("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp)):
                freeEl = np.loadtxt("TRAJ_%i.out/freeEl%i%s.dat"%(i,j,wp))
            else:
                print("   No hops in trajectory (option %i)"%j)
                continue

            if np.shape(freeEl) == (6,):
                freeEl = freeEl.reshape((1,6))

            steps = np.asarray(freeEl[:,0]*5, dtype=int)

            for k in inp:
                fn = "%i-%i-%i-%i.dat"%(k[0], k[1], k[2], k[3])

                with open("%s/%s"%(path,fn),"a") as f:
                    for l in steps:
                        a0  = np.asarray(xyz[(nat+2)*l+1+k[0]].split()[1:],
                                            dtype=float)
                        a1  = np.asarray(xyz[(nat+2)*l+1+k[1]].split()[1:],
                                            dtype=float)
                        a2  = np.asarray(xyz[(nat+2)*l+1+k[2]].split()[1:],
                                            dtype=float)
                        a3  = np.asarray(xyz[(nat+2)*l+1+k[3]].split()[1:],
                                            dtype=float)

                        d1, d2 = a0 - a1, a3 - a2
                        dx = a2 - a1

                        d1 -= np.dot(d1,dx)/la.norm(dx)**2*dx
                        d2 -= np.dot(d2,dx)/la.norm(dx)**2*dx

                        det = np.dot(dx, np.cross(d1,d2))
                        testa = np.arctan2(det, np.dot(d1,d2))

                        angle = 180/np.pi * (np.arccos(np.dot(d1, d2) /
                                            (la.norm(d1) * la.norm(d2))))
                        if testa > 0:
                            angle = 360 - angle

                        f.write("%5i %6.3f\n"%(l,angle))

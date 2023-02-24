#!/usr/bin/env python3

import os
import ast
import subprocess as sp
from configparser import ConfigParser

import numpy as np
import numpy.linalg as la

config = ConfigParser()
config.read('config.ini')
kDistr = config['Continuum']['kDistr']

################################################################################
#                                                                              #
#                       pre- and post-dynamics features                        #
#                                                                              #
################################################################################


def inputCheck(indict):
    anyOption = False

    if indict['submit']:
        nproc, tottime, qcmethod = indict['submit']
        writeSub(nproc, tottime, qcmethod)
        print("\nSlurm submit script HORTENSIAarray.sub generated\n")
        anyOption = True

    if indict['initg']:
        f1, f2 = indict['initg']
        wignerStructuresQChem(f1, f2)
        quit("structure.in files created in INITSTRUCT/")

    if indict['check']:
        if indict['ntrajs'] == None:
            quit("\nTotal number of trajectories (option --ntrajs) needed "
                 "as input!\n")
        trajCheck(indict['ntrajs'], indict['wavepacket'])
        anyOption = True

    if indict['electron']:
        if indict['ntrajs'] == None:
            quit("\nTotal number of trajectories (option --ntrajs) needed "
                 "as input!\n")
        freeElectron(indict['ntrajs'], indict['wavepacket'])
        print("\nfreeelectrons.dat written\n")
        anyOption = True

    if indict['dlist']:
        if indict['ntrajs'] == None:
            quit("\nTotal number of trajectories (option --ntrajs) needed "
                 "as input!\n")
        ionDistance(indict['ntrajs'], indict['dlist'], indict['wavepacket'])
        print("\nHopping distances calculated!\n")
        anyOption = True

    if indict['alist']:
        if indict['ntrajs'] == None:
            quit("\nTotal number of trajectories (option --ntrajs) needed "
                 "as input!\n")
        hoppingAngle(indict['ntrajs'], indict['alist'], indict['wavepacket'])
        print("\nHopping angles calculated!\n")
        anyOption = True

    if indict['dhlist']:
        if indict['ntrajs'] == None:
            quit("\nTotal number of trajectories (option --ntrajs) needed "
                 "as input!\n")
        hoppingDihedrals(indict['ntrajs'], indict['dhlist'],
                         indict['wavepacket'])
        print("\nHopping dihedrals calculated!\n")
        anyOption = True

    if indict['clean']:
        cleanFolder()
        quit("\nFolder cleaned from all unnecessary files!\n")

    return anyOption


def writeSub(nproc, tottime, qcmethod):
    """
    Generates an example slurm submit script according to how we use the
    program, quite definitely individual changes have to be made
    """

    output  = "#!/bin/bash\n"
    output += "#SBATCH --nodes=1\n"
    output += "#SBATCH --cpus-per-task=%s\n"%nproc
    output += "#SBATCH --mem-per-cpu=3gb\n"
    output += "#SBATCH -p batch\n"
    output += "#SBATCH --time=60-0\n\n"
    if qcmethod == "g09":
        output += "module load g09-source-local\n\n"
    elif qcmethod == "g16":
        output += "module load g16\n\n"
    else:
        output += "module load qchem\n\n"
    output += "name=TRAJ_$SLURM_ARRAY_TASK_ID\n\n"
    output += "cp -r ${SLURM_SUBMIT_DIR}/${name} "
    output += "/scratch/${SLURM_JOBID}/${name}\n"
    output += "cd /scratch/${SLURM_JOBID}/${name}\n\n"
    if qcmethod != "qchem":
        output += "export GAUSS_SCRDIR=$PWD\n"
    output += "export OMP_NUM_THREADS=%s\n"%nproc
    output += "export MKL_NUM_THREADS=%s\n\n"%nproc
    output += "hortensia --run 1> mist.out 2> error.out\n\n"
    output += "tottime=%s\n\n"%tottime
    output += "nowtime=$( tail -n 1 trajpop.dat | awk '{print$1}' )\n"
    output += "trajpop=$( tail -n 1 trajpop.dat | awk '{print$2}' )\n\n"
    output += "if [ $trajpop == 0 ] || [ $nowtime == $tottime ] \nthen\n    "
    output += "mv /scratch/${SLURM_JOBID}/${name} "
    output += "${SLURM_SUBMIT_DIR}/${name}.out\nelse\n    "
    output += "mv /scratch/${SLURM_JOBID}/${name} "
    output += "${SLURM_SUBMIT_DIR}/${name}.error\nfi"

    with open("HORTENSIAarray.sub","w") as f:
        f.write(output)


def wignerStructures(f1, f2):
    """
    Starts the Wigner program for the generation of initial conditions
    """

    if not os.path.exists("INITSTRUCT"):
        os.system("mkdir INITSTRUCT")
    else:
        os.system("rm -rf INITSTRUCT/*")

    os.system("cp %s %s INITSTRUCT"%(f1,f2))

    cwd = os.getcwd()
    os.chdir(cwd+"/INITSTRUCT")

    linflag = input("Is molecule linear? y/[n]  ")
    if linflag.lower() == "y":
        linflag = 1
    else:
        linflag = 0

    dist_type = input("Choose the desired distribution. "
                      "Wigner = \'W\', Husimi = \'H\', "
                      "|psi(x)|^2*|phi(q)|^2 for v=1: \'V1\'.    ")
    if (dist_type.upper() == "H") or (dist_type.upper() == "V1"):
        # excited modes only for Husimi distribution or V1
        tmpstr   = "Enter the index of normal modes in v=1. Syntax: 1 2 3 4.   "
        excmodes = np.asarray(input(tmpstr).split(), dtype=int)
    else:
        excmodes = []

    if dist_type == "W":
        # temperature effects only for Wigner distribution
        temp = float(input("Desired temperature in K.  "))
    else:
        temp = 0.1

    ncond = int(input("Number of initial conditions.  "))

    from hortensia_latest.wigner.WignerBase import WignerEnsemble
    from hortensia_latest.wigner.HarmonicGaussian import HarmonicGaussian

    DistClass = WignerEnsemble(HarmonicGaussian(linflag, f1, f2), temp, ncond,
                               "structure", dist_type, excmodes)
    DistClass.getWignerEnsemble()


def wignerStructuresQChem(f1, f2):
    """
    Starts the Wigner program for the generation of initial conditions
    """

    if not os.path.exists("INITSTRUCT"):
        os.system("mkdir INITSTRUCT")
    else:
        os.system("rm -rf INITSTRUCT/*")

    os.system("cp %s %s INITSTRUCT"%(f1,f2))

    cwd = os.getcwd()
    os.chdir(cwd+"/INITSTRUCT")

    linflag = input("Is molecule linear? y/[n]  ")
    if linflag.lower() == "y":
        linflag = 1
    else:
        linflag = 0

    dist_type = input("Choose the desired distribution. "
                      "Wigner = \'W\', Husimi = \'H\', "
                      "|psi(x)|^2*|phi(q)|^2 for v=1: \'V1\'.    ")
    if (dist_type.upper() == "H") or (dist_type.upper() == "V1"):
        # excited modes only for Husimi distribution or V1
        tmpstr   = "Enter the index of normal modes in v=1. Syntax: 1 2 3 4.   "
        excmodes = np.asarray(input(tmpstr).split(), dtype=int)
    else:
        excmodes = []

    if dist_type == "W":
        # temperature effects only for Wigner distribution
        temp = float(input("Desired temperature in K.  "))
    else:
        temp = 0.1

    ncond = int(input("Number of initial conditions.  "))

    from hortensia_latest.wigner.WignerBase import WignerEnsemble
    from hortensia_latest.wigner.HarmonicQChem import HarmonicQChem

    DistClass = WignerEnsemble(HarmonicQChem(linflag, f1, f2), temp, ncond,
                               "structure", dist_type, excmodes)
    DistClass.getWignerEnsemble()


def trajCheck(trajs, wavepacket):
    """
    Checks for erroneous trajectories and moves them to .error,
    evaluates the total hopping statistics and writes pop.dat files with the
    total population of every time step
    """

    if isinstance(kDistr, list):
        import hortensia_latest.analysis_list as analysis_list
        analysis_list.trajCheck(trajs, wavepacket)

        return

    out = "out.out"

    states = len(ast.literal_eval(config['States']['anion_states']))
    nrpertraj = int(config['Dynamics']['nrpertraj'])
    ntrajs = 0
    nerror = 0

    pop = np.zeros((int(config['Dynamics']['steps']) + 1, states + 1))

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

            bound = np.loadtxt("TRAJ_%i.out/boundPop.dat"%(i))
            temp  = np.loadtxt("TRAJ_%i.out/trajpop%s.dat"%(i,wp))
            bound = np.asarray(bound[:,1], dtype=int)
            temp  = np.asarray(temp[:,1] , dtype=int)
            print("   Hopped into continuum: %5i"%(nrpertraj-temp[-1]))
            print("   Stayed bound         : %5i"%temp[-1])

            for j,k in enumerate(temp):
                pop[j,bound[j]] += k
                pop[j,-1]       += nrpertraj - k
            pop[len(temp):,-1] += nrpertraj

    nhop   = pop[-1,-1]
    nnohop = np.sum(pop[-1,:-1])

    if nnohop + nhop == 0:
        quit("\nNo trajectories evaluated, quitting program")
    else:
        pop /= (nhop + nnohop)
        with open("pop%s.dat"%wp,"w") as f:
            for i,j in enumerate(pop):
                f.write("%5.2f "%(i * float(config['Dynamics']['timestep'])))
                for k in j:
                    f.write("%5.5f "%k)
                f.write("\n")

    print("\n DYNAMICS OVERALL \n------------------")
    print("Of %i nuclear trajectories, %i returned, %i errors"%(
           trajs, ntrajs, nerror))

    print("   Trajectories that hopped : %7i"%nhop)
    print("   Trajectories still bound : %7i"%nnohop)
    print("\n%5.1f %% OF RETURNED TRAJECTORIES HOPPED\n"%(
           100 * nhop / (nnohop+nhop)))


def freeElectron(trajs, wavepacket):
    """
    combines all freeEl.dat files into one
    """

    if isinstance(kDistr, list):
        import hortensia_latest.analysis_list as analysis_list
        analysis_list.freeElectron(trajs, wavepacket)

        return

    if wavepacket:
        wp = 'WP'
    else:
        wp = ''

    with open("freeelectrons%s.dat"%wp,"w") as f:
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

            if os.path.exists("TRAJ_%i.out/freeEl%s.dat"%(i,wp)):
                temp = np.loadtxt("TRAJ_%i.out/freeEl%s.dat"%(i,wp))
            else:
                continue

            if np.shape(temp) == (6,):
                temp = np.asarray([temp])

            with open("freeelectrons%s.dat"%(wp),"a") as f:
                for k in temp:
                    f.write("%8.2f %12.5e "%(k[0], k[1]))
                    f.write("%12.5e %12.5e "%(k[2], k[3]))
                    f.write("%14.7e %14.7e\n"%(k[4], k[5]))


def ionDistance(trajs, atomlist, wavepacket):
    """
    calculates atomic distances for given centers at all hopping instances and
    writes it to HoppingDistances/
    """

    if isinstance(kDistr, list):
        import hortensia_latest.analysis_list as analysis_list
        analysis_list.ionDistance(trajs, atomlist, wavepacket)

        return

    atomlist = ast.literal_eval(atomlist)
    if isinstance(atomlist[0], int):
        atomlist = [atomlist]

    wp = ('WP' if wavepacket else '')

    os.system("mkdir -p HoppingDistances%s"%wp)
    os.system("rm -rf HoppingDistances%s/*"%wp)
    os.system("mkdir -p AllDistances%s"%wp)
    os.system("rm -rf AllDistances%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        if os.path.exists("TRAJ_%i.out/freeEl%s.dat"%(i,wp)):
            freeEl = np.loadtxt("TRAJ_%i.out/freeEl%s.dat"%(i,wp))
        else:
            print("   No hops in trajectory")
            continue

        if np.shape(freeEl) == (6,):
            freeEl = freeEl.reshape((1,6))

        steps = np.asarray(freeEl[:,0]*5, dtype=int)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            l = f.readlines()
        nat = int(l[0])

        totsteps = len(l)//(nat+2)
        xyz = np.zeros((totsteps, nat, 3), dtype=float)

        l12 = []
        [l12.append([]) for _ in range(len(atomlist))]
        for t in range(totsteps):
            for x in range(nat):
                xyz[t, x] = l[t*(nat+2)+x+2].split()[1:4]
            
            for j, at in enumerate(atomlist):
                l12[j] += [la.norm(xyz[t, at[1]-1] - xyz[t, at[0]-1])]

        for j, at in enumerate(atomlist):
            with open("AllDistances%s/%i-%i.dat"%(wp,at[0],at[1]),"a") as f:
                for t, d in enumerate(l12[j]):
                    f.write("%5i %12.9f\n"%(t, d))

        for j in atomlist:
            with open("HoppingDistances%s/%i-%i.dat"%(wp,j[0],j[1]),"a") as f:
                for k in steps:
                    l12 = la.norm(xyz[k, j[1]-1] - xyz[k, j[0]-1])
                    f.write("%5i %12.9f\n"%(k, l12))


def hoppingAngle(trajs, inp, wavepacket):
    """
    calculates given angles at all hopping instances and writes it to
    HoppingAngles/
    """

    if isinstance(kDistr, list):
        import hortensia_latest.analysis_list as analysis_list
        analysis_list.hoppingAngle(trajs, inp, wavepacket)

        return

    import ast

    inp = ast.literal_eval(inp)
    if isinstance(inp[0], int):
        inp = [inp]

    wp = ('WP' if wavepacket else '')

    os.system("mkdir -p HoppingAngles%s"%wp)
    os.system("rm -rf HoppingAngles%s/*"%wp)
    os.system("mkdir -p AllAngles%s"%wp)
    os.system("rm -rf AllAngles%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        if os.path.exists("TRAJ_%i.out/freeEl%s.dat"%(i,wp)):
            freeEl = np.loadtxt("TRAJ_%i.out/freeEl%s.dat"%(i,wp))
        else:
            print("   No hops in trajectory")
            continue

        if np.shape(freeEl) == (6,):
            freeEl = freeEl.reshape((1,6))

        steps = np.asarray(freeEl[:,0]*5, dtype=int)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            l = f.readlines()
        nat = int(l[0])

        totsteps = len(l)//(nat+2)
        xyz = np.zeros((totsteps, nat, 3), dtype=float)

        angle = []
        [angle.append([]) for _ in range(len(inp))]
        for t in range(totsteps):
            for x in range(nat):
                xyz[t, x] = l[t*(nat+2)+x+2].split()[1:4]
            
            for j, at in enumerate(inp):
                d1 = xyz[t, at[0]-1] - xyz[t, at[1]-1]
                d2 = xyz[t, at[2]-1] - xyz[t, at[1]-1]

                a = 180 / np.pi * (np.arccos(np.dot(d1, d2) /
                                   (la.norm(d1) * la.norm(d2))))
                if a > 180:
                    a = 360 - a

                angle[j] += [a]

        for j, at in enumerate(inp):
            filename = "%i-%i-%i.dat"%(at[0], at[1], at[2])
            with open("AllAngles%s/%s"%(wp, filename), "a") as f:
                for t, a in enumerate(angle[j]):
                    f.write("%5i %6.3f\n"%(t, a))

        for j in inp:
            filename = "%i-%i-%i.dat"%(j[0], j[1], j[2])
            with open("HoppingAngles%s/%s"%(wp, filename),"a") as f:
                for k in steps:
                    d1 = xyz[k, at[0]-1] - xyz[k, at[1]-1]
                    d2 = xyz[k, at[2]-1] - xyz[k, at[1]-1]

                    angle = 180 / np.pi * (np.arccos(np.dot(d1, d2) /
                                          (la.norm(d1) * la.norm(d2))))
                    if angle > 180:
                        angle = 360 - angle

                    f.write("%5i %6.3f\n"%(k,angle))


def hoppingDihedrals(trajs, inp, wavepacket):
    """
    calculates dihedral angles at all hopping instances and writes it to
    HoppingDihedrals/
    """

    if isinstance(kDistr, list):
        import hortensia_latest.analysis_list as analysis_list
        analysis_list.hoppingDihedrals(trajs, inp, wavepacket)

        return

    import ast

    inp = ast.literal_eval(inp)
    if isinstance(inp[0], int):
        inp = [inp]

    wp = ('WP' if wavepacket else '')

    os.system("mkdir -p HoppingDihedrals%s"%wp)
    os.system("rm -rf HoppingDihedrals%s/*"%wp)
    os.system("mkdir -p AllDihedrals%s"%wp)
    os.system("rm -rf AllDihedrals%s/*"%wp)

    for i in range(trajs):
        if not os.path.isdir("TRAJ_%i.out"%i):
            continue

        print("Trajectory %i calculating"%i)

        if os.path.exists("TRAJ_%i.out/freeEl%s.dat"%(i,wp)):
            freeEl = np.loadtxt("TRAJ_%i.out/freeEl%s.dat"%(i,wp))
        else:
            print("   No hops in trajectory")
            continue

        if np.shape(freeEl) == (6,):
            freeEl = freeEl.reshape((1,6))

        steps = np.asarray(freeEl[:,0]*5, dtype=int)

        with open("TRAJ_%i.out/dynamics.xyz"%i,"r") as f:
            l = f.readlines()
        nat = int(l[0])

        totsteps = len(l)//(nat+2)
        xyz = np.zeros((totsteps, nat, 3), dtype=float)

        angle = []
        [angle.append([]) for _ in range(len(inp))]
        for t in range(totsteps):
            for x in range(nat):
                xyz[t, x] = l[t*(nat+2)+x+2].split()[1:4]
            
            for j, at in enumerate(inp):
                d1  = xyz[t, at[0]-1] - xyz[t, at[1]-1]
                d2  = xyz[t, at[3]-1] - xyz[t, at[2]-1]
                dx  = xyz[t, at[2]-1] - xyz[t, at[1]-1]

                d1 -= np.dot(d1, dx) / la.norm(dx)**2 * dx
                d2 -= np.dot(d2, dx) / la.norm(dx)**2 * dx

                det = np.dot(dx, np.cross(d1, d2))
                testa = np.arctan2(det, np.dot(d1, d2))

                a = 180/np.pi * (np.arccos(np.dot(d1, d2) /
                                           (la.norm(d1) * la.norm(d2))))
                if testa > 0:
                    a = 360 - a

                angle[j] += [a]

        for j, at in enumerate(inp):
            filename = "%i-%i-%i-%i.dat"%(at[0], at[1], at[2], at[3])
            with open("AllDihedrals%s/%s"%(wp, filename), "a") as f:
                for t, a in enumerate(angle[j]):
                    f.write("%5i %6.3f\n"%(t, a))

        for j in inp:
            filename = "%i-%i-%i-%i.dat"%(j[0], j[1], j[2], j[3])

            with open("HoppingDihedrals%s/%s"%(wp,filename),"a") as f:
                for k in steps:
                    d1  = xyz[k, at[0]-1] - xyz[k, at[1]-1]
                    d2  = xyz[k, at[3]-1] - xyz[k, at[2]-1]
                    dx  = xyz[k, at[2]-1] - xyz[k, at[1]-1]

                    d1 -= np.dot(d1,dx)/la.norm(dx)**2*dx
                    d2 -= np.dot(d2,dx)/la.norm(dx)**2*dx

                    det = np.dot(dx, np.cross(d1,d2))
                    testa = np.arctan2(det, np.dot(d1,d2))

                    angle = 180/np.pi * (np.arccos(np.dot(d1, d2) /
                                         (la.norm(d1) * la.norm(d2))))
                    if testa > 0:
                        angle = 360 - angle

                    f.write("%5i %6.3f\n"%(k,angle))


def cleanFolder():
    data = "neu.* an.* dynamics.xyz out.out Gau* *.dat geom*"
    os.system("rm -f %s"%data)

#!/usr/bin/env python3

import os
import sys
import subprocess as sp

import numpy as np

import hortensia_latest.misc as misc

sys.path.insert(0, os.getcwd())

def externalQC(stateSet, calcSet, molSet, fileNr, fo):
    """
    Generates input and carries out quantum-chemical calculations with the
    Gaussian program package. Reads and returns forces
    """

    states , state   = stateSet[0], stateSet[1]
    statenr, excited = stateSet[2], stateSet[3]

    Nproc  , func    = calcSet[0] , calcSet[1]
    cartBFs, conv    = calcSet[2] , calcSet[3]
    maxcyc , method  = calcSet[4] , calcSet[5].lower()

    coord  , atS     = molSet[0], molSet[1]
    charge , mult    = molSet[2], molSet[3]

    def makeInput(rerun=False):
        """
        Generates input for Gaussian
        """
        # Head
        inp     = "%Nproc="+"%i\n"%Nproc
        inp    += "%Mem="+"%iMb\n"%(2500*Nproc)
        inp    += "#p %s/GEN "%func
        inp    += len(cartBFs)*"%s"%tuple([i for i in cartBFs])
        inp    += "Iop(6/7=3) Iop(9/40=14) "
        if rerun:
            inp += "scf(maxcyc=%i,qc,conver=%i) nosym "%(maxcyc, conv)
        else:
            inp += "scf(maxcyc=%i,conver=%i) "%(maxcyc, conv)
            inp += "nosym Guess(Read,TCheck) "
        inp    += "GFInput "
        aninp   = "%Chk=an.chk\n"  + inp
        neuinp  = "%Chk=neu.chk\n" + inp

        if statenr == 0:
            aninp += "Force"

        aninp += "\n\nComment\n\n%i %i\n"%(charge, mult)

        if mult == 1:
            neuinp += "\n\nComment\n\n%i %i\n"%(charge-1, mult+1)
        else:
            neuinp += "\n\nComment\n\n%i %i\n"%(charge+1, mult-1)

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

        # Tail 2, if excited states are involved
        if excited:
            tail2  = "--Link1--\n%Chk=an.chk\n"
            tail2 += "%Nproc="+"%i\n"%Nproc
            tail2 += "%Mem="+"%iMb\n"%(2500*Nproc)
            tail2 += "#p TD(full,nstates=%i"%max(states)
            if state != 0:
                tail2 += ",root=%i) Force "%state
            else:
                tail2 += ") "
            tail2 += "%s/GEN "%func
            tail2 += len(cartBFs)*"%s"%tuple([i for i in cartBFs])
            tail2 += "Iop(6/7=3) Iop(9/40=14) "
            tail2 += "scf(maxcyc=%i,conver=%i) "%(maxcyc, conv)
            tail2 += "nosym Guess(Read) Geom=AllCheck\n\n"
            tail2 += tail
        else:
            tail2  = ""

        return aninp+struct+tail+tail2, neuinp+struct+tail


    aninp, neuinp = makeInput()

    with open("an.in"+fileNr,"w") as f:
        f.write(aninp)
    with open("neu.in"+fileNr,"w") as f:
        f.write(neuinp)

    # setting enviromental variable to ensure we calculate in '.'
    cwd = os.getcwd()
    os.environ['GAUSS_SCRDIR'] = cwd

    # !!!
    if os.path.exists("an.log2"):
        os.system("mv an.log2 an.log1")
    if os.path.exists("neu.log2"):
        os.system("mv neu.log2 neu.log1")

    # actual qc calculation
    fo.write("Calculating anionic system...\n")
    sp.call('%s an.in%s an.log%s'%(method, fileNr, fileNr),
            shell=True)
    fo.write("Calculating ionized system...\n")
    sp.call('%s < neu.in%s > neu.log%s'%(method, fileNr, fileNr),
            shell=True)

    # check for errors in calculation, namely scf convergence and
    # (if excited states are needed) negative excitation energies
    anfail, neufail = False,False
    cmd     = "grep 'Normal termination' an.log%s"%fileNr
    grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    out,err = grepout.communicate()
    if out == b"":
        anfail = True
        fo.write("anion calculation not converged, " +
                    "restarting calculation\n")
    if excited:
        with open("an.log%s"%fileNr,"r") as f:
            for i in f:
                if "Excited State" in i:
                    exEn1 = float(i.split()[4])
                    if exEn1 < 0:
                        anfail = True
                        fo.write("Negative excitation energy in anion " +
                                "calculation, restarting calculation\n")
                        break

    cmd     = "grep 'Normal termination' neu.log%s"%fileNr
    grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
    out,err = grepout.communicate()
    if out == b"":
        neufail = True
        fo.write("neutral calculation not converged, " +
                    "restarting calculation\n")

    # if qc calculation failed, calculation is restarted once with the
    # qc scf option and without .chk file of previous calculation
    if anfail:
        fo.write("an.chk removed, switched to QC SCF procedure\n")
        fo.write("Retrying calculation of anion...\n")
        os.system("rm an.chk")

        aninp, neuinp = makeInput(True)
        with open("an.in"+fileNr,"w") as f:
            f.write(aninp)

        sp.call('%s an.in%s an.log%s'%(method, fileNr, fileNr),
                shell=True)

        with open("an.log%s"%fileNr,"r") as f:
            for i in f:
                if "Convergence criterion not met" in i:
                    fo.write("anion calculation still not converged, " +
                                "aborting program\n")
                    quit("SCF for anion not converged!")

                if excited:
                    if "Excited State" in i:
                        exEn1 = float(i.split()[4])
                        if exEn1 < 0:
                            fo.write("Still negative excitation energy in "+
                                        "anion calculation, aborting "+
                                        "program\n")
                            quit("Negative excitation energy in anion!")
    if neufail:
        fo.write("neu.chk removed, switched to QC SCF procedure\n")
        fo.write("Retrying calculation of neutral molecule...\n")
        os.system("rm neu.chk")

        aninp, neuinp = makeInput(True)
        with open("neu.in"+fileNr,"w") as f:
            f.write(neuinp)

        sp.call('%s < neu.in%s > neu.log%s'%(method, fileNr, fileNr),
                shell=True)

        with open("neu.log%s"%fileNr,"r") as f:
            for i in f:
                if "Convergence criterion not met" in i:
                    fo.write("neutral calculation still not converged, " +
                                "aborting program\n")
                    quit("SCF for neutral not converged!")

    acc = readAcc(atS, "an.log%s"%fileNr)

    return acc


def readAcc(atS, log):
    """
    extracts atomic forces from Gaussian output file
    """

    nat = len(atS)
    acc = []

    with open(log,"r") as f:
        i = f.readline()
        while i:
            if "Forces (" in i:
                f.readline()
                f.readline()

                for j in range(nat):
                    acc.append(f.readline().split()[2:])
                break

            i = f.readline()

    acc = np.asarray(acc, dtype=float)

    for i, j in enumerate(acc):
        acc[i] /= misc.masses(atS[i])

    return acc


def stateEnergy(statenr, maxState):
    if statenr == 0:
        cmd     = "grep 'SCF D' an.log2 |awk '{print$5}'"
        grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        out,err = grepout.communicate()
        stdout  = float(out.split()[0])
    else:
        cmd     = "grep 'SCF D' an.log2 |awk '{print$5}'"
        grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        out,err = grepout.communicate()
        gsEn    = float(out.split()[0])
        cmd     = "grep 'Excited S' an.log2 |awk '{print$5}'"
        grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        out,err = grepout.communicate()
        excEn   = float(out.split()[statenr-1])/misc.hartree_to_eV
        stdout  = gsEn + excEn

    return stdout

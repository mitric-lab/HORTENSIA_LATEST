#!/usr/bin/env python3

import os
import subprocess as sp
import numpy as np

import hortensia_latest.misc as misc

def externalQC(stateSet, calcSet, molSet, fileNr, fo):
    """
    Generates input and carries out quantum-chemical calculations with the
    QChem program package. Reads and returns forces
    """

    states , state   = stateSet[0], stateSet[1]
    statenr, excited = stateSet[2], stateSet[3]

    Nproc  , func    = calcSet[0] , calcSet[1]
    cartBFs, conv    = calcSet[2] , calcSet[3]
    maxcyc           = calcSet[4]

    coord  , atS     = molSet[0], molSet[1]
    charge , mult    = molSet[2], molSet[3]

    if fileNr == "1":
        if os.path.isdir("tmpA"):
            os.system("rm -r tmpA/*")
        else:
            os.system("mkdir tmpA")
        if os.path.isdir("tmpN"):
            os.system("rm -r tmpN/*")
        else:
            os.system("mkdir tmpN")

        inp0A  = "$molecule\n  READ geomA\n$end\n\n$rem\n"
        inp0N  = "$molecule\n  READ geomN\n$end\n\n$rem\n"
        inp1   = "MEM_TOTAL        %i\n"%(2700*Nproc)
        inp1  += "MEM_STATIC       2000\n"
        inp1  += "AO2MO_DISK       50000\n"
        inp1  += "GUI              2\n\n"
        inp1  += "SCF_ALGORITHM    DIIS\n"
        inp1  += "SCF_CONVERGENCE  %i\n"%(conv)
        inp1  += "SCF_MAX_CYCLES   %i\n"%(maxcyc)
        inp1  += "SYMMETRY_IGNORE  TRUE\n\n"
        inp1  += "METHOD           %s\n"%func
        inp1  += "BASIS            GEN\n"
        inp1  += "PURECART         22222\n"
        inp1  += "UNRESTRICTED     TRUE\n"
        anInp  = "JOBTYPE          FORCE\n"
        if excited:
            anInp += "RPA              TRUE\n"
            anInp += "CIS_N_ROOTS      %i\n"%max(states)
            anInp += "CIS_STATE_DERIV  %i\n"%statenr

        neuInp = "JOBTYPE          SP\n"

        inp2   = "$end\n\n$basis\n"

        with open("basis","r") as f:
            basis = f.readlines()

        for i in basis:
            inp2 += i

        inp2  += "$end\n "
        guess  = "SCF_GUESS        READ\n"

        with open("an.in1","w") as f:
            f.write(inp0A + inp1 + anInp + inp2)
        with open("an.in2","w") as f:
            f.write(inp0A + inp1 + guess + anInp + inp2)
        with open("neu.in1","w") as f:
            f.write(inp0N + inp1 + neuInp + inp2)
        with open("neu.in2","w") as f:
            f.write(inp0N + inp1 + guess + neuInp + inp2)

    # Structure
    structA  = "$molecule\n%i %i\n"%(charge, mult)
    if mult == 1:
        structN  = "$molecule\n%i %i\n"%(charge+1, mult+1)
    else:
        structN  = "$molecule\n%i %i\n"%(charge+1, mult-1)
    coord_angs = coord * misc.bohr_to_angs
    for i,j in enumerate(coord_angs):
        structA += "%s %.9f %.9f %.9f\n"%(atS[i], j[0], j[1], j[2])
        structN += "%s %.9f %.9f %.9f\n"%(atS[i], j[0], j[1], j[2])
    structA += "$end"
    structN += "$end"

    with open("geomA","w") as f:
        f.write(structA)
    with open("geomN","w") as f:
        f.write(structN)

    if os.path.exists("an.log2"):
        os.system("mv an.log2 an.log1")
        os.system("mv an.fchk2 an.fchk1")
        os.system("mv neu.log2 neu.log1")
        os.system("mv neu.fchk2 neu.fchk1")

    cwd = os.getcwd()
    os.environ['QCSCRATCH']  = "%s"%cwd

    fo.write("Calculating anionic system...\n")
    cmdAn   = "qchem -nt %i -save "%Nproc
    cmdAn  += "an.in%s an.log%s tmpsaveA > mist.dat"%(fileNr, fileNr)
    os.environ['QCLOCALSCR'] = '%s/tmpA'%cwd
    sp.call(cmdAn, shell=True)

    fo.write("Calculating ionized system...\n")
    cmdNeu  = "qchem -nt %i -save "%Nproc
    cmdNeu += "neu.in%s neu.log%s tmpsaveN > mist.dat"%(fileNr, fileNr)
    os.environ['QCLOCALSCR'] = '%s/tmpN'%cwd
    sp.call(cmdNeu, shell=True)

    os.system("mv an.fchk an.fchk%s"%fileNr)
    os.system("mv neu.fchk neu.fchk%s"%fileNr)

    acc = readAcc(atS)

    return acc


def readAcc(atS):
    nat = len(atS)
    acc = []
    with open("tmpsaveA/GRAD","r") as f:
        lines = f.readlines()
        for i in range(nat):
            acc.append(lines[3+i].split())

    acc = -np.asarray(acc, dtype=float)

    for i, j in enumerate(acc):
        acc[i] /= misc.masses(atS[i])

    return acc


def stateEnergy(statenr, maxState):
    if statenr == 0:
        cmd     = "grep 'Total energy in' an.log2 |awk '{print$9}'"
        grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        out,err = grepout.communicate()
        stdout  = float(out)
    else:
        cmd     = "grep 'Total energy for' an.log2 |awk '{print$6}'"
        grepout = sp.Popen([cmd],shell=True,stdout=sp.PIPE,stderr=sp.STDOUT)
        out,err = grepout.communicate()
        out     = out.split()[-maxState:]
        stdout  = float(out[statenr-1])

    return stdout

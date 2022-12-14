#!/usr/bin/env python3

import numpy as np

from hortensia_latest.wigner.WignerBase import WignerEnsemble

if __name__ == "__main__":
    try:
        interface = int(input("Wigner distribution from normal modes. "+
                              "Quantum chemistry program used? "+
                              "1: Molpro, 2: Gaussian, "+
                              "3: Turbomole, 4: Mndo.  "))
    except ValueError:
        quit("Integer required.")

    linflag = input("For linear molecule enter 1.  ")
    if linflag == "":
        linflag = 0
    else:
        try:
            linflag=int(linflag)
        except ValueError:
            quit("Integer required.")

    intext1  = "Input filename 1 (Standard output of QC program; "
    intext1 += "for Gaussian: formatted checkpoint file).  "
    intext2  = "Input filename 2. Defaults: Molpro: \'\', Gaussian: "
    intext2 += "\'freq.log\', Turbomole: \'coord\', Mndo: \'molden.dat\'.  "
    infile1  = input(intext1)
    infile2  = input(intext2)

    intext1   = "Choose the desired distribution. Wigner = \'W\', "
    intext1  += "Husimi = \'H\', |psi(x)|^2*|phi(q)|^2 for v=1: \'V1\'.    "
    dist_type = input(intext1)
    if (dist_type == "H") or (dist_type == "V1"):
        # excited modes only for Husimi distribution or V1
        tmpstr   = "Enter the index of normal modes in v=1. Syntax: 1 2 3 4.   "
        excmodes = np.asarray(input(tmpstr).split(), dtype=int)
    else:
        excmodes = []

    if dist_type == "W":
        # temperature effects only for Wigner distribution
        try:
            temp = float(input("Desired temperature in K.  "))
        except ValueError:
            quit("Floating point number required.  ")
    else:
        temp = 0.1

    try:
        ncond = int(input("Number of initial conditions.  "))
    except ValueError:
        quit("Integer required.")

    outfile = "structure"

    if interface == 1:
        from HarmonicMolpro import HarmonicMolpro as HarmonicInterface
    elif interface == 2:
        from HarmonicGaussian import HarmonicGaussian as HarmonicInterface
        if infile2 == "":
            infile2 = "freq.log"
    elif interface == 3:
        from HarmonicTurbomole import HarmonicTurbomole as HarmonicInterface
        if infile2 == "":
            infile2 = "coord"
    elif interface == 4:
        from HarmonicMndo import HarmonicMndo as HarmonicInterface
        if infile2 == "":
            infile2 = "molden.dat"

    # Initialisation of WignerEnsemble class
    #
    # 1) HarmonicInterface according to the quantum chemistry program used.
    #    Interface is initialised with flag for linear/nonlinear molecule (0/1)
    #    and by providing names of input files
    #    - MOLPRO: Standard output of frequency calculation
    #    - GAUSSIAN: (i) formatted checkpoint file + standard outputfile
    #                (ii) [AT THE MOMENT THIS DOES NOT WORK]
    #                     only standard outputfile 
    #                     (second filename should be "")
    #    - TURBOMOLE: aoforce outputfile, "coord" file
    #    - MNDO: Standard output of frequency calculation, "molden.dat" file
    # 2) Temperature in K
    # 3) Number of initial conditions to be generated
    # 4) Name string of outputfiles (e.g. "dynamics" or "structure")
    # 5) Type of distribution: Wigner or Husimi
    # 6) Index of vibrationally excited modes

    test = WignerEnsemble(HarmonicInterface(linflag, infile1, infile2),
                          temp, ncond, outfile, dist_type, excmodes)
    test.getWignerEnsemble()

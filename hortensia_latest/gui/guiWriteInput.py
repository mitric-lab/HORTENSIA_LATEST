#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hortensia_latest.info as info

def write(data):
    output  = "### Input configuration file for the "
    output += "HORTENSIA program package\n\n"

    output += "[Dynamics]\n\n"
    output += "### time step for nuclear dynamics in fs\n"
    output += "timestep       = %s\n"%data[16]
    output += "### number of nuclear dynamics steps\n"
    output += "steps          = %s\n"%data[17]
    output += "### number of population dynamics steps per "
    output += "nuclear dynamics step\n"
    output += "popsteps       = %s\n"%data[18]
    output += "### number of times hopping is checked per step "
    output += "and trajectory\n"
    output += "nrpertraj      = %s\n"%data[19]
    output += "### option for restarting "
    output += "(needs dynamics.xyz, velocities.dat, trajpop.dat)\n"
    output += "restart        = %s\n"%data[20]
    output += "### controls the amount of output\n"
    output += "### -1: output every skipstep'th time step, no norm.dat file\n"
    output += "###  0: standard output\n"
    output += "###  1: more detailed printed output\n"
    output += "printlevel     = %s\n"%data[21]
    if int(data[21]) < 0:
        output += "### if printlevel < 0, "
        output += "skipstep = n only prints every nth step\n"
        output += "skipstep       = %s\n"%data[22]
    output += "### type of output for the different output files\n"
    output += "### possible options: - full   : "
    output += "prints values for every free electron state\n"
    output += "###                   - summed : "
    output += "prints sum of all free electron states\n"
    output += "###                   - none   : "
    output += "no output is printed/ energies: an + neu\n"
    output += "coupprint      = %s\n"%data[23][0]
    output += "coefprint      = %s\n"%data[23][1]
    output += "nadiaprint     = %s\n"%data[23][2]
    output += "probprint      = %s\n"%data[23][3]
    output += "enerprint      = %s\n\n\n"%data[23][4]

    output += "[QC]\n\n"
    output += "### Quantum chemistry program "
    output += "(possible options: g09, g16, qchem)\n"
    output += "qcmethod       = %s\n"%data[24]
    output += "### Molecular charge of anion\n"
    output += "charge         = %s\n"%data[25]
    output += "### Multiplicity of anion\n"
    output += "mult           = %s\n"%data[26]
    output += "### Functional used for calculation\n"
    output += "func           = %s\n"%data[27]
    output += "### Amount of threads the Gaussian calculation should run on\n"
    output += "Nproc          = %s\n"%data[28]
    output += "### SCF convergence criterion (10^-N)\n"
    output += "convergence    = %s\n"%data[29]
    output += "### Maximum number of SCF iterations\n"
    output += "maxiter        = %s\n\n\n"%data[30]

    output += "[States]\n\n"
    output += "### list with considered bound states for dynamics\n"
    output += "anion_states   = [%s]\n"%data[9]
    output += "### state in which the dynamics starts, "
    output += "index of anion_states list\n"
    output += "starting_state = %s\n"%data[10]
    output += "### precision with which states/DOs are constructed "
    output += "from micro states (in %)\n"
    output += "precision      = %s\n"%data[11]
    output += "### maximum number of micro states considered for "
    output += "Dyson orbitals\n"
    output += "max_micro      = %s\n"%data[12]
    output += "## energy shift of the ionized ground state (in eV)\n"
    output += "Eshift         = %s\n"%data[13]
    output += "### lifetime of resonance in fs, where exp(-t/(2*tau))\n"
    output += "tau            = %s\n"%data[14]
    output += "### type of orthogonalization (options: 'state', 'mo')\n"
    output += "orthotype      = %s\n\n\n"%data[15]

    output += "[Continuum]\n\n"
    output += "### distribution of k vectors, possible options:\n"
    output += "###   -   fib : Fibonacci spheres for arbitrary "
    output += "numbers of nkfib per energy\n"
    output += "###   -  snub : Snub cube with 24 equally distributed "
    output += "k per energy\n"
    output += "###   - cubic : cuboids with specific nk and Ek in "
    output += "each dimension\n"
    output += "kDistr         = %s\n"%(data[0])
    if data[1] != "":
        output += "### maximum plane wave energy (in eV)\n"
        output += "maxEk          = %s\n"%data[1]
        output += "### number of different plane wave energies\n"
        output += "nEk            = %s\n"%data[2]
    if data[3] != "":
        output += "### number of k per energy\n"
        output += "nkfib          = %s\n"%data[3]
    if data[4] != "":
        output += "### maximum k-value in each direction "
        output += "(in a.u., corresponding to Emax = |k|^2/2)\n"
        output += "maxkx          = %s\n"%data[4][0]
        output += "maxky          = %s\n"%data[4][1]
        output += "maxkz          = %s\n"%data[4][2]
        output += "### number of k-values in each direction\n"
        output += "nkx            = %s\n"%data[5][0]
        output += "nky            = %s\n"%data[5][1]
        output += "nkz            = %s\n"%data[5][2]
    output += "\n\n"

    output += "[2e-Integrals]\n\n"
    output += "### if <PW(k)i|jl> integrals should be evaluated exactly "
    output += "(False: constant PW approximation)\n"
    output += "exact2elInt    = %s\n"%data[6]
    if data[6] == "True":
        output += "### number of k vectors per k energy "
        output += "(possible options: 2/6/8/12)\n"
        output += "intNk          = %s\n"%data[7]
        output += "### calculate integrals for every nth k energy "
        output += "(plus lowest/highest energy)\n"
        output += "intkSkip       = %s\n\n\n"%data[8]
    else:
        output += "\n\n"

    # Tutorial stuff
    output += "#####################################################\n"
    output += "#####  calculation of auto-ionization dynamics  #####\n"
    output += "#####################################################\n\n"
    output += "# 1. Files needed to start the calculations are\n"
    output += "#    - config.ini\n"
    output += "#    - structure.in\n"
    output += "#    - basis\n"
    output += "# 2. If you wish to restart a trajectory,\n"
    output += "#    additional mandatory files are\n"
    output += "#    - dynamics.xyz\n"
    output += "#    - velocities.dat\n"
    output += "#    - tempcoef.dat\n"
    output += "#    - boundPop.dat\n"
    output += "#    - trajpop.dat\n#\n"
    output += "# The program also has some useful functionalities\n"
    output += "# beyond the dynamics simulation. An overview is\n"
    output += "# available through the help function called with\n"
    output += "#    'hortensia --help'\n"

    with open("config.ini","w") as f:
        f.write(output)

def writeSub(nproc, tottime, qcmethod):
    output  = "#!/bin/bash\n"
    output += "#SBATCH --nodes=1\n"
    output += "#SBATCH --cpus-per-task=%i\n"%nproc
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
    output += "export OMP_NUM_THREADS=%i\n"%nproc
    output += "export MKL_NUM_THREADS=%i\n\n"%nproc
    output += "hortensia --run 1> mist.out 2> error.out\n\n"
    output += "tottime=%.2f\n\n"%(tottime)
    output += "nowtime=$( tail -n 1 trajpop.dat | awk '{print$1}' )\n"
    output += "trajpop=$( tail -n 1 trajpop.dat | awk '{print$2}' )\n\n"
    output += "if [ $trajpop == 0 ] || [ $nowtime == $tottime ] \nthen\n    "
    output += "mv /scratch/${SLURM_JOBID}/${name} "
    output += "${SLURM_SUBMIT_DIR}/${name}.out\nelse\n    "
    output += "mv /scratch/${SLURM_JOBID}/${name} "
    output += "${SLURM_SUBMIT_DIR}/${name}.error\nfi"

    with open("HORTENSIAarray.sub","w") as f:
        f.write(output)

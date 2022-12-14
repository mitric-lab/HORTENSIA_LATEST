#!/usr/bin/env python3

#########################################
###  settings for free electron grid  ###
#########################################

### distribution of k vectors, possible options:
###   -   'fib' : Fibonacci spheres for arbitrary numbers of nkfib per energy
###   -  'snub' : Snub cube with 24 equally distributed k per energy
###   - 'cubic' : cuboids with specific nk and Ek in each dimension
kDistr         = ['snub', 'fib', 'cubic', 'snub', 'fib', 'cubic']
kDistr         = 'fib'

##### if kDistr != 'cubic' -> spherical distributions  #####
### maximum plane wave energy (in eV)
maxEk          = [1.5, 1.5, 1.5, 1.5]
maxEk          = 1.5

### number of different plane wave energies
nEk            = [100, 100, 1000, 400]
nEk            = 50

### only if kDistr == 'fib'; number of k per energy
nkfib          = [100, 200]
nkfib          = 1000

##### if kDistr == 'cubic' -> cuboidal distribution  #####
### maximum k-value in each direction (in a.u., corresponding to Emax = |k|^2/2)
maxkx          = [0.2, 0.2]
maxky          = [0.2, 0.2]
maxkz          = [0.2, 0.2]
#maxkx          = 0.2
#maxky          = 0.2
#maxkz          = 0.2

### number of k-values in each direction
nkx            = [8, 9]
nky            = [8, 9]
nkz            = [8, 9]
#nkx            = 30
#nky            = 30
#nkz            = 30


######################################
###  2-electron integral settings  ###
######################################

### if <PW(k)i|jl> integrals should be evaluated exactly (False: constant PW approximation)
exact2elInt    = False
### number of k vectors per k energy (possible options: 2/6/8/12)
intNk          = 12
### calculate integrals for every nth k energy (plus lowest/highest energy)
intkSkip       = 1


############################
###  setting for states  ###
############################

### list with considered bound states for dynamics
anion_states   = [0]
### state in which the dynamics starts, index of anion_states list
starting_state = 0
### precision with which states/DOs are constructed from micro states (in %)
precision      = 99.5
### maximum number of micro states considered for Dyson orbitals
max_micro      = 100
### energy shift of the ionized ground state (in eV)
#Eshift         = -0.0698
Eshift         = 4.0
### lifetime of resonance in fs, where exp(-t/(2*tau))
tau            = 1.0
### type of orthogonalization (options: "state", "mo")
orthotype      = "state"


###########################
###  dynamics settings  ###
###########################

### time step for nuclear dynamics in fs
timestep       = 0.2
### number of nuclear dynamics steps
steps          = 5
### number of population dynamics steps per nuclear dynamics step
popsteps       = 100
### number of times hopping is checked per step and trajectory
nrpertraj      = 100
### option for restarting (needs dynamics.xyz, velocities.dat, trajpop.dat)
restart        = False
### controls the amount of output
### -1: output every skipstep'th time step, no norm.dat file
###  0: standard output
###  1: more detailed printed output
printlevel     = 0
### if printlevel < 0, skipstep = n only prints every nth step
skipstep       = 5
### type of output for the different output files
### possible options: - 'full'   : prints values for every free electron state
###                   - 'summed' : prints sum of all free electron states
###                   - 'none'   : no output is printed/ energies: an + neu
coupprint      = 'summed'
coefprint      = 'full'
nadiaprint     = 'none'
probprint      = 'summed'
enerprint      = 'none'


####################################
###  quantum chemistry settings  ###
####################################

### Quantum chemistry program (possible options: 'g09', 'qchem')
qcmethod       = 'g09'
### Molecular charge of anion
charge         = -1
### Multiplicity of anion
mult           = 2
### Functional used for calculation
func           = "PBEPBE"
### Amount of threads the Gaussian calculation should run on
Nproc          = 2
### SCF convergence criterion (10^-N)
convergence    = 8


#################################################
###  calculation of auto-ionization dynamics  ###
#################################################

# 1. Files needed to start the calculations are
#    - AID.py
#    - configAID.py
#    - structure.in
#    - basis
# 2. Your basis set needs to be decontracted (if this is not the case,
#    './AID.py basis' can decontract your existing basis file)
# 3. For the g09 calculations, the module chem/g09-source from Alexander's
#    module files has to be loaded
# 5. If you wish to restart a trajectory, additional mandatory files are:
#    - dynamics.xyz
#    - velocities.dat
#    - tempcoef.dat
#    - trajpop.dat
#
#
# The program also has some useful functions beyond the dynamics simulation.
# They can be used with the command './AID.py [option] [further options]
# To date these options are:
#
# 'clean'    : clean the current working directory from all output generated by
#              the program (be careful, since ALL .dat files are deleted)
#
# 'basis'    : decontract the basis set given in basis and write it in the
#              basis file (the old basis set is moved to basis.old)
#
# 'init'     : create starting conditions in 'INITSTRUCT/' by calling
#              Jens' WignerEnsemble.py
#
# 'distance' : calculate the atomic distances of selected atoms at hopping steps
#              further options: [number of trajs] [atomlist] [output file name]
#              format of atomlist in a list format with atom numbers as follows:
#              [[1,2],[3,4],...]
#              this will produce output files in the folder 'HoppingDistances'
#
# 'angle'    : calculate selected molecular angles at the hopping steps
#              further options: [number of trajs] [anglelist] [output file name]
#              anglelist is of analogous format as atomlist, so:
#              [[1,2,3],[2,3,4],...]
#              this will produce output files in the folder 'HoppingAngles'
#
# 'electron' : evaluate the k vectors and energy of the free electrons the
#              trajectories hopped in
#              further options: [number of trajs] [output file name]
#
# 'check'    : checks the output trajectories, prints errors and evaluates the
#              overall population
#              further options: [number of trajs] [output file name]

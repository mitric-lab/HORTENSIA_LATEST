#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Frame 1
q0 = "Steps for integration of electronic coefficients per nuclear time step"
q1 = "Initial number of attempted hops per time step and nuclear trajectory."+\
     "\nAfter each successful hop, this number is internally reduced by one."
q2 = "full\t: Couplings between populated and every other state are printed\n"+\
     "summed\t: Couplings between populated and other bound states "+\
     "as well as sum\n\t  of absolute between populated and free-electron "+\
     "states are printed\n"+\
     "none\t: No couplings are printed"
q4 = "full\t: Absolute of electronic coefficients of all states are printed\n"+\
     "summed\t: Absolute of electronic coefficients of bound states and "+\
     "absolute\n\t  of sum of free-electron states are printed\n"+\
     "none\t: No coefficients are printed"
q5 = "full\t: Hopping probabilities of every state are printed\n"+\
     "summed\t: Hopping probabilities of bound states and sum of hopping "+\
     "probabilities\n\t  of free-electron states are printed\n"+\
     "none\t: No coefficients are printed"
q6 = "full\t: All potential state energies are printed\n"+\
     "none\t: All bound state energies and ionized state energy are printed"
q7 = "-1: output written every [Skipped steps]'th time step, no norm.dat\n"+\
     " 0: standard output\n 1: more detailed output for dyson orbitals"
# Frame 3
q8 = "0 is the anion ground state\nIf any state >0, excited state "+\
     "calculations are performed"
q9 = "Precision threshold for the construction of electronic states from "+\
     "single\nSlater determinants\nIf precision threshold is not reached "+\
     "within the maximum number of microstates,\n"+\
     "precision threshold is ignored!"
q10 = "Maximum number of Slater determinants considered for the "+\
      "construction of\nelectronic states\nIf precision threshold is not "+\
      "reached within the maximum number of microstates,\n"+\
      "precision threshold is ignored!"
q11 = "Shifts the ionized system up by this amount, e.g. to match ionization "+\
      "energies\n\nIt is generally advised against this option!"
q12 = "In case of negative vertical detachment energies (VDE), "+\
      "electronic population\ndecays exponentially with this time constant "+\
      "to simulate electronic resonances"
# Frame 4
q13 = "Integrals of the form <Psi(k)i|jl> are evaluated exactly for\n"+\
      "chosen k and then interpolated for the whole k-grid\n"+\
      "As of now, it is advised against this option, since it is "+\
      "highly experimental and\nextremely time-consuming and the "+\
      "interpolation is probably quite inaccurate"
q24 = "Type of orthogonalization between anion and neutral/free electron " + \
      "states\nalways of the form |J_new> = N * (|J> - sum_L |L> <L|J>)"
q25 = "The plane wave functions are orthogonalized with " + \
      "respect to the \noccupied anion MOs\n"
q26 = "The free electrons state are orthogonalized with " + \
      "respect to all anion states\n"
q27 = "The plane wave functions are orthogonalized with respect to the " + \
      "Dyson orbitals"
q28 = "Plane waves are not orthogonalized and instead it is assumed that " + \
      "the \noverlap between anion and free electron states is zero"
# Frame 5
q14 = "It is possible to parallelly compute dynamics with different " + \
      "k-grid distributions.\nThe underlying nuclear dynamics is the " + \
      "same, but hopping is calculated\nindividually"
q15 = "Fibonacci: Distributes k-vectors approximately equally on a sphere " + \
      "according to\n\t  the Fibonacci sphere algorithm\n" + \
      "Snub\t: 24 k-vectors are distributed equally on a sphere to " + \
      "form a Snub cube\n" + \
      "Cubic\t: k-vectors are distributed with equal spacing dx, dy, dz"
q16 = "Maximum kinetic energy of free electron,\n" + \
      "k-vectors are energetically equally spaced"
q29 = "If included, only hops are allowed which contain a minimum kinetic " + \
      "energy equal \nto the zero-point energy. \nAlso if included, the " + \
      "maximum plane wave energy is then given as the \nsum of the " + \
      "zero-point energy difference and the vibrational excitation energy"
q30 = "Zero-point energy difference between the initial and the neutral " + \
      "ground state"
q31 = "Excitation energy of mode excited in the generation of initial " + \
      "conditions"
# Frame 6
q17 = "Number of nuclear trajectories, on which the continuum hopping" + \
      " 'trajectories'\nare based"
q18 = "Generates folders named 'TRAJ_' in the current directory"
q19 = "Copies the files 'basis', 'config.ini' and initial conditions " + \
      "from\n'INITSTRUCT' to the folders 'TRAJ_'"
q20 = "Writes a bash script compatible with Slurm"
# Frame 7
q21 = "Needs the two files specified below in the current working " + \
      "directory.\nCreates the folder 'INITSTRUCT', copies the two files " + \
      "to 'INITSTRUCT' and\ncalls Jens' Wigner program for the generation " + \
      "of initial conditions for\nnuclear trajectories"
q22 = "|\u03a8(x)|\u00B2|\u03a8(q)|\u00B2\t: distribution according to " + \
      "the product of the square of the\n\t\t  momentum and position " + \
      "distribution\n" + \
      "Husimi\t\t: Husimi distribution\n" + \
      "Wigner\t\t: Wigner distribution"
q23 = "Vibrational modes sorted by rising energy, starting with index 1\n\n" + \
      "Temperature in K"

qb = "Generates input files and exits the program.\n\n"+\
     "A prompt will open that lets you generate input according to the "+\
     "provided\nsettings."

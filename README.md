# HORTENSIA

Program for the surface-hopping simulation of non-adiabatic vibration-induced autoionization of molecular anions

## Installation

#### Prerequisites:
- [Cython](https://github.com/cython/cython) is needed to compile modules used in the calculation of diabatic couplings (necessary) and exact two-electron integrals with one plane wave (highly experimental, default: off)
- A version of [Gaussian09/Gaussian16](https://gaussian.com/) or [QChem](https://www.q-chem.com/) is needed for the quantum-chemical calculations of each timestep.

#### Dependencies:
The following packages are needed to run the program and downloaded during installation:
- [scipy](https://github.com/scipy/scipy), mainly for the integration of the electronic Schrödinger equation
- [pyscf](https://github.com/pyscf/pyscf) for the calculation of two-electron integrals
- [joblib](https://github.com/joblib/joblib) for the parallelization of the calculation of diabatic couplings
- [matplotlib](https://github.com/matplotlib/matplotlib) for the plots produced by aid_latest/gui/guiAnalysis.py

#### Actual installation

In the program folder use the commands

    $ python cysetup.py build_ext --inplace  
    $ pip install .  

to first compile the Cython files and then pip-install the program.

#### Deinstallation

    $ pip uninstall hortensia_latest  



## How to use

After installation, the dynamics simulation can be started with with the command  

    $ hortensia --run  

Files needed to start the calculations are
- config.ini
- structure.in
- basis (atomic basis in Gaussian format)

If you wish to restart a trajectory, additional mandatory files are:
- dynamics.xyz
- velocities.dat
- tempcoef.dat
- trajpop.dat

#### Options

The program also has some useful functions beyond the dynamics simulation, a detailed overlook of which can be obtained by calling the help function with 

    $ hortensia --help  

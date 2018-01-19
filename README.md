# hmcwithlammps

# Python scripts for performing hybrid Monte Carlo (HMC) with LAMMMPS.  

Updated: January 19, 2018

Notes:
- A manuscript documenting the code development and validation strategy is currently in preparation
- The code was developed and tested using LAMMPS-17Nov16 release  
- The Python scripts require LAMMPS to be built as a shared library (see http://lammps.sandia.gov/doc/Section_start.html)

Installing LAMMPS as a shared library (LAMMPS-17Nov16 release): 
1.  Build LAMMPS as a shared library
    - cd /pathtolammps/src
    - make machine mode=shlib, where machine specifies the appropriate makefile located at /pathtolammps/src/MAKE/Makefile.machine.The Makefile.machine should specify the appropriate FFTW and MPICH libraries; FFTW and MPICH must also be built as shared libraries
2. Set environmental variables to add /pathtolammps/python/lammps.py and /pathtolammps/src/liblammps.so to the appropriate paths. For bash shell, this can be accomplished by add the following lines to ~/.bashrc:
    - export PYTHONPATH=/pathtolammps/python:$PYTHONPATH
    - export LD_LIBRARY_PATH=/pathtolammps/src:$LD_LIBRARY_PATH

Examples:
- LJ_Argon_132K_NVT:  Canonical ensemble HMC simulation of Lennard-Jones argon at 132 K and a density of 0.451 g/cc
- LJ_Argon_125K_100atm_NPT: Isothermal-isobaric ensemble HMC simulation of Lennard-Jones argon at 125 K and 100 atm 
- LJ_Argon_Checkensemble_NVT: Canonical ensemble HMC simulations of Lennard-Jones argon performed at 132 K and 145.5 K using a fixed density of 0.451 g/cc.  Data from the simulations was used to check the ensemble consistency of the code via the potential energy test described in  J. Chem. Theory Comput. 9, 909-926 (2013). The consistency test was performed using the checkensemble package [https://github.com/shirtsgroup/checkensemble]. The directory includes the original simulation input and output data, scripts to format the output data for compatibility with checkensemble, and the results from the checkensemble test. 
- mW_Water_WLbias:  Biased HMC simulation of the mW water model [Molinero and Moore, J. Phys. Chem. B 113, 4008 (2008)] at 220 K and 1 bar.  The Wang-Landau algorithm is used to apply an adaptive bias to the system to drive homogeneous ice nucleation.  The bias is applied to the size of the largest crystalline cluster defined using the criteria of Reinhardt and Doye [J. Chem. Phys. 136, 054501 (2012)]. The cluster analysis is performed using a supporting Fortran library (oplib.f90); the Python interface to the Fortran library must be built using F2PY [https://github.com/pearu/f2py]. 


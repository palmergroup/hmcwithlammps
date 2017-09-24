# hmcwithlammps

# Python scripts for performing hybrid Monte Carlo (HMC) in the isothermal-isobaric ensemble with LAMMMPS.  

Updated: September 24, 2017

Notes:
- A manuscript documenting the code development and validation strategy is currently in preparation
- The code was developed and tested using LAMMPS-17Nov16 release  
- The Python scripts require LAMMPS to be built as a library (see http://lammps.sandia.gov/doc/Section_start.html)

Examples:
- Example LJ_Argon_132K_NVT:  Canonical ensemble HMC simulation of Lennard-Jones argon at 132 K and a density of 0.451 g/cc
- Example LJ_Argon_125K_100atm_NPT: Isothermal-isobaric ensemble HMC simulation of Lennard-Jones argon at 125 K and 100 atm 
- Example mW_Water_WLbias:  Biased HMC simulation of the mW water model [Molinero and Moore, J. Phys. Chem. B 113, 4008 (2008)] at 220 K and 1 bar.  The Wang-Landau algorithm is used to apply an adaptive bias to the system to drive homogeneous ice nucleation.  The bias is applied to the size of the largest crystalline cluster defined using the criteria of Reinhardt and Doye [J. Chem. Phys. 136, 054501 (2012)]. The cluster analysis is performed using a supporting Fortran library (oplib.f90); the Python interface to the Fortran library must be built using F2PY [https://github.com/pearu/f2py]. 


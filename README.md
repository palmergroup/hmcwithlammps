# hmcwithlammps

# Python scripts for performing hybrid Monte Carlo (HMC) with LAMMMPS.  

Updated: September 24, 2017

Notes:
- A manuscript documenting the code development and validation strategy is currently in preparation
- The code was developed and tested using LAMMPS-17Nov16 release  
- The Python scripts require LAMMPS to be built as a library (see http://lammps.sandia.gov/doc/Section_start.html)

Examples:
- LJ_Argon_132K_NVT:  Canonical ensemble HMC simulation of Lennard-Jones argon at 132 K and a density of 0.451 g/cc
- LJ_Argon_125K_100atm_NPT: Isothermal-isobaric ensemble HMC simulation of Lennard-Jones argon at 125 K and 100 atm 
- LJ_Argon_Checkensemble_NVT: Canonical ensemble HMC simulations of Lennard-Jones argon performed at 132 K and 145.5 K using a fixed density of 0.451 g/cc.  Data from the simulations was used to check the ensemble consistency of the code via the potential energy test described in  J. Chem. Theory Comput. 9, 909-926 (2013). The consistency test was performed using the checkensemble package [https://github.com/shirtsgroup/checkensemble]. The directory includes the original simulation input and output data, scripts to format the output data for compatibility with checkensemble, and the results from the checkensemble test. 
- mW_Water_WLbias:  Biased HMC simulation of the mW water model [Molinero and Moore, J. Phys. Chem. B 113, 4008 (2008)] at 220 K and 1 bar.  The Wang-Landau algorithm is used to apply an adaptive bias to the system to drive homogeneous ice nucleation.  The bias is applied to the size of the largest crystalline cluster defined using the criteria of Reinhardt and Doye [J. Chem. Phys. 136, 054501 (2012)]. The cluster analysis is performed using a supporting Fortran library (oplib.f90); the Python interface to the Fortran library must be built using F2PY [https://github.com/pearu/f2py]. 


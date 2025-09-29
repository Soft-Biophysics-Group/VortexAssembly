# VortexAssembly
Codes and executable to run defect engineered vortex assemblies introduced in https://arxiv.org/abs/2504.13073


This package contains
- A notebook VortexAssembly which can be used to run and visualize simulations
- A C++ executable Sys which runs a simulated annealing from a parameter files and produces textfiles containing the system configuration after annealing
- Extra python files LatticeTools.py System.py and SideClasses.py that help read and write the output of the executable
- A folder Lattice/ containing the text files to relate the interaction maps indexing with the relative orientations of two particles.
- A folder Results/ where the simulation results are saved

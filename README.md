# VortexAssembly_2D
Codes and executable to run defect engineered vortex assemblies introduced in https://arxiv.org/abs/2504.13073

This package contains
- A notebook VortexAssembly.ipynb which can be used to run and visualize simulations
- C++ files (.hpp and .cpp) and a makefile to compile to obtain an executable Sys. This executable runs a simulated annealing from a parameter files and produces textfiles containing the system configuration after annealing
- Extra python files LatticeTools.py System.py and SideClasses.py that help read and write the output of the executable
- A folder Lattice/ containing the text files to relate the interaction maps indexing with the relative orientations of two particles.
- A folder Results/ where the simulation results are saved


To use it:
1/ Download all the files and folders (make sure that the subolders Lattice/ and Results/ are there 
2/ Compile the C++ files by using the command "make" from your terminal, in the main folder. This creates the Sys executable
3/ Use the VortexAssembly.ipynb to run simulation and visualize the results

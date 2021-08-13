README for hybrid RDME-gCME-ODE minimal cell program and simulations

These programs require the computer to be equipped with GPUs that are supported by CUDA 10 or greater.  The simulations by Thornburg et al., 2021 were run on an NVIDIA Titan V and NVIDIA Tesla Volta V100

To run these simulations, when installing Lattice Microbes (https://github.com/zanert2/Lattice_Microbes/) the following line needs to be executed at the step of making with cmake:

cmake path/Lattice-Microbes/src/ -D MPD_GLOBAL_T_MATRIX=True -D MPD_GLOBAL_R_MATRIX=True

Once Lattice Microbes has been installed, install odecell in the same environment

Folders:
model_data: Contains all of the kinetic parameters and initial condition files to initialize spatial Minimal Cell simulations

program: Contains the jupyter notebook to run spatial Minimal Cell simulations as well as the python modules used in the hybrid algorithm.

simulations: Contains results of spatial Minimal Cell simulations and the analysis jupyter notebook

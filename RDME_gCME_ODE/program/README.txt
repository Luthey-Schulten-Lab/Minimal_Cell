Guide to hybrid RDME program files

MinCell_RDME_CME_ODE.ipynb: Jupyter notebook that executes the hybrid simulation. The notebook will run when opened and simply "Restart and Run All". The results get saved in ../simulations/ and there is an option to name the sim folder that will be created at the top of the notebook. 

defMetRxns.py: Defines metabolic reactions for ODE part of the simulation

diffusion.py: Defines basic diffusion ruleas for the simulation and defines functions to calculate diffusion rates

GIP_rates.py: Defnies functions to calculate elongation rates for genetic information processing reactions

hook.py: Contains the hybrid algorithm that executes all solvers and calls all communication functions

in_out.py: Defines the communication functions between methods

integrate.py: Calls numeric integrator for ODE part of the simulation

lipid_patch.py: Modifications to lipid reaction parameters. (These are the parameters reported in the supplementary table of the associated paper (Thornburg et al., 2021))

MC_CME.py: Writes the reactions for the global CME part of the simulation

MC_RDME.py: Adds all particles, their diffusion rules, and all reactions to the spatial simulation

PPP_patch.py: Modifications to PPP pathway parameters and a few other reactions. (These are the parameters reported in the supplementary table of the associated paper (Thornburg et al., 2021))

regions_and_complexes.py: Constructs cellular architecture and adds degradosomes, ribosomes, and SecY

run_CME.py: Executes global CME step of the hybrid algorithm

Rxns_ODE.py: Defines functions used to write the rates of the metabolic ODE reactions

setICs.py: Sets initial conditions for metabolites

Simp_ODE.py: Initializes ODE simulation that will be passed to the integrator

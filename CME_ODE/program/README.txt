Date: 6/14/21
Author: David Bianchi
Luthey-Schulten Lab (University of Illinois Urbana-Champaign)

README for running hybrid CME-ODE JCVI-syn3A simulations - with multiple replication initiation events - restart simulations (dynamic gene expression rates)
-------------------------------------------------------------------------------------------------------------------------------------------------------------

NOTE: First ensure that you have the odecell module, used to write the ODE equations for metabolism, installed via ../odecell (one level up).

REMINDER: Useful Lattice Microbes Software Suite Documentation and Installation Guide is Available at - https://github.com/zanert2/Lattice_Microbes


Single Cell Notebook: Sample Jupyter Notebook to run and analyze a simulation of a single cell
===================================================================================================
serialReplicate-ExampleNotebook.ipynb - runs and performs an example analysis of a single simulated in silico JCVI-syn3A cell.



To Run Simulations in Parallel
-----------------------------------
Install Required Packages (to your conda environment)
======================================================
conda install -c anaconda mpi4py (Necessary to install dependencies required for parallelization routine)

Bash run file: mpi_runs.sh (to run 5 cell replicates)
======================================================
To run give the command: ./mpi_runs.sh

Which will run the following:

mpirun -np 5 python3 ~/projects/minCell_CMEODE/restart_poly/mpi_wrapper.py -st cme-ode -t 125 -rs 1

where:

-np: number of cells/replicates: an int

-st: simulation type = cme-ode (default)

-t: simulation time (minutes) = 125 (default)

-rs: restart time (minutes) = 1 (default) - the time period at which gene expression reaction rates are reset with updated metabolite pools (i.e. NTP pools)


Command to run simulations:
=============================

"./mpi_runs > log.log" 

nohup: "nohup ./mpi_runs > nohup.log &"


Output of simulations:
=============================
Upon running a simulation, the following files will be generated in the /simulations directory:

1) a .lm file (HDF5) : this file contains the simulation data (chemical species time traces in HDF5 format), these files can be opened and analyzed via PYSTDLM OpenLMFile functionality (Documentation available at: https://github.com/zanert2/Lattice_Microbes)

2) a .csv file (viewable in Excel): this file contains for each timestep a comma separated list of all of the species and their particle values (at each time, i.e. 1 min., 2 min., ...)

3) a .log file (text file): this file contains the simulation log data for the first timestep of simulation time. This can be useful in debugging if issues occur. The verbosity of the details written to this file can be adjusted via changing the LM.LoggerLevel() command (Explained in the Lattice Microbes Documentation available at: https://github.com/zanert2/Lattice_Microbes).

4) a flux file (located in simulations/fluxes/): this file contains the simulation reaction flux data for all metabolic reactions at each timestep. ***This file can be analyzed*** - TODO: Add an example of flux analysis.



Program Structure:
==============================



Main Simulation Input File(s):
-------------------------------

MinCell_CMEODE_mpi_.py: runs the initial simulation setup and first minute of simulation time. (Saves .lm file containing simulation data, inluding number of species and reactions etc.)

MCrestartLoop.py: Iteratively restarts the simulation at an interval of 1 minute so that CME Gene Expression
Reactions are updated by ODE MEtabolism output.



Hook File(s):
------------------------------
hook_(re)start.py - Communicates/hooks the CME Simulation of Genetic Information Process to the ODE Simulation of Metabolic Reactions.



Define Reactions:
-----------------------------

defLipCentNuc_two.py: defines the metabolic reactions to be integrated as ODEs, which are input via the odecell module.



I/O handling:
----------------------------

in_out_two.py: Handles Input and Output to .csv files of simulation particle numbers and reaction fluxes as well as cell growth tracking and particle number to concetration conversion.

Modifications to Reaction Parameters:

PPP_patch.py: Modifies .tsv original reaction parameters for the Pentose Phosphate Pathway reactions.

lipid_patch.py: Modifies .tsv original reaction parameters for lipid metabolism/cell growth reactions.



Replication Initiation:
--------------------------

rep_start.py: Defines the initial setup and conditions for reactions for multiple replication initiation events.

rep_restart.py : Defines the gene replication reactions for multiple replication initiation events (with cell volume aware dnaA filament formation for replication initiation).



Integration:
-------------------------

integrate.py: Defines the ODE Integration Model with fixed timestepping.

Patches:
------------------------
PPP_patch.py - "Patch" that updates certain Pentose Phosphate Pathway parameters from those obtained through the conventional kinetic parameter gathering pipeline.

lipid_patch.py - "Patch" that updates certain Lipid Metabolism parameters from those obtained through the conventional kinetic parameter gathering pipeline.

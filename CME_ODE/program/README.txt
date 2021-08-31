Date: 6/14/21
Author: David Bianchi
Luthey-Schulten Lab (University of Illinois Urbana-Champaign)

README for running hybrid CME-ODE JCVI-syn3A simulations - with multiple replication initiation events - restart simulations (dynamic gene expression rates)
-------------------------------------------------------------------------------------------------------------------------------------------------------------

Bash run file: mpi_runs.sh
============================
mpirun -np 65 python3 ~/projects/minCell_CMEODE/restart_poly/mpi_wrapper.py -st cme-ode -t 125 -rs 1

where:

-np: number of cells/replicates: an int

-st: simulation type = cme-ode (default)

-t: simulation time (minutes) = 125 (default)

-rs: restart time (minutes) = 1 (default) - the time period at which gene expression reaction rates are reset with updated metabolite pools (i.e. NTP pools)


Command to run simulations:
=============================

"./mpi_runs > log.log" 

nohup: "nohup ./mpi_runs > nohup.log &"


Program Structure:
==============================


Main Simulation Input File(s):
-------------------------------

MinCell_CMEODE_mpi_.py: runs the initial simulation setup and first minute of simulation time. (Saves .lm file containing simulation data, inluding number of species and reactions etc.)

MCrestartLoop.py: Iteratively restarts the simulation at an interval of 1 minute so that CME Gene Expression
Reactions are updated by ODE MEtabolism output.


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

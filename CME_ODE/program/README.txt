6/14/21
Author: David Bianchi

README for running "Restart" Simulations - WITH MULTIPLE REPLICATION INITIATION EVENTS - of the CME-ODE Model for JCVI-Syn3A
----------------------------------------------------------------------------------

Bash run file: mpi_runs.sh
============================
mpirun -np 65 python3 ~/projects/minCell_CMEODE/restart_poly/mpi_wrapper.py -st cme-ode -t 125 -rs 1

where:

-np: number of cells/replicates: an int

-st: simulation type = cme-ode (default)

-t: simulation time (minutes) = 125 (default)

-rs: restart time (minutes) = 1 (default)


Command to run simulations:
=============================

"./mpi_runs > log.log" 

nohup: "nohup ./mpi_runs > nohup.log &"


Program Structure:
==============================

Main Simulation Input File(s):

MinCell_CMEODE_mpi_.py: runs the initial simulation setup and first minute of simulation time. (Saves .lm file containing simulation data, inluding number of species and reactions etc.)

MCrestartLoop_.py: Iteratively restarts the simulation at an interval of 1 minute so that CME Gene Expression
Reactions are updated by ODE MEtabolism output.

Define Reactions:

defLipCentNuc_two.py: defines the metabolic reactions to be integrated as ODEs, which are input via the odecell module.

I/O handling:

in_out_two.py: Handles Input and Output to .csv files of simulation particle numbers and reaction fluxes as well as cell growth tracking and particle number to concetration conversion.

Modifications to Reaction Parameters:

PPP_patch.py: Modifies .tsv original reaction parameters for the Pentose Phosphate Pathway reactions.

lipid_patch.py: Modifies .tsv original reaction parameters for lipid metabolism/cell growth reactions.

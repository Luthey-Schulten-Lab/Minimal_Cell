# Minimal_Cell

This directory contains the simualtion files for the well-stirred and spatial models of the Minimal Cell JCVI-syn3A by Thornburg et al., 2021. Also included is the directory to install the package odecell which is used in the simulations to write the ODE equations for the metabolism.

CME_ODE: Contains program for well-stirred model of JCVI-syn3A capable of simualting a whole cell cycle.  The modified FBAm is included in model_data/FBA/

RDME_CME_ODE: Contains program for spatially resolved model of the first 20 minutes of the cell cycle for JCVI-syn3A

odecell: Package required for metabolic reactions.  Writes the reaction rates into a system of ODEs that can be given to a solver

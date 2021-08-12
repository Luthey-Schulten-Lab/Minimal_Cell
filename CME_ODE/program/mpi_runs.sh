#!/bin/bash

## Try 5 CME-ODE replicates for 4 minutes
mpirun -np 30 python3 ~/projects/minCell_CMEODE/restart_poly_dbl/mpi_wrapper.py -st cme-ode -t 122 -rs 1

#wait

#!/bin/bash

## Test of 50 replicates
mpirun -np 60 python3 ~/projects/minCell-CMEODE/central_hybrid/mpi_wrapper.py -st cme-ode

wait

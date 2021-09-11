#!/bin/bash

## Try 5 CME-ODE replicates for 5 minutes with a restart time of 1 minute
mpirun -np 5 python3 mpi_wrapper.py -st cme-ode -t 5 -rs 1

#wait

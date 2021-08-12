"""
Author: David Bianchi
Date: 12/10/2020
"""

"""
A wrapper script for MPI jobs on compute clusters
"""

# Import necessary modules

import os

from contextlib import redirect_stdout

import subprocess

from mpi4py import MPI

import argparse

import time

import sys

import lm

# Set PATH to CWD
sys.path.append(os.getcwd())


# Runs a CME-ODE simulation
def runCMEODE(rank, simTime, restartTime):

        # NOTE: AVOID USING A SUBPROCESS HERE, INSTEAD USE OS TO CALL THE BASH SCRIPT FOR DESIRED NUMBER OF REPLICATES

        # Run a hybrid CME-ODE simulation
        #run = subprocess.Popen(["python3","LipCentNucGIP_CMEODE_counters.py","-r",("{0}").format(int(rank)+1)], shell=False)
        #run = subprocess.Popen(["python3","MinCell_CMEODE.py","-r",("{0}").format(int(rank)+1)], shell=False)

        # Wait for the process to complete
        #run.wait()
        #os.("./mpi_test")

        # Launch the process: It doesn't actually provide 2 replication forks at the moment so this is fine
        run = subprocess.run(['python3',"MinCell_CMEODE_mpi_two_TwoRep.py","-procid",("{0}").format(int(rank)+1),"-t",("{0}").format("1")],shell=False)

        #run.wait()
        #run.close()

        runTwo = subprocess.run(['python3',"MCrestartLoop_twoRep.py","-procid",("{0}").format(int(rank)+1),"-t",("{0}").format(simTime),"-iter","1","-rs",("{0}").format(restartTime)],shell=False)
        # Wait for the process to complete - Not necessary
        #runTwo.close()

        return


# Main Simulation Function
# Take in the user input and enact simulation launch
def main():

        # Parse simulation arguments/parameters
        ap = argparse.ArgumentParser()

        ap.add_argument('-st', '--simType', required = True)
        ap.add_argument('-t', '--simTime', required= True)
        ap.add_argument('-rs', '--restartTime', required= True)
        args = ap.parse_args()

        # Get the MPI rank
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        # Get the simulation time
        try:
            simTime = int(args.simTime)
            restartTime = int(args.restartTime)
        except:
            print("Please enter simulation time (-t) as the integer valued simulation time in minutes.")

        # Setup and run a hybrid CME-ODE simulation
        if (str(args.simType) == "cme-ode"):
                runCMEODE(rank,simTime,restartTime)

        # There has been an error in user input
        else:
                print("Incorrect simulation type input")
                print("Enter 'cme-ode' for a hybrid CME-ODE simulation")

# Run main() if not importing as a module
if __name__ == "__main__":
        main()

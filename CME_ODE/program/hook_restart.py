"""
The hook simulation driver for a Hybrid CME-ODE JCVI Syn3A simulation.

Author: David Bianchi
"""

import Simp as Simp
import Rxns as Rxns
import integrate as integrate
import copy as copy
import in_out as in_out
import sys
import pandas as pd
import numpy as np
import time as timer
import lm as lm

### Define our own hybrid CME-ODE solver class derived from the LM Gillespie Direct Solver:
class MyOwnSolver(lm.GillespieDSolver):


        """
        Initialize the ODE hook solver

        @param self The object pointer
        @param f Array containing the species counts of the ODE species
        @param delt The communication timestep between the hook and main LM simulation
        @param odestep The maximum stepsize given to the Adaptive Timestepping ODE solver
        @param speciesCount The instance of SpeciesCount Class used to pass species counts 
        """
        #def initializeSolver(self, f, delt, ode_step, speciesCount):
        def __init__(self, f, delt, ode_step, speciesCount,cythonBool,totalTime,procID, iteration):

            """
            Initialize the ODE hook solver

            Parameters:
            self, the object pointer
            f (numpy array), contains the species counts of the ODE species
            delt (float), communication timestep between hook simulation and main LM simulation
            odestep (float), the maximum stepsize given an Adaptive Timestepping ODE solver
            speciesCount (species Count), instance of SpeciesCount Class used to pass species count data
            cythonBool (bool), Should ODE Reaction Solver be compiled with Cython (True, False)
            procID (str), The Process ID for each simulated "cell".

            Returns:
            None
            """

            # Call the constructor for the derived class
            # Not necessary to use MESolverFactory.setSolver('lm::cme::GillespieDSolver)?
            lm.GillespieDSolver.__init__(self)

            # Save the initial conditions, for restarting the solver upon a new replicate
            self.ic = (f,delt,ode_step,speciesCount,cythonBool,totalTime)

            # The time a which hook solver has been stepped into, initial value = 0
            self.oldtime = 0.0

            # The process ID for creating flux log files etc.
            self.procID = str(procID)

            self.iter = str(iteration)

            print("initializing solver")

            # Set the initial conditions
            self.restart()

        """
        Get the same initial conditions for a new replicate

        @param self The object pointer
        """
        def restart(self):

            """
            Get the same initial conditions for a new simulation replicate (Restart the Hook)

            Parameters:
            self, the object pointer

            Returns:
            None
            """
            
            # Set the previous time to be 0, we are starting the simulation
            self.oldtime = 0.0

            # Deep Copy of all of the initial conditions
            self.f = copy.deepcopy(self.ic[0])
            self.delt = copy.deepcopy(self.ic[1])
            self.odestep = copy.deepcopy(self.ic[2])
            self.species = copy.deepcopy(self.ic[3])
            self.cythonBool = copy.deepcopy(self.ic[4])
            self.totalTime = copy.deepcopy(self.ic[5])

            # Update need enzyme Counts in the particle map
            self.species.update(self)

            print("Done with restart")

        
        
        """
        The hookSimulation method defined here will be called at every frame write
        time.  The return value is either 0 or 1, which will indicate if we
        changed the state or not and need the lattice to be copied back to the GPU
        (In the case of the RDME) before continuing.  If you do not return 1, 
        your changes will not be reflected.

        @param self The object pointer
        @param time The current simulation time
        """
        def hookSimulation(self, time):

                """
                The hookSimulation method defined here will be called at every frame write
                time.  The return value is either 0 or 1, which will indicate if we
                changed the state or not and need the lattice to be copied back to the GPU
                (In the case of the RDME) before continuing.  If you do not return 1, 
                your changes will not be reflected.

                Parameters:
                self, the object pointer
                time, the current simulation time

                Returns:

                1 (int), if changes should be passed to the main LM Simulation
                0 (int), if changes should not be passed to the main lm Simulation
                """

                # We have reached the simulation start time, if doing multiple replicates
                # No need to update
                if (time==0.0):
                    print("New Replicate", flush=True)
                    self.restart()
                    minute = 0
                    return 0

                # We are at a CME-ODE communication timestep
                else:

                    # At the first timestep update the needed protein counts
                    if ((time > self.delt) and (time < (self.delt*2.0))):
                        self.species.update(self)
                        #Simp.upIC(self.species)
                        #Simp.upIC(self.species)

                    # Update to current solver species counts
                    start = timer.time()
                    print("Updating species: ", start)
                    self.species.update(self)
                    end = timer.time()
                    print("Finished update: ",end)
                    print("Time is: ",time)

                    # Initialize and define the reaction model
                    model = Simp.initModel(self.species)


                    ### Want to get the current values, not necessarily the initial values
                    initVals=integrate.getInitVals(model)

                    ### Boolean control of cython compilation, versus scipy ODE solvers
                    cythonBool = self.cythonBool

                    if (cythonBool == True):
                        solver=integrate.setSolver(model)
                    
                    else:
                        solver=integrate.noCythonSetSolver(model)

                    ### Run the integrator: But instead of passing self.delt pass self.oldtime
                    res = integrate.runODE(initVals,time,self.oldtime,self.odestep,solver,model)

                    resFinal = res[-1,:]
                    
                    resStart = res[0,:]
                    
                    if (int(time)/100).is_integer():
                        print('Progress: ' + str(int(time)) + ' out of ' + str(int(self.totalTime)))
                    
                    if (int(time)/60).is_integer():
                        minute = int(int(time)/60)
                        currentFluxes = solver.calcFlux(0, resStart )

                        # Create list of reactions and fluxes
                        fluxList = []
                        for indx,rxn in enumerate(model.getRxnList()):
                            fluxList.append( (rxn.getID(), currentFluxes[indx]) )

                        fluxDF = pd.DataFrame(fluxList)

                        fluxFileName = './simulations/fluxes/' + 'rep-' + self.procID + '-fluxDF.csv' #'/fluxDF_'+str(self.iter)+'min_start.csv'

                        fluxDF.to_csv(fluxFileName,header=False,mode='a')
                        
                        minute = int(int(time)/60)
                        currentFluxes = solver.calcFlux(0, resFinal )

                        # Create list of reactions and fluxes
                        fluxList = []
                        for indx,rxn in enumerate(model.getRxnList()):
                            fluxList.append( (rxn.getID(), currentFluxes[indx]) )

                        fluxDF = pd.DataFrame(fluxList)

                        fluxFileName = './simulations/fluxes/' + 'rep-' + self.procID + '-fluxDF-end.csv' #'/fluxDF_'+str(self.iter)+'min_end.csv'

                        fluxDF.to_csv(fluxFileName,header=False,mode='a')
                       
 
                        fluxFileName = './simulations/fluxes/' + 'rep-' + self.procID + '-fluxDF.csv' #'/fluxDF_'+str(self.iter)+'min.csv'

                        fluxDF.to_csv(fluxFileName,header=False,mode='a')

                        print('Saved fluxes at ' + str(minute) + ' minutes.')
                    
                        print('Saved final fluxes.')


                    if time > (self.totalTime-self.delt):
                        print(time)
                        minute = int(int(time)/60)
                        finalFluxes = solver.calcFlux(0, resFinal )

                        # Create list of reactions and fluxes
                        fluxList = []
                        for indx,rxn in enumerate(model.getRxnList()):
                            fluxList.append( (rxn.getID(), finalFluxes[indx]) )

                        fluxDF = pd.DataFrame(fluxList)
                        fnStr='./simulations/fluxes/'+ 'rep-' + self.procID + '-fluxDF_final.csv' #'/' + str(self.iter) + 'fluxDF_final.csv'
                        print("Writing Final Fluxes and Csvs for Restart")
                        fluxDF.to_csv(fnStr,index=False,header=False,mode='a')


                    # Get the previous time in minutes
                    minute = int(int(time)/60)
                    # Set the previous time to the current time
                    self.oldtime = time


                    # Write the results
                    in_out.writeResults(self.species,model,resFinal,time,self.procID)


                    # Update the system with changes
                    return 1

                return 0

                

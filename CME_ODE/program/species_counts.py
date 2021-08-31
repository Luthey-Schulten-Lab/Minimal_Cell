"""
Author: David Bianchi
Date: 5/19/17

A module for creating a class that will allow for easy passing of
species count data to and from the solvers used in CME and
hybrid CME-ODE simulations
"""

import numpy as np

"""
A class to allow for easy access to species counts in simulation.
"""
class SpeciesCounts():
        
        """
        Constructor
       
        Parameters:

               self - The object pointer
               sim - The simulation object
        """
        def __init__(self, sim):
                
                # A numpy array will hold the counts, initially empty
                self.count_array = None

                # A map of species -> particle number stored in lattice microbes
                self.particleMap = sim.particleMap


        """
        Updates the pointer to speciesCounts to the current system particle vales

        Parameters:
        
               self - The object pointer
               solver - The solver being used by LM-CME (Gillespie Direct etc.)
        """
        def update(self, solver):
            
                # Get a numpy array containing particle counts for each species
                self.count_array = solver.getSpeciesCountView()

        """
        A getter for the specie count of a give specie.
        Overloading of rhs indexing.
       
        Parameters:

               self - The object pointer
               key - The name/identifier of a specie
        """
        def __getitem__(self,key):
                
                # Adjustment since indices start at 1 in LM main
                idx = self.particleMap[key] - 1
                
                # Return a pointer to the select specie
                return self.count_array[idx]

        """
        A setter for the specie count of a given specie.
        Overloading lhs indexing.
        
        Parameters:

               key - The name of specie
               val - The amount of the new count
        """
        def __setitem__(self,key,val):

                # Adjustment since indices start at 1 in LM main
                idx = self.particleMap[key] - 1

                # Set the value of the pointer to a given species
                self.count_array[idx] = val

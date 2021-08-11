##### HOOK #####

import lm as lm
import os
import random
import time as tm
from pyLM import CME
# from pyLM import RDME
from pyLM import LMLogger
from pyLM.units import *
import pySTDLM.PostProcessing as pp
import numpy as np

from jLM.RDME import Sim as RDMESim
from jLM.RDME import File as RDMEFile

from lm import GillespieDSolver

from pyLM import LMLogger
import logging
LMLogger.setLMLogConsole(logging.INFO)
lmlog=LMLogger.LMLogger

from pyLM import *
from pyLM.units import *
import math as math

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
import importlib
from collections import defaultdict, OrderedDict

# import matplotlib.pyplot as plt

from tqdm import tqdm

# from GIP import *
from diffusion import *
import MC_CME_polysomes
from MC_RDME_polysomes import * #_polysomes
# from cellgeometry import *
# from GIP_rates import *
import Simp_ODE as Simp
import Rxns_ODE as Rxns
import integrate as integrate
import copy as copy
import in_out_polysomes as in_out #_polysomes
import sys
import pandas as pd

### Define our own solver class derived from IntMpdRdmeSolver
class MyOwnSolver:
    
    def __init__(self, lmFile, simFolder, delt, odestep, cythonBool, pmap, totalTime, geneEnds, geneStarts, singleStatePtnDict, multiStatePtnDict, degDict, tRNAstateDict, RDME_species_list, PartIdxMap, rtRNA_ID_dict, ordered_poly_ribo):
        
        super(MyOwnSolver, self).__init__()

        # Call the constructor for the derived class
        # Not necessary to use MESolverFactory.setSolver('lm::cme::GillespieDSolver)?
#             lm.IntMpdRdmeSolver.__init__(self)

        # Save the initial conditions, for restarting the solver upon a new replicate
        self.ic = (delt, odestep, cythonBool, pmap, totalTime, geneEnds, geneStarts, singleStatePtnDict, multiStatePtnDict, degDict, tRNAstateDict)
        
        if isinstance(lmFile, (RDMEFile, RDMESim)):
            self.rdme = lmFile
        else:
            self.rdme = RDME.File(lmFile)

        # The time a which hook solver has been stepped into, initial value = 0
        self.oldtime = 0
        self.geneEnds = geneEnds
        self.geneStarts = geneStarts
        self.singleStatePtnDict = singleStatePtnDict
        self.multiStatePtnDict = multiStatePtnDict
        self.degDict = degDict
        self.tRNAstateDict = tRNAstateDict
        self.RDME_species_list = RDME_species_list
        self.PartIdxMap = PartIdxMap
        self.rtRNA_ID_dict = rtRNA_ID_dict
        self.pmap = pmap
        self.simFolder = simFolder
        self.ordered_poly_ribo = ordered_poly_ribo

        # Set the initial conditions
        self.restart()

        print("Initialized solver")
        
    def restart(self):

#             """
#             Get the same initial conditions for a new simulation replicate (Restart the Hook)

#             Parameters:
#             self, the object pointer

#             Returns:
#             None
#             """
            
#             # Set the previous time to be 0, we are starting the simulation
            self.oldtime = 0

#             # Deep Copy of all of the initial conditions
#             self.f = copy.deepcopy(self.ic[0])
            self.delt = copy.deepcopy(self.ic[0])
            self.odestep = copy.deepcopy(self.ic[1])
#             self.species = copy.deepcopy(self.ic[3])
            self.cythonBool = copy.deepcopy(self.ic[2])
            self.totalTime = copy.deepcopy(self.ic[4])

#             # Update need enzyme Counts in the particle map
#             #self.species.update(self)
#             #Simp.upIC(self.species)
#             #Simp.upIC(self.species)

            print("Done with restart")

    
    # Timestep Variables
#     tsnum=0
#     curt=0.0
#     delt=1.0

#     pbar = None

    # The hookSimulation method defined here will be called at every frame 
    # write time.  The return value is either 0 or 1, which will indicate 
    # if we changed the state or not and need the lattice to be copied back 
    # to the GPU before continuing.  If you do not return 1, your changes 
    # will not be reflected.
    def hookSimulation(self, t, lattice):
        
        time = t
        print('Delt',self.delt)
        print('ODEstep',self.odestep)
        
        geneEnds = self.geneEnds
        geneStarts = self.geneStarts
        singleStatePtnDict = self.singleStatePtnDict
        multiStatePtnDict = self.multiStatePtnDict
        degDict = self.degDict
        tRNAstateDict = self.tRNAstateDict
        RDME_species_list = self.RDME_species_list
        PartIdxMap = self.PartIdxMap
        rtRNA_ID_dict = self.rtRNA_ID_dict
        pmap = self.pmap
        simFolder = self.simFolder
        
        ordered_poly_ribo = self.ordered_poly_ribo
        
        if np.rint(time) == 0:
            
            print(simFolder)
            
            try:
                csim_folder=simFolder+'cme_sims/'
                print(csim_folder)
                os.makedirs(csim_folder)
                print('Created global CME directory')
            except:
                print('CME sim directory already exists')
                
            try:
                flux_folder = simFolder+'fluxes/'
                print(flux_folder)
                os.makedirs(flux_folder)
                print('Created fluxes directory')
            except:
                print('Fluxes directory already exists')
        
        lmlog.info("Hook at time: %f sec"%time)
        lmlog.info("Creating CME simulation...")
#         curtime=time
        
        #### Update pmap based on RDME counts
        
        RDME_cts = self.rdme.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())
        
        in_out.calc_RDME_costs(pmap, RDME_cts, self.rdme, degDict, singleStatePtnDict, multiStatePtnDict)
        
#         RDME_cts = self.rdme.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())
        
        in_out.RDME_to_pmap(pmap, RDME_cts, self.rdme, multiStatePtnDict, tRNAstateDict, RDME_species_list)
        
#         RDME_cts = self.rdme.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())

        print('Updated particle counts from RDME')
        
        in_out.updatePolysomes(pmap,self.rdme,lattice,PartIdxMap,ordered_poly_ribo)
        
        print('Updated polysomes')
        
        # Create simulation
        csim=CME.CMESimulation()
        
        MC_CME_polysomes.constructCME(csim,pmap)
        
        # Set time data
        csim.setWriteInterval(0.1)
        csim.setSimulationTime(1.0)
        
        print('Starting Global CME')

        # Save and run simulation
        CSIMfilename= simFolder + 'cme_sims/cmeSim.%d.lm'%np.rint(time)
        print(CSIMfilename)
#         CSIMfilename='./cmeSim.lm'
        lmlog.info("Saving %s..."%CSIMfilename)
#         tm.sleep(1)
        try:
            os.remove(CSIMfilename)
        except:
            print('Nothing to delete')
        csim.save(CSIMfilename)
        lmlog.info("Running CME simulation...")
#         csim.run(CSIMfilename, "lm::cme::GillespieDSolver", 1)
#         os.system("lm -r 1 -ws -sl lm::cme::GillespieDSolver -f %s"%CSIMfilename)
#         lm.runSolver(CSIMfilename, 1, solver=GillespieDSolver(), cudaDevices=[0], checkpointInterval=0)

        os.system("python run_CME.py %s"%CSIMfilename)

#         
#         tm.sleep(1)
#         self.tsnum += int(1)

        print('Finished Global CME')

        # Read CME state
        lmlog.info("Postprocessing...")
        fHandle=pp.openLMFile(CSIMfilename)
        
        in_out.updateGlobalCmeCnts(pmap,CSIMfilename)
            
#         for key, val in pmap.items():
            
#             try:
#                 count_trace = pp.getSpecieTrace(fHandle, key)

#                 pcount = count_trace[-1]

#                 pmap[key] = pcount
                
#             except:
                
#                 continue
                
        print('Updated particle counts from global CME')
                
#         ts = pp.getTimesteps(fHandle)
#         tsShifted=[]
#         for i in range(len(ts)):
#             tsShifted.append(ts[i]+self.curt)
#         self.curt = curtime
#         self.times.extend(tsShifted)
        pp.closeLMFile(fHandle)
    
        try:
            os.remove(CSIMfilename)
            print('Deleted global CME file')
        except:
            print('Nothing to delete')
    
#         RDME_cts = self.rdme.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())
        
        ### Run ODE ###
        
        # Initialize and define the reaction model
        #model = lipSimp.initModel(passpMap)
        model = Simp.initModel(pmap)
        print('Initialized ODE simulation')

        ### Want to get the current values, not necessarily the initial values
        initVals=integrate.getInitVals(model)

        ### Boolean control of cython compilation, versus scipy ODE solvers
        cythonBool = self.cythonBool

        if (cythonBool == True):
            solver=integrate.setSolver(model)

        else:
            solver=integrate.noCythonSetSolver(model)

        ### Run the integrator
        res = integrate.runODE(initVals,time,self.delt,self.odestep,solver,model)

        resFinal = res[-1,:]

        resStart = res[0,:]
        
        print('Time:',time,np.rint(time)/100,np.rint(time)/60)

        if (np.rint(time)/100).is_integer():
            print('Progress: ' + str(np.rint(time)) + ' out of ' + str(int(self.totalTime)))

        in_out.calcODEfluxes(time,solver,model,resStart,resFinal,self.totalTime,self.delt,simFolder)
            
        print('Finished ODE simulation')
        # Set the previous time to the current time
        self.oldtime = time      
        
        ### Place new RNA transcribed from global CME into RDME ###
        in_out.placeNewRNA(pmap,self.rdme,lattice,geneEnds,geneStarts,PartIdxMap,rtRNA_ID_dict)

        ### Update Counts from ODE ###
        in_out.writeOdeResults(pmap,model,resFinal,self.rdme,lattice,geneEnds,geneStarts,PartIdxMap)
        
#         RDME_cts = self.rdme.particleStatistics(particleLattice=lattice.getParticleLatticeView(), siteLattice=lattice.getSiteLatticeView())
        print('Number of species',len(pmap))
        self.pmap = pmap
        
        print('Updated particle counts from ODE')
        
        if np.rint(time) == 0:
            
            particleDF = pd.DataFrame()
            
            spec_IDs = []
            
            new_counts = []
            
            for specID, count in pmap.items():
                
                spec_IDs.append(specID)
                
                new_counts.append(count)
                
            particleDF['Time']= spec_IDs
            particleDF[np.rint(time)] = new_counts
            
#             print(pmap)
            
            particleDF.to_csv(simFolder+'particle_counts.csv',index=False)
            
        else:
            
#             print(pmap)
            
            new_counts = []
            
            particleDF = pd.read_csv(simFolder+'particle_counts.csv')
            
            for specID, count in pmap.items():
            
                new_counts.append(count)
                
            particleDF[np.rint(time)] = new_counts
            
            particleDF.to_csv(simFolder+'particle_counts.csv',index=False)
        
        print('Saved global particle counts')
        
        
        in_out.updateGipCosts(pmap,model,resFinal,self.rdme,lattice,geneEnds,geneStarts,PartIdxMap)
        
        print('Updated GIP costs')

#         pbar.update(1)
        
        lmlog.info("Resuming RDME simulation...")
        return 1


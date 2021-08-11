"""
.. module:: paropt
   :platform: Unix
   :synopsis: A collection of functions to optimize parameters from ODE/FBA models.

.. moduleauthor:: Marcelo C. R. Melo <melomcr@gmail.com>

"""

from odecell import solver

from scipy.optimize import differential_evolution

import csv, collections, random, time, sys
import numpy as np
from math import ceil, isnan

from scipy import integrate
from pycvodes import integrate_adaptive

import cython
cimport numpy as np
ctypedef np.double_t DTYPEDBL_t

class Target():
    
    def __init__(self):
        
        self.__type = ""
        self.__ids = []
        self.__multipliers = []
        self.__idMultPairs = []
        self.__value = 0
        self.__error = 1
        self.__indxs = []
        self.__norm = False
        self.__weight = 1
    
    def __init__(self, newType, newIDs, newMult, newValue, newError, newNorm, newWeight):
        
        self.__type = ""
        self.__ids = []
        self.__multipliers = []
        self.__idMultPairs = []
        self.__value = 0
        self.__error = 1
        self.__indxs = []
        self.__norm = False
        self.__weight = 1
        
        self.loadData(newType, newIDs, newMult, newValue, newError, newNorm, newWeight)
    
    
    def loadData(self, newType, newIDs, newMult, newValue, newError, newNorm, newWeight):
        
        #print("New target read.")
        #print(newType, newValue, newError)
        #print(newIDs)
        #print(newMult)
        #print()
        
        if not newType in ["odeRxn","odeMet","fbaRxn"]:
            print("ERROR) Unknown target type:", newType)
            return -1
        
        self.__type = newType
        self.__value = float(newValue)
        self.__error = float(newError)
        self.__norm = newNorm
        self.__weight = float(newWeight)
        
        for ID in newIDs.split(";"):
            self.__ids.append(ID.strip())
        
        if not newMult:
            newMult = "1"
        
        for mult in newMult.split(";"):
            self.__multipliers.append( float(mult.strip()) )
        
        for i in range(len(self.__ids)):
            self.__idMultPairs.append( (self.__ids[i], self.__multipliers[i]) )
        
        #print("New target read.")
        #print(self.__type, self.__value, self.__error)
        #print(self.__ids)
        #print(self.__multipliers)
        #print()
    
    def getType(self):
        return self.__type
    
    def getVal(self):
        return self.__value
    
    def getError(self):
        return self.__error
    
    def getNorm(self):
        return self.__norm
    
    def getIDs(self):
        return self.__ids
    
    def getMultipliers(self):
        return self.__multipliers
    
    def getIDsMultsPairs(self):
        return self.__idMultPairs
    
    def getIndxs(self):
        return self.__indxs
        
    def getWeight(self):
        return self.__weight
    
    def addIndx(self, newIndx):
        self.__indxs.append(newIndx)
    
class ParOpt():
    
    def __init__(self, model, trgGrp="", normRxnID="", stdMethod="",
                 loopTime=0, dt=0, totalLoops=0, tol=0, useJac=False,
                 useFunctor= False, useCython=False, noBuild=False,
                 integrator="", vebose=0, buildVerb=0, logFN=""):
        
        self.__model = 0
        
        # 10 minutes
        self.__totalTime = loopTime
        
        self.__dt = dt
        
        self.__maxODEIter = totalLoops
        
        self.__tol = tol
        
        self.__withJac = useJac
        
        self.__useFunctor = useFunctor
        
        self.__cythonBuild = useCython
        
        self.__noBuild = noBuild
        
        self.__verbose = vebose
        
        # Verbosity level for the model build
        self.__buildVerb = buildVerb
        
        
        
        self.__recordFile = 0
        
        self.__OptLogFileName = "optEvolution.csv"
        
        self.__odeintFormat = 0
        
        self.__normRxnID = ""
        
        self.__trgGrpDict = dict()
        self.__trgList = list()
        
        self.__bounds = []
        
        self.__integratorID = 0
        
        self.__targetStand = 1
        
        self.__errReturn = 1E10
        
        self.__rxnIdList = 0
        
        self.__initParamVals = 0
        
        ##############################
        
        if (loopTime == 0 or
            dt == 0 or
            totalLoops == 0 or
            tol == 0 ):
            print("Invalid Input!")
            sys.exit(10)
        
        if self.__cythonBuild and vebose:
            print("\n\nATTENTION: Cython build could take more than a minute!\n\n")
        
        self.setIntegrator(integrator)
        
        self.setModel(model)
        
        if trgGrp:
            self.loadTargetGroup(trgGrp, normRxnID)
        
        if stdMethod:
            self.setTargetStandardizationMethod(stdMethod)
        
        self.setOptLogFileName(logFN)
        
    
    def setModel(self, newModel):
        
        self.__model = newModel
        
        if self.__useFunctor:
            self.solvFunctor = solver.ModelSolver(self.__model)
            self.solvFunctor.prepareFunctor()
            
            self.__rxnIdList = self.solvFunctor.buildCall(odeint=self.__odeintFormat, 
                                                   useJac=self.__withJac, 
                                                   verbose=self.__buildVerb, 
                                                   cythonBuild=self.__cythonBuild,
                                                   noBuild=self.__noBuild,
                                                   functor=True)
        
        self.__modelOptSpace, self.__initParamVals = self.__model.getOptSpace()
        
        #return self.__modelOptSpace, self.__initParamVals
    
    def setIntegrator(self, integrator):
        
        if integrator.lower() == "ode":
            self.__integratorID = 1
            self.__odeintFormat = False
        
        elif integrator.lower() == "odeint":
            self.__integratorID = 2
            self.__odeintFormat = True
            
        elif integrator.lower() == "cvodes":
            self.__integratorID = 3
            self.__odeintFormat = False
            
        else:
            print("Unknown integrator!!")
            sys.exit(1)
    
    def getIntegrator(self):
        return self.__integratorID
    
    def setTotalODETime(self, newTime):
        self.__totalTime = newTime
    
    def getTotalODETime(self):
        return self.__totalTime
    
    def setODEdT(self, newdT):
        self.__dt = newdT
    
    def getODEdT(self):
        return self.__dt
    
    def setMaxODEIter(self, newMaxIter):
        self.__maxODEIter = newMaxIter
    
    def getMaxODEIter(self):
        return self.__maxODEIter
    
    def setTargetStandardizationMethod(self, method):
        """ Defines how will model values be standardized for target comparison.
        
        Fluxes can be STANDARDIZED in two different ways when comparing them to
        target values.
        "error" will divide the squared difference between flux value and 
                target value by the "error" value provided for the target.
        "target" will cause the difference between flux value and target value
                to be divided by the target value itself, and the quotient 
                will be squared.
        
        In either case, if a reaction ID is supplied to NORMALIZE some or all
        reaction fluxes by that flux, then this NORMALIZATION will be carried
        out before any comparison to target values.
        
        """
        
        if method == "error":
            self.__targetStand = 1
        elif method == "target":
            self.__targetStand = 2
        else:
            print("\nWARNING: Unknown option for standardization method.")
    
    def setVerbose(self, newV):
        self.__verbose = newV
    
    def getVerbose(self):
        return self.__verbose
    
    def setUseJac(self, newUseJac):
        self.__withJac = newUseJac
    
    def getUseJac(self):
        return self.__withJac
    
    def setBuildVerb(self, verb):
        self.__buildVerb = verb
    
    def getBuildVerb(self):
        return self.__buildVerb
    
    def setNormRxnID(self, newRxnID):
        if isinstance(newRxnID, str):
            self.__normRxnID = [newRxnID]
        else:
            self.__normRxnID = newRxnID
    
    def getNormRxnID(self):
        return self.__normRxnID
    
    def setOptLogFileName(self, OptLogFileName):
        self.__OptLogFileName = OptLogFileName
    
    def getOptLogFileName(self):
        return self.__OptLogFileName
    
    def setCythonize(self, newCytho):
        self.__cythonBuild = newCytho
    
    def getCythonize(self):
        return self.__cythonBuild
    
    def getUseFunctor(self):
        return self.__useFunctor
    
    def getInitParams(self):
        return self.__initParamVals
    
    def clearTargets(self):
        self.__trgList = []
    
    def loadTarget(self, targetType, targetID, targetMult, targetVal, targetErr, 
                    targetNorm=True, targetWeight=1):
        """ Loads data for targets of optimization
            
            The available target types are odeRxn, odeMet, and fbaRxn.
            
        """
        
        self.__trgList.append( Target(targetType, 
                                      targetID, 
                                      targetMult, 
                                      targetVal, 
                                      targetErr,
                                      targetNorm,
                                      targetWeight) )
    
    def loadTargetGroup(self, fileName, normRxnID=""):
        
        #print("Parsing target file " + fileName + " ...")
        
        if normRxnID != "":
            self.setNormRxnID(normRxnID)
            normTargetFlux = True
        else:
            self.__normRxnID = 0
            normTargetFlux = False
        
        with open(fileName) as csvfile:
            
            # Defines the dictionarry creator.
            parDictReader = csv.DictReader(csvfile, delimiter=",")
            
            for row in parDictReader:
                
                targetType = row["type"]
                targetID = row["name"]
                targetMult = row["multiplier"]
                targetVal = row["value"]
                targetErr = row["error"]
                
                targetNorm = normTargetFlux
                if "normTarget" in row.keys():
                    parOpt = row["normTarget"].lower()
                    
                    if (parOpt == "on") or (parOpt == "yes"):
                        targetNorm = True
                    elif (parOpt == "off") or (parOpt == "no"):
                        targetNorm = False
                    else:
                        targetNorm = normTargetFlux
                
                if "weight" in row.keys():
                    targetWeight = row["weight"]
                else:
                    targetWeight = 1
                
                if targetID == "---":
                    continue
                
                if float(targetErr) == 0:
                    print("\n---> ERROR: Target error value cannot be zero!")
                    sys.exit(2)
                
                self.__trgList.append( Target(targetType, 
                                              targetID, 
                                              targetMult, 
                                              targetVal, 
                                              targetErr,
                                              targetNorm,
                                              targetWeight) )
        
        if self.__model:
            self.checkTargets()
    
    def checkTargets(self):
        
        for trgt in self.__trgList:
            
            #if trgt.getIDs()[0] == "---":
                #print("Ignoring target ---")
                #continue
            
            if trgt.getType() == "odeRxn":
                
                for id in trgt.getIDs():
                    
                    if self.__model.getReaction(id) == 1:
                        print("ERROR) odeRxn not found in model:", id)
                        return 1
            
            if trgt.getType() == "odeMet":
                
                for id in trgt.getIDs():
                    
                    metID, metIndex = self.__model.parseMetIdentifier(id)
                    
                    if metID == "":
                        print("ERROR) odeMet not found in model:", id)
                        return 1
                    else:
                        trgt.addIndx(metIndex)
    
    
    def createBounds(self):
        
        self.__bounds = []
        
        for optItem in self.__modelOptSpace:
            self.__bounds.append((float(optItem.lb), float(optItem.ub)))
        
        return self.__bounds
    
    def createBoundsUpLow(self):
        
        self.__boundsL = []
        self.__boundsU = []
        
        for optItem in self.__modelOptSpace:
            self.__boundsL.append(float(optItem.lb))
            self.__boundsU.append(float(optItem.ub))
        
        return self.__boundsL, self.__boundsU
    
    def recordEvolution(self, xk, convergence):
        
        self.__recordFile.write(str(xk) + "," + str(convergence) + str("\n") )
        self.__recordFile.flush()
        print(xk, convergence)
        
        return False
        #if convergence < 0.01:
            #return True
        #else:
            #return False
    
    def optimize(self, maxiter=1000, popsize=15, 
                 useJac = False, seed=0, verbosity = 0,
                 tol=0.01, mutation=(0.5,1), strategy="best1bin" ):
        
        self.setUseJac(useJac)
        self.setBuildVerb(verbosity)
        
        self.createBounds()
        
        self.__recordFile = open(self.__OptLogFileName,"w")
        self.__recordFile.write("xk,val\n")
        self.__recordFile.flush()
        
        if seed:
            sendSeed = seed
        else:
            sendSeed = None
        
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html
        res = differential_evolution(func=self.calcObjective,
                                      bounds=self.__bounds, 
                                      disp=True,
                                      seed = sendSeed,
                                      maxiter = maxiter,
                                      popsize = popsize,
                                      callback= None,
                                      mutation=mutation,
                                      tol=tol)
        
        self.__recordFile.close()
        
        return res
    
    def calcObjective(self, x):
        
        start_time = time.time()
        
        initVals = np.asarray(self.__model.getInitVals(), dtype=np.double)
        finalRes = np.empty(len(initVals), dtype=np.double)
        
        if self.__useFunctor:
            
            solv = self.solvFunctor.functor( np.asarray(x) )
            
        else:
            self.__model.applyOpt(x)
            
            solv = solver.ModelSolver(self.__model)
        
        if self.__verbose:
            if self.__verbose > 1:
                print("Parameters:",x)
                print("Initial Concentrations:",initVals)
                
            if self.getUseJac():
                print("Running WITH Jacobian.")
            else:
                print("Running WITHOUT Jacobian.")
        
        
        converged = False
        maxTime = self.__totalTime
        maxIter = self.__maxODEIter
        iter = 1
        
        t0 = 0
        numsteps = ceil(self.__totalTime/self.__dt)
        
        if self.__verbose > 2:
            ixprVal = True
        else:
            ixprVal = False
        
        # scpipy ODE
        if self.__integratorID == 1:
            
            if not self.__useFunctor:
                self.__rxnIdList = solv.buildCall(verbose=self.__buildVerb, useJac=self.getUseJac(),
                       transpJac=self.getUseJac(), nocheck=True, odeint=self.__odeintFormat, 
                       cythonBuild=self.__cythonBuild, noBuild=self.__noBuild)
            
            solvAlg = "lsoda"
            #solvAlg = "dopri5"
            
            if self.getUseJac():
                integrator = integrate.ode(solv,jac=solv.calcJac).set_integrator(solvAlg)
            else:
                integrator = integrate.ode(solv,jac=None).set_integrator(solvAlg)
            
            integrator.set_initial_value(initVals, 0)
            
        
        # scpipy ODEINT
        elif self.__integratorID == 2:
            
            # If we are asked to use a jacobian, we will also use the
            # transposed version for better speed. This is only an
            # option in "integrate.odeint" (not in "integrate.ode").
            if not self.__useFunctor:
                self.__rxnIdList = solv.buildCall(verbose=self.__buildVerb, useJac=self.getUseJac(), 
                           transpJac=self.getUseJac(), nocheck=True, odeint=self.__odeintFormat,
                           cythonBuild=self.__cythonBuild, noBuild=self.__noBuild)
                
                transpJac = solv.getTranspJac()
                
            else:
                transpJac = self.solvFunctor.getTranspJac()
        
        # pycvodes integrate_adaptive
        elif self.__integratorID == 3:
            
            if not self.__useFunctor:
                self.__rxnIdList = solv.buildCall(verbose=self.__buildVerb, useJac=self.getUseJac(),
                       transpJac=False, nocheck=True, odeint=self.__odeintFormat, 
                       cythonBuild=self.__cythonBuild, noBuild=self.__noBuild)
            
            def f_wrap(t, y, dydt):
                result = solv(t, np.asarray(y))
                dydt[:] = result[:]
            
            if self.getUseJac():
                def j_wrap(t, y, Jmat, dfdt=None, fy=None):
                    result = solv.calcJac(t, np.asarray(y))
                    Jmat[:,:] = result[:,:]
                    
            else:
                j_wrap = None
            
            # Absolute tolerance
            atol=1e-4
            # Relative tolerance
            rtol=1e-8
            
            ##########################################################
            
#             import ctypes
#             ctypes.CDLL("libblas.so",ctypes.RTLD_GLOBAL)
#             ctypes.CDLL("liblapack.so",ctypes.RTLD_GLOBAL)
#             
#             from assimulo.solvers import CVode
#             from assimulo.problem import Explicit_Problem
#             
#             def func_wrap(t,y):
#                 return solv(t,y)
#             
#             #Define an Assimulo problem
#             exp_mod = Explicit_Problem(func_wrap, y0=initVals, name="ODE problem.")
#             
#             if self.getUseJac():
#                 def jac_wrap(t, y):
#                     return solv.calcJac(t, y)
#                 
#                 exp_mod.jac = jac_wrap
#             
#             #Define an explicit solver
#             exp_sim = CVode(exp_mod) #Create a CVode solver
#             
#             exp_sim.verbosity = 10
#             
#             #Sets the parameters
#             exp_sim.iter  = 'Newton' #Default 'FixedPoint'
#             exp_sim.discr = 'BDF' #Default 'Adams'
#             #exp_sim.discr = 'Adams'
#             exp_sim.atol = [1e-4] #Default 1e-6
#             exp_sim.rtol = 1e-8 #Default 1e-6
#             exp_sim.sensmethod = 'SIMULTANEOUS'
#             exp_sim.maxsteps = numsteps
#             exp_sim.maxh = self.__dt
#             exp_sim.inith = self.__dt
#             
#             print(exp_sim.get_options())
            
            ##########################################################
            
            
            
        else:
            print("No integrator was set in the solver object!")
            sys.exit(5)
        
        
        while (not converged) and (iter <= maxIter):
            
            if self.__verbose > 1:
                print(iter,") T zero:",t0,"; T max:",maxTime,"; Max ODE Iter:",maxIter,
                      "; Number of steps:",numsteps)
            
            # scpipy ODE
            if self.__integratorID == 1:
                
                try:
                    while integrator.successful() and integrator.t < maxTime:
                        integrator.integrate(integrator.t+self.__dt)
                        #print(integrator.t, integrator.y)
                except Exception as inst:
                    print("Caught exception:",inst)          # __str__ allows args to be printed directly,
                    print("--> Exception type:",type(inst))    # the exception instance
                    print("--> Exception args:",inst.args)     # arguments stored in .args
                    #if self.__verbose:
                        #print("Caught exception.",inst)
                    objctv = self.__errReturn
                    print("Returning error fitness value:",objctv)
                    return objctv
                
                # If the integration was not successfull, stop and return
                if not integrator.successful():
                    if self.__verbose:
                        print("integrator unsuccessful.")
                        
                        elapsed_time = time.time() - start_time
                        print("Elapsed time:", elapsed_time,"seconds\n")
                    
                    return self.__errReturn
                
                #Fast copy
                np.copyto(finalRes,integrator.y)
            
            # scpipy ODEINT
            elif self.__integratorID == 2:
            
                tArray = np.linspace(t0, maxTime, numsteps )
                
                try:
                    if self.getUseJac():
                        results = integrate.odeint(func=solv,
                                                   Dfun=solv.calcJac,
                                                   col_deriv=transpJac,
                                                   y0=initVals, 
                                                   t=tArray, 
                                                   full_output=False, 
                                                   ixpr=ixprVal)
                    else:
                        results = integrate.odeint(func=solv,
                                                   Dfun=None,
                                                   y0=initVals, 
                                                   t=tArray, 
                                                   full_output=False, 
                                                   ixpr=ixprVal)
                except Exception as inst:
                    print("Caught exception:",inst)          # __str__ allows args to be printed directly,
                    print("--> Exception type:",type(inst))    # the exception instance
                    print("--> Exception args:",inst.args)     # arguments stored in .args
                    #if self.__verbose:
                        #print("Caught exception.",inst)
                    objctv = self.__errReturn
                    print("Returning error fitness value:",objctv)
                    return objctv
                
                #Fast copy
                np.copyto(finalRes,results[numsteps-1,:])
                
            # pycvodes integrate_adaptive
            elif self.__integratorID == 3:
                
                ###### pycvodes
                
                try:
                    xout, results, info = integrate_adaptive(f_wrap, j_wrap, initVals, 
                                     t0, maxTime, atol, rtol, nsteps=numsteps, 
                                     method=None, return_on_error=True)
                    
                except Exception as inst:
                   print("Caught exception:",inst)          # __str__ allows args to be printed directly,
                   print("--> Exception type:",type(inst))    # the exception instance
                   print("--> Exception args:",inst.args)     # arguments stored in .args
                   #if self.__verbose:
                       #print("Caught exception.",inst)
                   objctv = self.__errReturn
                   print("Returning error fitness value:",objctv)
                   return objctv
                
                
                # If the integration was not successfull, stop and return
                if not info["success"]:
                   if self.__verbose:
                       print("integrator unsuccessful.")
                       #print(info)
                       
                       elapsed_time = time.time() - start_time
                       print("Elapsed time:", elapsed_time,"seconds\n")
                   
                   return self.__errReturn
                
                
                ###### assimulo
                
#                 t1, results = exp_sim.simulate(maxTime,1)
                
                ######
                
                #Fast copy
                np.copyto(finalRes,results[-1])
            ###
            
            #if self.__verbose:
                #print("NumMet > 0:",len(np.where(finalRes > 0)[0]),"(expected",len(finalRes),")")
                #print("metabolite final concentrations:")
                #metList = [met.getID() for met in self.__model.getMetList()]
                #for i in range(len(metList)):
                    #print(metList[i],"\t",finalRes[i])
            
            # If a metabolite concentration is below zero, stop the iteration
            # and return.
            # Tests ONLY the final metabolite concentrations
            negZeroCount = len(np.where(finalRes <= 0)[0])
            if negZeroCount > 0:
                
                objctv = self.__errReturn
                
                if self.__verbose:
                    print("Error in simulation. Metabolite(s) with NEGATIVE or ZERO concentration. iter:",iter)
                    print("NumMet <= 0:",negZeroCount)
                    #print("Metabolite concentrations:")
                    #print(finalRes)
                    
                    print("Objective:", objctv)
                    elapsed_time = time.time() - start_time
                    print("Elapsed time:", elapsed_time,"seconds\n")
                
                return objctv
            
            
            if self.__odeintFormat:
                solvRes = solv(finalRes,0)
            else:
                solvRes = solv(0,finalRes)
            
            if len(np.where(abs(np.asarray(solvRes)) < self.__tol)[0]) == len(finalRes):
                converged = True
                if self.__verbose:
                    print("\tCONVERGED! iter:",iter)
            
            t0 += self.__totalTime
            maxTime += self.__totalTime
            iter += 1
            
            # Fast copy
            np.copyto(initVals, finalRes)
        
        #if self.__verbose and (not converged):
        #    print("Ran out of iterations after",iter,"loops.")
        if (not converged):
            
            objctv = self.__errReturn
            
            if self.__verbose:
                print("Ran out of iterations after",iter,"loops.")
                print("Objective:", objctv)
                elapsed_time = time.time() - start_time
                print("Elapsed time:", elapsed_time,"seconds\n")
                
            return objctv
        
        rxnFlux = solv.calcFlux(0, finalRes)
        
        rxnIdDict = dict()
        index = 0
        for rxnID in self.__rxnIdList:
            rxnIdDict[rxnID] = index
            index += 1
        
        # Normalizes all flux going through the network by the
        # glucose uptake.
        normRefFlux = 1
        if self.__normRxnID != 0:
            totalNormFlux = 0
            for rxnID in self.__normRxnID:
                if self.__verbose:
                    print("Normalizing fluxes with rxnID: \"",rxnID,
                     "\" and flux: ",rxnFlux[rxnIdDict[rxnID]])
                totalNormFlux += rxnFlux[rxnIdDict[rxnID]]
            normRefFlux = 100.0/totalNormFlux
        
        objctv = 0 
        
        if self.__verbose > 1:
            print("Calculating objective...")
            
        for trgt in self.__trgList:
            
            if self.__verbose > 1:
                print("\n",trgt.getType(), trgt.getVal(), trgt.getError())
            
            currTrg = 0
            
            if trgt.getType() == "odeRxn":
                
                # If this target should be compared against normalized
                # reaction fluxes, we apply the multiplier. Otherwise,
                # multiply by 1.
                if trgt.getNorm():
                    currNormMult = normRefFlux
                else:
                    currNormMult = 1
                
                for idMult in trgt.getIDsMultsPairs():
                    if self.__verbose > 1:
                        print(idMult[0], idMult[1], rxnFlux[rxnIdDict[idMult[0]]],
                              currNormMult)
                    currTrg += rxnFlux[rxnIdDict[idMult[0]]]*idMult[1]*currNormMult
                    
                if self.__verbose > 1:
                    if self.__targetStand == 1:
                        print("Total:",currTrg,"; Diff:",((currTrg - trgt.getVal())**2)/trgt.getError())
                    elif self.__targetStand == 2:
                        print("Total:",currTrg,"; Diff:",((currTrg - trgt.getVal())/trgt.getVal())**2)
                
            elif trgt.getType() == "odeMet":
                
                idMult = trgt.getIDsMultsPairs()
                indxs = trgt.getIndxs()
                
                for i in range(len(indxs)):
                    
                    currTrg += finalRes[indxs[0]]*idMult[1]
            
            if self.__targetStand == 1:
                objctv += trgt.getWeight()*((currTrg - trgt.getVal())**2)/trgt.getError()
            elif self.__targetStand == 2:
                objctv += trgt.getWeight()*(((currTrg - trgt.getVal())/trgt.getVal())**2)
        
        if isnan(objctv):
            objctv = 1E10
            #objctv = float("inf")
        
        if self.__verbose:
            print("Objective:", objctv)
            elapsed_time = time.time() - start_time
            print("Elapsed time:", elapsed_time,"seconds\n")
        
        return objctv
    
    
    def calcObjectiveODEINT_benchmark(self, x):
        
        start_time = time.time()
        
        start_time_applyOpt = time.time()
        self.__model.applyOpt(x)
        elapsed_time_applyOpt = time.time() - start_time_applyOpt
        
        start_time_CreateSolver = time.time()
        solv = solver.ModelSolver(self.__model)
        # If we are asked to use a jacobian, we will also use the
        # transposed version for better speed. This is only an
        # option in "integrate.odeint" (not in "integrate.ode").
        solv.buildCall(verbose=self.__buildVerb, useJac=self.getUseJac(),
                       transpJac=self.getUseJac(), nocheck=True, cythonBuild=self.__cythonBuild)
        elapsed_time_CreateSolver = time.time() - start_time_CreateSolver
        
        start_time_GetInitVal = time.time()
        initVals = np.asarray(self.__model.getInitVals())
        finalRes = np.empty(len(initVals))
        elapsed_time_GetInitVal = time.time() - start_time_GetInitVal
        
        start_time_VerbPrint1 = time.time()
        if self.__verbose:
            if self.__verbose > 1:
                print("Parameters:",x)
                print("Initial Concentrations:",initVals)
                
            if self.getUseJac():
                print("Running WITH Jacobian.")
            else:
                print("Running WITHOUT Jacobian.")
        elapsed_time_VerbPrint1 = time.time() - start_time_VerbPrint1
        
        # 10 minutes
        #self.__totalTime = 1/60
        # (1/10000000) = 0.36 ms
        #self.__dt = (1/10000000)
        #numbSteps = 10000
        #dt = totalTime/numbSteps
        
        maxTime = self.__totalTime
        maxIter = self.__maxODEIter
        iter = 1
        t0 = 0
        numsteps = ceil(self.__totalTime/self.__dt)
        
        objctv = 0
        
        if self.__verbose > 1:
            ixprVal = True
        else:
            ixprVal = False
        
        start_time_AllLoops = time.time()
        while iter <= maxIter:
            
            start_time_FullLoop = time.time()
            
            start_time_CreateTimeVec = time.time()
            tArray = np.linspace(t0, maxTime, numsteps )
            elapsed_time_CreateTimeVec = time.time() - start_time_CreateTimeVec
            
            if self.__verbose > 1:
                print(iter,") T zero:",t0,"; T max:",maxTime,"; Max ODE Iter:",maxIter,
                      "; Number of steps:",numsteps)
            
            start_time_ODEINT = time.time()
            #try:
            if self.getUseJac():
                results = integrate.odeint(func=solv,
                                           Dfun=solv.calcJac,
                                           col_deriv=solv.getTranspJac(),
                                           y0=initVals, 
                                           t=tArray, 
                                           full_output=False, 
                                           ixpr=ixprVal)
            else:
                results = integrate.odeint(func=solv,
                                           Dfun=None,
                                           y0=initVals, 
                                           t=tArray, 
                                           full_output=False, 
                                           ixpr=ixprVal)
            #except Exception as inst:
                ##print(type(inst))    # the exception instance
                ##print(inst.args)     # arguments stored in .args
                ##print(inst)          # __str__ allows args to be printed directly,
                #if self.__verbose:
                    #print("Caught exception.",inst)
                #return self.__errReturn
            
            elapsed_time_ODEINT = time.time() - start_time_ODEINT
            
            # If a metabolite concentration is below zero, stop the iteration
            # and return
            # Tests ALL metabolite concentrations throughout the simulations (SLOW)
            #for timeStep in range(len(results)):
                
                #if len(np.where(results[timeStep,:] < 0)[0]) > 0:
                    
                    #if self.__verbose:
                        #print("metabolite negative concentration.")
                    
                    #return self.__errReturn
            
            start_time_GetRes = time.time()
            #finalRes = results[numsteps-1,:]
            np.copyto(finalRes,results[numsteps-1,:])
            elapsed_time_GetRes = time.time() - start_time_GetRes
            
            #if self.__verbose:
                #print("NumMet > 0:",len(np.where(finalRes > 0)[0]),"(expected",len(finalRes),")")
                #print("metabolite final concentrations:")
                #metList = [met.getID() for met in self.__model.getMetList()]
                #for i in range(len(metList)):
                    #print(metList[i],"\t",finalRes[i])
            
            start_time_ConcTest = time.time()
            # If a metabolite concentration is below zero, stop the iteration
            # and return.
            # Tests ONLY the final metabolite concentrations
            if len(np.where(finalRes < 0)[0]) > 0:
                
                objctv = self.__errReturn
                
                if self.__verbose:
                    print("Error in simulation. Metabolite(s) with NEGATIVE concentration. iter:",iter)
                    print("NumMet < 0:",len(np.where(finalRes < 0)[0]))
                    print("Metabolite concentrations:")
                    print(finalRes)
                    
                    print("Objective:", objctv)
                    elapsed_time = time.time() - start_time
                    print("Elapsed time:", elapsed_time,"seconds\n")
                
                return objctv
            
            elapsed_time_ConcTest = time.time() - start_time_ConcTest
            
            start_time_ZeroTest = time.time()
            
            if len(np.where(finalRes > 0)[0]) != len(finalRes):
                
                objctv = self.__errReturn
                
                if self.__verbose:
                    print("Error in simulation. Metabolite(s) with ZERO concentration. iter:",iter)
                    print("NumMet > 0:",len(np.where(finalRes > 0)[0]),"(expected",len(finalRes),")")
                    print("Metabolite concentrations:")
                    print(finalRes)
                    
                    print("Objective:", objctv)
                    elapsed_time = time.time() - start_time
                    print("Elapsed time:", elapsed_time,"seconds\n")
                    
                return objctv
            
            elapsed_time_ZeroTest = time.time() - start_time_ZeroTest
            
            start_time_ConvrgTest = time.time()
            
            if len(np.where(abs(np.asarray(solv(finalRes,0))) < self.__tol)[0]) == len(finalRes):
                if self.__verbose:
                    print("\tCONVERGED! iter:",iter)
                break
            
            t0 += self.__totalTime
            maxTime += self.__totalTime
            iter += 1
            
            elapsed_time_ConvrgTest = time.time() - start_time_ConvrgTest
            
            start_time_CopyInitVal = time.time()
            np.copyto(initVals, finalRes)
            elapsed_time_CopyInitVal = time.time() - start_time_CopyInitVal
            
            elapsed_time_FullLoop = time.time() - start_time_FullLoop
            
            #print("Elapsed time FullLoop:", elapsed_time_FullLoop,"seconds")
            #print("Elapsed time CreateTimeVec:", elapsed_time_CreateTimeVec,"seconds (",elapsed_time_CreateTimeVec/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time ODEINT:", elapsed_time_ODEINT,"seconds (",elapsed_time_ODEINT/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time GetRes:", elapsed_time_GetRes,"seconds (",elapsed_time_GetRes/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time ConcTest:", elapsed_time_ConcTest,"seconds (",elapsed_time_ConcTest/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time ZeroTest:", elapsed_time_ZeroTest,"seconds (",elapsed_time_ZeroTest/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time ConvergTest:", elapsed_time_ConvrgTest,"seconds (",elapsed_time_ConvrgTest/elapsed_time_FullLoop*100,"%)")
            #print("Elapsed time CopyInitVal:", elapsed_time_CopyInitVal,"seconds (",elapsed_time_CopyInitVal/elapsed_time_FullLoop*100,"%)")
            #print("")
            
        
        elapsed_time_AllLoops = time.time() - start_time_AllLoops
        
        if self.__verbose and (iter == maxIter):
            print("Ran out of iterations.")
        
        start_time_CalcFlux = time.time()
        rxnFlux = solv.calcFlux(0, finalRes)
        elapsed_time_CalcFlux = time.time() - start_time_CalcFlux
        
        start_time_NormFlux = time.time()
        rxnIdDict = dict()
        index = 0
        for rxnID in self.__rxnIdList:
            rxnIdDict[rxnID] = index
            index += 1
            
        # Normalizes all flux going through the network by the
        # glucose uptake.
        normRefFlux = 1
        if self.__normRxnID != 0:
            totalNormFlux = 0
            for rxnID in self.__normRxnID:
                if self.__verbose:
                    print("Normalizing fluxes with rxnID: \"",rxnID,
                     "\" and flux: ",rxnFlux[rxnIdDict[rxnID]])
                totalNormFlux += rxnFlux[rxnIdDict[rxnID]]
            normRefFlux = 100.0/totalNormFlux
        elapsed_time_NormFlux = time.time() - start_time_NormFlux
        
        start_time_CalcObj = time.time()
        if self.__verbose > 1:
            print("Calculating objective...")
        
        for trgt in self.__trgList:
            
            if self.__verbose > 1:
                print("\n",trgt.getType(), trgt.getVal(), trgt.getError())
            
            currTrg = 0
            
            if trgt.getType() == "odeRxn":
                
                # If this target should be compared against normalized
                # reaction fluxes, we apply the multiplier. Otherwise,
                # multiply by 1.
                if trgt.getNorm():
                    currNormMult = normRefFlux
                else:
                    currNormMult = 1
                
                for idMult in trgt.getIDsMultsPairs():
                    if self.__verbose > 1:
                        print(idMult[0], idMult[1], rxnFlux[idMult[0]],currNormMult)
                    currTrg += rxnFlux[idMult[0]]*idMult[1]*currNormMult
                    
                if self.__verbose > 1:
                    if self.__targetStand == 1:
                        print("Total:",currTrg,"; Diff:",((currTrg - trgt.getVal())**2)/trgt.getError())
                    elif self.__targetStand == 2:
                        print("Total:",currTrg,"; Diff:",((currTrg - trgt.getVal())/trgt.getVal())**2)
                
            elif trgt.getType() == "odeMet":
                
                idMult = trgt.getIDsMultsPairs()
                indxs = trgt.getIndxs()
                
                for i in range(len(indxs)):
                    
                    currTrg += finalRes[indxs[0]]*idMult[1]
            
            if self.__targetStand == 1:
                objctv += ((currTrg - trgt.getVal())**2)/trgt.getError()
            elif self.__targetStand == 2:
                objctv += ((currTrg - trgt.getVal())/trgt.getVal())**2
        elapsed_time_CalcObj = time.time() - start_time_CalcObj
        
        start_time_CheckNan = time.time()
        if isnan(objctv):
            objctv = 1E10
            #objctv = float("inf")
        elapsed_time_CheckNan = time.time() - start_time_CheckNan
        
        if self.__verbose:
            print("Objective:", objctv)
            elapsed_time = time.time() - start_time
            print("Elapsed time:", elapsed_time,"seconds\n")
            
            print("Elapsed time ApplyOpt:", elapsed_time_applyOpt,"seconds (",elapsed_time_applyOpt/elapsed_time*100,"%)")
            print("Elapsed time CreateSolver:", elapsed_time_CreateSolver,"seconds (",elapsed_time_CreateSolver/elapsed_time*100,"%)")
            print("Elapsed time GetInitVal:", elapsed_time_GetInitVal,"seconds (",elapsed_time_GetInitVal/elapsed_time*100,"%)")
            print("Elapsed time VerbosePrint:", elapsed_time_VerbPrint1,"seconds (",elapsed_time_VerbPrint1/elapsed_time*100,"%)")
            print("Elapsed time AllLoops:", elapsed_time_AllLoops,"seconds (",elapsed_time_AllLoops/elapsed_time*100,"%)")
            print("Elapsed time CalcFlux:", elapsed_time_CalcFlux,"seconds (",elapsed_time_CalcFlux/elapsed_time*100,"%)")
            print("Elapsed time NormFlux:", elapsed_time_NormFlux,"seconds (",elapsed_time_NormFlux/elapsed_time*100,"%)")
            print("Elapsed time Calcobj:", elapsed_time_CalcObj,"seconds (",elapsed_time_CalcObj/elapsed_time*100,"%)")
            print("Elapsed time CheckNan:", elapsed_time_CheckNan,"seconds (",elapsed_time_CheckNan/elapsed_time*100,"%)")
            print(" ------------------------------------------------------ \n")
        
        return objctv

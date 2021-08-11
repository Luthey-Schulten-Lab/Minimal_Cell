"""
.. module:: desolver
   :platform: Unix
   :synopsis: The package provides an interface for ODE model solving.

.. moduleauthor:: Marcelo C. R. Melo <melomcr@gmail.com>

"""

from odecell import modelbuilder

import types, csv, collections, sys
from scipy import integrate
import numpy as np
from string import Template
from copy import deepcopy
from functools import reduce

import cobra
from pickle import load, dump

import cython
import importlib, subprocess
import pyximport, os
pyximport.install(build_dir=os.getcwd()+'/pyxbld',setup_args={'include_dirs':np.get_include()}, reload_support=True)

class ModelSolver():
    
    def __init__(self, model):
        self.model = model
        
        self.__transpJac = False
        
        self.__currParVals = 0
        
        self.localDic = dict()
        self.compiledCodeCall = 0
        self.compiledCodeJac = 0
        self.compiledCodeFlux = 0
        
    def __call__(self,t,y):
        return self.myCall(t,y)
    
    def getTranspJac(self):
        return self.__transpJac
    
    #def calcODEINT(self,y,t):
    #    return 0
    
    def calcJac(self,t,y):
        return 0
    
    #def calcJacODEINT(self,y,t):
    #    return 0
    
    def calcFlux(self,t,y):
        return 0
    
    def functor(self, parameters):
        return 0
    
    ## Optimizes the string-keys for parameters that will be optimized.
    #
    def prepareFunctor(self):
        
        optSpace, self.__currParVals = self.model.getOptSpace()
        
        indx = 0
        newVals = []
        for optItem in optSpace:
            newVals.append("self.params[" + str(indx) + "]")
            indx += 1
        
        self.model.applyOpt( newVals )
    
    def buildCall(self, verbose=0, useJac=False,
                  transpJac=False, nocheck=False, 
                  odeint=True, cythonBuild=False, 
                  functor=False, noBuild=False):
        
        #print("\nBuilding call...")
        
        self.model.prepModel()
        
        if not nocheck:
            if self.model.checkModel() != 0:
                print("--> Call function cannot be built!")
                return 1
            
            if useJac:
                if self.model.checkModelJac() != 0:
                    print("--> Jacobian function cannot be built!")
                    return 1
        
        
        self.__transpJac = transpJac
        
        ident1 = "    "
        ident = ident1
        
        if functor:
            if not cythonBuild:
                print("ERROR) Functor build is only available with cython!")
                sys.exit(2)
                
            #if not useJac:
                #print("ERROR) Functor build is only available with Jacobian definitions!")
                #sys.exit(2)
        
        cythVarPref = ""
        if cythonBuild:
            
            if verbose:
                print("Cythonizing the call!\n")
            
            cythVarPref = "cdef double "
            
            header = ""
            
            header += "import cython\n"
            
            #if useJac:
            header += "import numpy as np\n"
            header += "cimport numpy as np\n\n"
            header += "ctypedef np.double_t DTYPEDBL_t\n\n"
            
            header += "if cython.compiled:\n"
            header += "    print(\"Yep, I'm compiled.\")\n"
            header += "else:\n"
            header += "    print(\"Just a lowly interpreted script.\")\n"
            
            # Only used in case we build the solver as a functor.
            classDef = ""
            
            if functor:
                
                # Raises identation level for all default scope
                ident = ident1*2
                
                classDef = "cdef class stepClass:\n\n"
                
                classDef += ident1 + "cdef public params\n\n"
                
                classDef += ident1 + "def __init__(self, np.ndarray[DTYPEDBL_t, ndim=1] newParams):\n"
                classDef += ident1*2 + "cdef int i\n"
                classDef += ident1*2 + "cdef int n = len(newParams)\n"
                classDef += ident1*2 + "self.params = np.zeros([n], dtype=np.double)\n"
                classDef += ident1*2 + "for i in range(n):\n"
                classDef += ident1*3 + "self.params[i] = newParams[i]\n\n"
                
                calcFluxDef = ident1 + "def calcFlux(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n"
                calcFluxDef += ident1*2 + "return self.calcFlux_c(t, y) \n\n"
                
                calcFluxDef += ident1 + "cdef np.ndarray[DTYPEDBL_t, ndim=1] calcFlux_c(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
                
                
                #np.ndarray[DTYPEDBL_t, ndim=1]
                if odeint:
                    
                    callDef = ident1 + "def __call__(self, np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n"
                    callDef += ident1*2 + "return self.compiledCall(y, t)\n\n"
                    
                    # Call with parameters arranged for scipy.integrate.odeint
                    callDef += ident1 + "cdef compiledCall(self, np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n\n"
                    
                    
                    callJacDef = ident1 + "def calcJac(self, np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n"
                    callJacDef += ident1*2 + "return self.compiledCallJac(y, t)\n\n"
                    
                    # Call for Jacobian calculation in scipy.integrate.odeint
                    callJacDef += ident1 + "cdef compiledCallJac(self, np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n\n"
                    
                else:
                    callDef = ident1 + "def __call__(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n"
                    callDef += ident1*2 + "return self.compiledCall(t, y)\n\n"
                    
                    # Call with parameters arranged for scipy.integrate.ode
                    callDef += ident1 + "cdef compiledCall(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
                    
                    
                    # Call for Jacobian calculation.
                    callJacDef = ident1 + "def calcJac(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n"
                    callJacDef += ident1*2 + "return self.compiledCallJac(t, y)\n\n"
                    
                    # Call for Jacobian calculation in scipy.integrate.odeint
                    callJacDef += ident1 + "cdef compiledCallJac(self, float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
                    
                
            else:
                #np.ndarray[DTYPEDBL_t, ndim=1]
                if odeint:
                    # Call with parameters arranged for scipy.integrate.odeint
                    callDef = "cpdef compiledCall(np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n\n"
                    
                    # Call for Jacobian calculation in scipy.integrate.odeint
                    callJacDef = "cpdef compiledCallJac(np.ndarray[DTYPEDBL_t, ndim=1] y, float t):\n\n"
                else:
                    # Call with parameters arranged for scipy.integrate.ode
                    callDef = "cpdef compiledCall(float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
                    
                    # Call for Jacobian calculation.
                    callJacDef = "cpdef compiledCallJac(float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
                
                calcFluxDef = "def calcFlux(float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n"
                calcFluxDef += "    return calcFlux_c(t, y) \n\n"
                
                calcFluxDef += "cdef np.ndarray[DTYPEDBL_t, ndim=1] calcFlux_c(float t, np.ndarray[DTYPEDBL_t, ndim=1] y):\n\n"
            
        else:
        
            # Go ask scipy why ode and odeint ask for different function calls!
            if odeint:
                # Call with parameters arranged for scipy.integrate.odeint
                callDef = "def compiledCall(self, y,t):\n\n"
                
                # Call for Jacobian calculation in scipy.integrate.odeint
                callJacDef = "def compiledCallJac(self, y,t):\n\n"
            else:
                # Call with parameters arranged for scipy.integrate.ode
                callDef = "def compiledCall(self, t,y):\n\n"
                
                # Call for Jacobian calculation.
                callJacDef = "def compiledCallJac(self, t,y):\n\n"
            
            calcFluxDef = "def calcFlux(self, t,y=[]):\n\n"
        
        
        paramLines = ""
        
        # Writes all explicit parameters in the final function.
        for parName,parObj in self.model.getExplicitParameters().items():
            
            paramLines += ident + "# " + parName + "\n"
            paramLines += ident + cythVarPref + parObj.formKey + " = " + str(parObj.val) + "\n\n"
        
        paramLines += "\n\n"
        
        if not cythonBuild:
            callDef += paramLines
        callJacDef += paramLines
        calcFluxDef += paramLines
        
        metInitLines = ""
        
        # Loops over all metabolites whose concentrations are tracked
        # and, for each one, sets the metabolite name to the input vector value.
        metIndx = 0
        for met in self.model.getMetList():
            
            metID = met.getID()
            metName = met.getName()
            
            metInitLines += ident + "# (" + metID + ") " + metName + "\n"
            metInitLines += ident + cythVarPref + metID + " = y[" + str(metIndx) + "]\n\n"
            
            metIndx += 1
            
        metInitLines += "\n\n"
        
        if not cythonBuild:
            callDef += metInitLines
        calcFluxDef += metInitLines
        callJacDef += metInitLines
        
        if cythonBuild:
            
            if functor:
                callDef += ident + "cdef np.ndarray[DTYPEDBL_t, ndim=1] fluxes = self.calcFlux_c(t,y)\n\n"
            else:
                callDef += ident + "cdef np.ndarray[DTYPEDBL_t, ndim=1] fluxes = calcFlux_c(t,y)\n\n"
        
        ##########################
        ## Writes all final rate forms for model reactions
        
        definedGrads = set()
        
        metIDList = [met.getID() for met in self.model.getMetList()]
        metIDList = np.asarray(metIDList)
        numMets = len(metIDList)
        
        depMetIDSet = reduce(set.union,[met.getDependMets() for met in self.model.getMetList()])
        
        depMetMapping = dict()
        for met in self.model.getMetList():
            for depMetID in met.getDependMets():
                if depMetID in depMetMapping.keys():
                    print("ERROR! Dependent Metabolite cannot be dependent on more than one free metabolite.")
                    sys.error(2)
                depMetMapping[depMetID] = met.getID()
        
        if verbose > 1:
            print("Dependent metabolites:",depMetIDSet)
            print("Dependencies:",depMetMapping)
        
        rxnLines = ""
        rxnLinesJac = "" # So we keep only auxiliary rxns.
        
        rxnRetunVals = []
        rxnRetunValID = []
        
        prefix = "V"
        
        # Loops over all reactions, compiles their reaction rates,
        # and writes in the function.
        rxnCounter = 0
        for rxn in self.model.getRxnList() :
            
            rxnID = rxn.getID()
            rxnName = rxn.getName()
            
            rxnLines += ident + "# (" + rxnID + ") " + rxnName + "\n"
            
            rxnLines += ident + cythVarPref + rxn.getFinalRate(prefix) + "\n\n"
            
            rxnRetunVals.append(rxn.getResult())
            rxnRetunValID.append(rxnID)
            
            if rxn.getNumber() == 0:
                
                if (verbose > 1):
                    print("\n#######################################\n")
                    print(rxnID,"-->",rxn.getResult())
                    print(rxn.getFinalRate(prefix))
                
                rxnLinesJac += ident + "# (" + rxnID + ") " + rxnName + "\n"
                rxnLinesJac += ident + cythVarPref + rxn.getFinalRate(prefix) + "\n\n"
                
                if (verbose > 1):
                    print("ReactionKeys:",rxn.getKeysVals())
                    print("Dependent Rxn Keys:",rxn.getDependentKeysVals())
                
                # For each reaction, get its gradient components and the respective
                # metabolite ID (that is, the metDI for the variable with respect to
                # which the rate form was differentiated).
                for gradDat in rxn.getGrad():
                    
                    if (verbose > 1):
                        print(gradDat,np.where(metIDList == gradDat[0]))
                    
                    # In case there is a gradient component defined but its key is
                    # bound to a fixed parameter, rather than to a metabolite, we skip it.
                    # One example is the transport and modification of a molecule 
                    # from the extracellular environment (e.g, PTS complex), 
                    # the extracellular could be considered constant, but not the
                    # intracellular concentration of the molecule.
                    if len(np.where(metIDList == gradDat[0])[0]) < 1 and not (gradDat[0] in depMetIDSet):
                        continue
                    
                    # For dependent reactions, we will write its gradient with respect to a free metabolite
                    # and not with respect to the dependent metabolite it may be defined with.
                    # The rate froms must reflect this definition and have, in its gradient components,
                    # direct dependencies on the dependent metabolite derivatives.
                    if gradDat[0] in depMetIDSet:
                        targMet =  depMetMapping[gradDat[0]] 
                    else:
                        targMet = gradDat[0]
                    
                    partName = "d" + rxn.getResult() + "_d" + targMet
                    
                    if not partName in definedGrads:
                        rxnLinesJac += ident + "# d/d" + targMet+ "(" + rxn.getID() + ") [" + gradDat[0] + "] " + rxn.getName() + "\n"
                        rxnLinesJac += ident + cythVarPref + partName + " = " + gradDat[1]
                        rxnLinesJac += "\n\n"
                        definedGrads.add(partName)
                        
            elif cythonBuild:
                callDef += ident + "# (" + rxnID + ") " + rxnName + "\n"
                callDef += ident + cythVarPref + rxn.getResult() + " = fluxes[" + str(rxnCounter) + "]\n\n"
            
            rxnCounter += 1
        
        rxnLines += "\n\n"
        rxnLinesJac += "\n\n"
        
        if not cythonBuild:
            callDef += rxnLines
        calcFluxDef += rxnLines
        callJacDef += rxnLinesJac
        
        ##########################
        ## Writes all final forms for model jacobians
        
        # Prepares and pre-populates a 2D list of all IvJ results.
        jacResIJ = []
        for i in range(numMets):
            jacResIJ.append(["0"]*numMets)
        
        # Loops over all metabolites whose concentrations are tracked
        # and, for each one, determine all jacobian elements.
        for met in self.model.getMetList():
            
            metID = met.getID()
            
            if (verbose > 1):
                print("\n#######################################\n")
                print(metID)
            
            # The first "[0]" gets the first array of indices that "np.where" 
            # returns (it is also the only array of indices since it is a 1D list).
            # The second "[0]" gets the first (and only) index (since metabolite
            # IDs cannot be repeated).
            iIndx = np.where(metIDList == metID)[0][0]
            
            # Each metabolite may have non-zero jacobian elements with respect
            # to different metabolites, depending on the reaction with which it
            # is asociated.
            for rxnIndex in met.getReactions():
                
                rxn = self.model.getReaction(rxnIndex)
                
                if (verbose > 1):
                    print("\n",rxn)
                
                # Gets the stoichiometry for the current metabolite in this reacion
                stoic = rxn.getStoichiometry(metID)
                # Prepares the string for the final call function.
                # If the stoichiometry is positive (for a product), or zero,
                # adds a "plus" sign to build the final function.
                if stoic >= 0 :
                    stoiStr = "+" + str(stoic)
                else :
                    stoiStr = str(stoic)
                
                if (verbose > 1):
                    print("\nGradient Components:")
                
                # For each reaction, get its gradient components and the respective
                # metabolite ID (that is, the metDI for the variable with respect to
                # which the rate form was differentiated).
                for gradDat in rxn.getGrad():
                    
                    # In case there is a gradient component defined but its key is
                    # bound to a fixed parameter, rather than to a metabolite, we skip it.
                    # One example is the transport and modification of a molecule 
                    # from the extracellular environment (e.g, PTS complex), 
                    # the extracellular could be considered constant, but not the
                    # intracellular concentration of the molecule.
                    if len(np.where(metIDList == gradDat[0])[0]) < 1 and not (gradDat[0] in depMetIDSet):
                        continue
                    
                    if (verbose > 1):
                        #print(gradDat,np.where(metIDList == gradDat[0]))
                        print(gradDat)
                    
                    # For dependent reactions, we will write its gradient with respect to a free metabolite
                    # and not with respect to the dependent metabolite it may be defined with.
                    # The rate froms must reflect this definition and have, in its gradient components,
                    # direct dependencies on the dependent metabolite derivatives.
                    if gradDat[0] in depMetIDSet:
                        targMet =  depMetMapping[gradDat[0]] 
                        freeTarget = False
                    else:
                        targMet = gradDat[0]
                        freeTarget = True
                    
                    jIndx = np.where(metIDList == targMet)[0][0]
                    
                    partName = "dV" + str(rxn.getNumber()) + "_d" + targMet
                    
                    if not partName in definedGrads:
                        callJacDef += ident + "# (" + rxn.getID() + ") [" + gradDat[0] + "] " + rxn.getName() + "\n"
                        if freeTarget:
                            callJacDef += ident + cythVarPref + partName + " = " + gradDat[1]
                        else:
                            if (verbose > 1):
                                print("\tAutomatic Chain Rule for reaction",rxn.getID(), "target metabolite:",targMet)
                            callJacDef += ident + cythVarPref + partName + " = (" + gradDat[1] + ")*d" + gradDat[0] + "_d" + targMet
                        callJacDef += "\n\n"
                        definedGrads.add(partName)
                    
                    if jacResIJ[iIndx][jIndx] == "0":
                       jacResIJ[iIndx][jIndx] = stoiStr + "*" + partName
                    else:
                        jacResIJ[iIndx][jIndx] +=  " + " + stoiStr + "*" + partName
                    
                if (verbose > 1):
                    print("---")
                
        convFact = self.model.getUnitConvFactor()
        modelGR = self.model.getGR()
        
        for i in range(numMets):
            for j in range(numMets):
                
                if jacResIJ[i][j] != "0":
                    
                    # Adds term for unit conversion
                    if (convFact != 1):
                        jacResIJ[i][j] = "(" + jacResIJ[i][j] + ")*" + str(convFact)
                    
                    # Adds the term for dillution
                    if (i == j) and (modelGR != 0):
                        if not self.model.getMetList(i).getMode() == "odeonly":
                            jacResIJ[i][j] += " - " + str(modelGR)
                    
                elif (i == j) and (modelGR != 0):
                    if not self.model.getMetList(i).getMode() == "odeonly":
                        jacResIJ[i][j] = " - " + str(modelGR)
            
        ## Creates the return line for the calcJac function.
        returnJacList = []
        
        for i in range(numMets):
            
            line = jacResIJ[i][0]
            
            for j in range(1, numMets):
                
                line += "," + jacResIJ[i][j]
            
            line = "[" + line + "]"
            
            returnJacList.append(line)
        
        
        
        # Assembles the return statement for the call function.
        returnJacStr = "[" + ",\n ".join(returnJacList) + "]"
        
        callJacDef += "\n\n"
        
        #callJacDef += ident + "print(\"Calculating jacobian\")\n"
        
        # writes the return statement with the list of deltas.
        if self.__transpJac:
            callJacDef += ident + "return np.asmatrix(" + returnJacStr + ").T\n\n"
        else:
            callJacDef += ident + "return np.asmatrix(" + returnJacStr + ")\n\n"
        
        ##########################
        
        
        
        
        ## Creates the return line for the __call__ function.
        returnValsList = []
        
        # Loops over all metabolites whose concentrations are tracked
        # and, for each one, add together the results of relevant reaction rates.
        for met in self.model.getMetList():
            
            metID = met.getID()
            
            callDef += ident + cythVarPref + "Delta" + metID + " = "
            
            splitLineCounter = 1
            for rxnIndex in met.getReactions():
                
                # Counts terms to split lines:
                if splitLineCounter == 10:
                    callDef += " \n" + ident + "Delta" + metID + " += "
                    splitLineCounter = 1
                else:
                    splitLineCounter += 1
                
                rxn = self.model.getReaction(rxnIndex)
                
                # Gets the stoichiometry for the current metabolite in this reacion
                stoic = rxn.getStoichiometry(metID)
                # Prepares the string for the final call function.
                # If the stoichiometry is positive (for a product), or zero,
                # adds a "plus" sign to build the final function.
                if stoic >= 0 :
                    stoiStr = "+" + str(stoic)
                else :
                    stoiStr = str(stoic)
                
                callDef += stoiStr + "*" + prefix + str(rxn.getNumber()) + " "
            
            if (len(met.getReactions()) == 0 and met.getConnFlux() == 0 ):
                callDef += "0"
            else:
                callDef += "+ " + str(met.getConnFlux()) 
            
            # Adds line for unit conversion
            if (self.model.getUnitConvFactor() != 1):
                callDef += "\n" + ident + "Delta" + metID + " *= " + str(self.model.getUnitConvFactor())
            
            # Adds the term for dillution
            if (self.model.getGR() != 0):
                if not met.getMode() == "odeonly":
                    callDef += "\n" + ident + "Delta" + metID + " -= " + str(self.model.getGR()) + "*" + metID
            
            callDef += "\n"
            
            returnValsList.append("Delta" + metID)
        
        # Assembles the return statement for the call function.
        if cythonBuild:
            returnValStr = "np.asarray([" + ", ".join(returnValsList) + "])"
        else:
            returnValStr = "[" + ", ".join(returnValsList) + "]"
        
        callDef += "\n\n"
        
        # writes the return statement with the list of deltas.
        callDef += ident + "return " + returnValStr + "\n\n"
        
        ##########################
        ## Creates the return line for the calcFlux function.
        
        ## Possibility of printing values
        #calcFluxDef += ident + "rxnFlluxDict = collections.OrderedDict()\n"
        #for key,val in rxnRetunVals.items():
            #calcFluxDef += ident + "rxnFlluxDict[\"" + key + "\"] = " + val +  "\n"
        
        ## Assembles the return statement for the call function.
        #returnValFluxStr = "rxnFlluxDict"
        
        # Assembles the return statement for the call function.
        if cythonBuild:
            returnValFluxStr = "np.asarray([" + ", ".join(rxnRetunVals) + "])"
        else:
            returnValFluxStr = "[" + ", ".join(rxnRetunVals) + "]"
        
        # writes the return statement with the list of deltas.
        calcFluxDef += ident + "return " + returnValFluxStr + "\n\n"
        
        ##########################
        
        # Check the compilation of new call function.
        if (verbose > 1):
            #print("\n####\nDefinition of main function:\n####\n\n",callDef)
            #print("\n####\nDefinition of jacobian function:\n####\n\n",callJacDef)
            #print("\n####\nDefinition of reaction flux function:\n####\n\n",calcFluxDef)
            
            with open("testCompiledFunctions.py", "w") as outfile:
                outfile.write("\n####\n#Definition of main function:\n####\n\n")
                outfile.write(callDef)
                if useJac:
                    outfile.write("\n####\n#Definition of jacobian function:\n####\n\n")
                    outfile.write(callJacDef)
                outfile.write("\n####\n#Definition of reaction flux function:\n####\n\n")
                outfile.write(calcFluxDef)
        
        ## Implements new function into object
        
        if cythonBuild:
            
            with open("cythonCompiledFunctions.pyx", "w") as outfile:
                outfile.write(header)
                outfile.write("\n####\n#Definition of Functor class:\n####\n\n")
                outfile.write(classDef)
                outfile.write("\n####\n#Definition of reaction flux function:\n####\n\n")
                outfile.write(calcFluxDef)
                outfile.write("\n####\n#Definition of main function:\n####\n\n")
                outfile.write(callDef)
                if useJac:
                    outfile.write("\n####\n#Definition of jacobian function:\n####\n\n")
                    outfile.write(callJacDef)
            
            import time
            start_time = time.time()
            
            if noBuild:
                if verbose:
                    print("WARNING: noBuild =",noBuild)
            else:
                
                with open("setup_tmp.py", "w") as tmpFile:
                    tmpFile.write('from distutils.core import setup \n\
from Cython.Build import cythonize \n\
import numpy \n\
setup( ext_modules=cythonize("cythonCompiledFunctions.pyx", compiler_directives={"language_level": 3, "boundscheck": False }), include_dirs=[numpy.get_include()] ) \n\
')
                commandList = ["python3","setup_tmp.py","build_ext","--inplace"]
                print( subprocess.Popen(commandList, shell=False, stdout=subprocess.PIPE).stdout.read() )
            
            import cythonCompiledFunctions
            
            importlib.reload(cythonCompiledFunctions)
            elapsed_time = time.time() - start_time
            if verbose:
                print("Cython Build/Load time:", elapsed_time,"seconds\n")
            
            if functor:
                
                # Overloads "functor" with a class from the imported library
                self.functor = cythonCompiledFunctions.stepClass
                
                #self.myCall = cythonCompiledFunctions.stepClass
                #self.myCall = types.MethodType(cythonCompiledFunctions.compiledCall, self)
                
                #self.calcFlux = cythonCompiledFunctions.stepClass.calcFlux
                
                #self.calcJac = cythonCompiledFunctions.stepClass.compiledCallJac
                
            else:
                
                self.myCall = cythonCompiledFunctions.compiledCall
                #self.myCall = types.MethodType(cythonCompiledFunctions.compiledCall, self)
                
                self.calcFlux = cythonCompiledFunctions.calcFlux
                
                if useJac:
                    self.calcJac = cythonCompiledFunctions.compiledCallJac
        else:
            
            self.localDic.clear()
            
            # The following lines add new methods dynamically to the solver object.
            
            # Compiles the string into python code
            self.compiledCodeCall = compile(callDef, "<string>", "exec")
            # Executes the code. This DOES NOT execute the core of the function, it
            # DEFINES a new function itself, so that it can be called later.
            # The function is stored in a local dictionary that can be edited (calling
            # locals() returns a dictionary that cannot be edited form here).
            exec(self.compiledCodeCall, globals(), self.localDic)
            # Assigns the newly created function to a name of an existing function.
            # "types.MethodType" is used to cast the compiled funciton code into a
            # method of an object.
            self.myCall = types.MethodType(self.localDic['compiledCall'], self)
            
            if useJac:
                self.compiledCodeJac = compile(callJacDef, "<string>", "exec")
                exec(self.compiledCodeJac, globals(), self.localDic)
                self.calcJac = types.MethodType(self.localDic['compiledCallJac'], self)
            
            self.compiledCodeFlux = compile(calcFluxDef, "<string>", "exec")
            exec(self.compiledCodeFlux, globals(), self.localDic)
            self.calcFlux = types.MethodType(self.localDic['calcFlux'], self)
        
        return list(rxnRetunValID)


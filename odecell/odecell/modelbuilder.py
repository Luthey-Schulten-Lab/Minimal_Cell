"""
.. module:: desolver
   :platform: Unix
   :synopsis: The package provides an interface for ODE model building and connecting with an FBA model.

.. moduleauthor:: Marcelo C. R. Melo <melomcr@gmail.com>

"""

import types, csv, collections, sys
from scipy import integrate
from string import Template
from copy import deepcopy
from functools import reduce
import cobra
from pickle import load, dump

## Base class for defining rate forms.
#
# This class allows one to define new rate forms and automates the creation
# of specific strings defining a reaction law. The laws can then be combined
# to form a kinetic model.
class RateForm():
    
    def __init__(self, newbase = "$Vmax * $Sub1 /($Km + $Sub1)"):
        """ Constructor of the class.
        
        Sets the rate from for the newly created object and produces a
        set containing all keys in the rate form.
        
        Args:
            newbase (str or RateForm): If string, it must be the template string 
                defining the form of the new type of rate. If another rate
                form object, its template string will be copied to the current object.
        Returns:
            none
        
        """
        
        if isinstance(newbase, RateForm):
            newbase = newbase.getBaseRate()
        
        self.__baseRateTemplate = Template(newbase)
        
        self.__keySet = set()
        
        ## Gradient dictionary. Stores the differentiated form with respect 
        ## to its keys.
        self.__gradDict = dict()
        
        for item in Template.pattern.findall(newbase):
            self.__keySet.add(item[1] + item[2])
    
    # Checks if objects are equal.
    # __ne__ is implemented by default in Py3 as "not __eq__"
    def __eq__(self, other): 
        return self.__dict__ == other.__dict__
    
    ## Returns the template form of this type of rate.
    # 
    # @param self The object pointer.
    def getBaseRate(self):
        return self.__baseRateTemplate.template
    
    ### Replaces identifiers with user-supplied values and returns a string.
    #
    # This function uses a dictionary of parameter values and substrate names
    # to substitute the identifiers in the template reaction form, returning 
    # a string with the rate for a specific reaction.
    # @param self The object pointer.
    # @param subsDict The dictionary with parameter values and metabolite names.
    def getRate(self, subsDict):
        return self.__baseRateTemplate.substitute(subsDict)
    
    ## Returns the set of keys for the rate form.
    # 
    # @param self The object pointer.
    def getKeys(self):
        return self.__keySet
    
    def setGradComp(self, key, form):
        self.__gradDict[Template(key)] = Template(form)
    
    def getGradDict(self, subsDict):
        
        retDict = dict()
        
        for key,val in self.__gradDict.items():
            retDict[key.substitute(subsDict)] = val.substitute(subsDict)
        
        return retDict
    
    ## Another solution could be to have Closures. Substitute 
    ## the parameters for the values of a particular reaction, and leave
    ## the concentrations as function parameters. This function would then
    ## return a pair: the function and a list of indexes for the "y" list
    ## with concentrations that comes from the solver.
    ## The __call__ would then iterate over a a list of lambdas, calling 
    ## them with their respective metabolite concentrations (determined 
    ## using their index-vectors), and storing all results
    ## in a list, wich would be used in another loop to determine the deltas
    ## for each metabolite.
    ## 
    ##    def newFunc(par1,par2):
    ##        def finalFunc(x,y): 
    ##            return par1*x/(par1 + y/par2)
    ##        return finalFunc
    ##    
    ##    func(*shortConcentrationsList)
    

class Metabolite():
    
    ## Constructor for class Metabolite
    #
    # @param self The object pointer.
    # @param metID The ID string for this metabolite.
    # @param fbaMetID The ID used in an FBA model connected to the ODE model.
    # @param metName The name string for this metabolite.
    # @param initVal The value of the initial condition of the metabolite.
    def __init__(self, metID, metName = "", initVal = 0, fbaMetID = "", metMode=""):
        self.__ID = metID
        self.__FBAID = fbaMetID
        self.__name = metName
        self.__initVal = initVal
        self.__currVal = initVal
        self.__mode = metMode
        
        # Set of reactions where the metabolite is used.
        self.__rxnSet = set()
        
        self.__connRxns = []
        
        self.__connFlux = 0
        
        self.__dependentMets = set()
    
    def __str__(self):
        returnStr = self.__ID + ": " + self.__name + \
            ", Concentration: " + str(self.__currVal) + "\n"
        if len(self.__rxnSet) > 0:
            returnStr += "Metabolite associated with reaction(s) " + \
                str(self.__rxnSet)
        else :
            returnStr += "Metabolite not associated with any reactions."
        return returnStr
    
    # Metabolite Name
    def getName(self):
        return self.__name
    
    def getID(self):
        return self.__ID
    
    def getFBAID(self):
        return self.__FBAID
    
    def setMode(self, newMode):
        self.__mode = newMode
    
    def getMode(self):
        return self.__mode
    
    def addReaction(self, rxnIndx):
        self.__rxnSet.add(rxnIndx)
    
    def rmReaction(self, rxnIndx):
        self.__rxnSet.remove(rxnIndx)
    
    def getReactions(self):
        return self.__rxnSet
    
    def addDependMet(self, met):
        return self.__dependentMets.add(met)
    
    def getDependMets(self):
        return self.__dependentMets
    
    def getInitVal(self):
        return self.__initVal
    
    def getCurrVal(self):
        return self.__currVal
    
    def setInitVal(self, newVal):
        self.__initVal = newVal
    
    def setCurrVal(self, newVal):
        self.__currVal = newVal
    
    def setConnFlux(self, newFlux):
        self.__connFlux = newFlux
    
    def getConnFlux(self):
        return self.__connFlux
    
    def cleanConnections(self):
        self.__connRxns = []
    
    ## Add reaction index and the stoichiometry for this metabolite
    #
    # @param self The object pointer.
    # @param rxnIndx The index for the connecting reaction,
    # @param rxnStoich The stoichiometry for this metabolite in the reaction.
    def addFBAConnection(self, rxnIndx, rxnStoich):
        self.__connRxns.append( (rxnIndx, rxnStoich) )
        return(0)
    
    def getConnRxns(self):
        return self.__connRxns
    
    def calcConnFlux(self, fbsSolution, fbamodel=0):
        
        self.__connFlux = 0
        
        if fbamodel != 0:
            print("Connecting fluxes for metabolite:",self.__ID)
        
        for connRxn in self.__connRxns:
            if fbamodel != 0 and fbsSolution.x[connRxn[0]] != 0:
                print(connRxn[0], fbamodel.reactions[connRxn[0]].id, fbsSolution.x[connRxn[0]], connRxn[1])
            self.__connFlux += connRxn[1]*fbsSolution.x[connRxn[0]]
        
        return(0)
    
class Reaction():
    
    ## Constructor for class Reaction
    #
    # @param self The object pointer.
    # @param rxnID The ID string for this reaction.
    # @param rxnName The name string for this reaction.
    def __init__(self,rxnID, rxnName):
        self.__ID = rxnID
        self.__name = rxnName
        
        # formKey VS IDs
        self.__substrates = {}
        # formKey VS IDs
        self.__products = {}
        # formKey VS value
        self.__parameters = {}
        # Reaction rate form
        self.__rateForm = 0
        # Rate form name
        self.__rateFormName = ""
        # Result to be placed at the RHS of the final equation.
        self.__result = 0
        # Reaction number used to construct a "__result" when one is not 
        # supplied by the user during model construction.
        self.__number = 0
        # Equivalent reaction(s) in the FBA model.
        # The ODE flux will be imposed on these reactions
        # in a hybrid ODE/FBA simulation.
        self.__fbaEquivRxnsSet = set()
        
        self.__stoich = dict()
        
        self.__dependentRxnsIndxs = set()
        self.__depKeysVals = dict()
        
        # Indicates wether this reaction should be checked for presence of
        # Substrate(s) and Product(s).
        self.__checkReaction = True
        
    def __str__(self):
        returnStr = self.__ID + ": " + self.__name + "\n"
        returnStr += "Reaction Index: " + str(self.__number) + "\n"
        returnStr += "Rate form name: " + self.__rateFormName + "\n"
        returnStr += "Rate form: " + self.getBaseRateForm() + "\n"
        
        if len(self.__substrates) > 0:
            returnStr += "Substrates: \n" # + str(self.__substrates) + "\n"
            for key,val in self.__substrates.items():
                returnStr += "\t"+ key + ": " + val + " (" + str(self.getStoichiometry(val)) + ")\n"
        else:
            returnStr += "WARNING: No substrates have been defined!\n"
        
        if len(self.__products) > 0:
            returnStr += "Products: \n" # + str(self.__products) + "\n"
            for key,val in self.__products.items():
                returnStr += "\t"+ key + ": " + val + " (" + str(self.getStoichiometry(val)) + ")\n"
        else:
            returnStr += "WARNING: No products have been defined!\n"
        
        if len(self.__parameters) > 0:
            returnStr += "Parameters: " + str(self.__parameters) + "\n"
        
        diffSet = self.getUnboundKeys()
        
        if len(self.getUnboundKeys()) > 0:
            returnStr += "Unbound keys: " + str(diffSet)
        else:
            returnStr += "All keys are bound."
        return returnStr
    
    def getID(self):
        return self.__ID
    
    def getName(self):
        return self.__name
    
    def setRateForm(self, newRateForm, newRrateFormName):
        self.__rateForm = newRateForm
        self.__rateFormName = newRrateFormName
    #    rate = rateForms[self.rateType].getRate(self.parameters)
    
    ## Returns the set of keys for the rate form.
    # 
    # @param self The object pointer.
    def getKeys(self):
        return self.__rateForm.getKeys()
    
    ## Returns the set of unbound keys in the reaction object.
    # 
    # @param self The object pointer.
    def getUnboundKeys(self):
        
        definedKeys = set(self.__substrates.keys())
        definedKeys.update(set(self.__products.keys()))
        definedKeys.update(set(self.__parameters.keys()))
        diffSet = set(self.getKeys()).difference(definedKeys)
        
        return diffSet
    
    ## Returns the rate form name.
    # 
    # @param self The object pointer.
    def getRateFormName(self):
        return self.__rateFormName
    
    ## Returns the base rate form.
    #
    # @param self The object pointer.
    def getBaseRateForm(self):
        return self.__rateForm.getBaseRate()
    
    def getSubstrates(self):
        return self.__substrates.values()
    
    def getProducts(self):
        return self.__products.values()
    
    def setResult(self, newRes):
        """ Sets a fixed known variable name for the result (LHS) of the reaction.
        
        This function will set a fixed variable name for the result of the
        reaction, that is, to the symbol at the LHS of the equation, so that 
        this result can be used as a parameter or metabolite in other reactions.
        
        Args:
            newRes (str): Fixed result for the reaction.
        Returns:
            none
        """
        
        self.__result = newRes
        self.__number = 0
    
    def getResult(self):
        return self.__result
    
    def getStoichiometry(self, metID):
        
        stoic = 0
        
        if metID in self.__stoich.keys():
            stoic = self.__stoich[metID]
        else:
            subStoic = list(self.__substrates.values()).count(metID)
            prodStoic = list(self.__products.values()).count(metID)
            stoic = prodStoic - subStoic
            
        return stoic
    
    def getKeysVals(self):
        finalDict = deepcopy(self.__substrates)
        finalDict.update(self.__products)
        finalDict.update(self.__parameters)
        
        return finalDict
    
    def addDepKeysVals(self, depKeysValsDict):
        self.__depKeysVals.update(depKeysValsDict)
    
    def getDependentKeysVals(self):
        
        return self.__depKeysVals
    
    def getRHS(self):
        """ Returns a strign with the RHS of the reaction equation.
        
        This function uses the substrates, porducts and parameters previously 
        defined to build a dictionary of rate form keys. This dictonary is then
        passed to the RateForm object associated with this reaction object, and
        a final rate form is returned.
        
        Args:
            None
        Returns:
            Right Hand Side of the final rate form, with parameter values 
            and metabolite names as defined by the user.
         
        """
        
        return self.__rateForm.getRate(self.getKeysVals())
        
    def getFinalRate(self, prefix = "V"):
        """ Returns a strign with the final rate law for the reaction.
        
        For standardization of the LHS of model equations, the LHS is automatically
        created by concatenating a prefix (by default "V")
        to the index of the current reaction object, which will depend on the 
        order with which reactionswere added to the model.
        
        Args:
            prefix (str): The string prefix to be added to the LHS of the equation.
        Returns:
            Final rate form, with parameter values and metabolite names.
            
            The 18th reaction defined in the model would have (by default),
            the following format:
            
            V18 = RHS
            
            Where the RHS is defined by substituting the keys in the RateForm 
            associated with this reaction. One can get the RHS by calling the 
            function getRHS().
        """
        
        if self.__result == 0:
            finalResult = str(prefix) + str(self.__number)
            self.__result = finalResult
        else:
            finalResult = self.__result
        
        return finalResult + " = " + self.getRHS()
    
    def getGrad(self):
        """ Returns the diferentiated forms of this reaction.
        
        Taking dy/dj, this function determines which differentiated are available
        and returns a list of the form
        [(dy/dj1), ..., (dy/djn)]
        
        Args:
            none
        Returns:
            List of differentiated forms with respect to the dependent metabolites.
        """
        
        unionDic = deepcopy(self.getKeysVals())
        unionDic.update(self.getDependentKeysVals())
        
        return [ (key,val) for key,val in self.__rateForm.getGradDict(unionDic).items() ]
    
    def setNumber(self, newNumber):
        self.__number = newNumber
    
    def getNumber(self):
        return self.__number
    
    def setCheckRxn(self, newCheckRxn):
        self.__checkReaction = newCheckRxn
    
    def getCheckRxn(self):
        return self.__checkReaction
    
    def isSubstrate(self, metID):
        return metID in self.__substrates.values()
    
    def isProduct(self, metID):
        return metID in self.__products.values()
    
    def addSubstrate(self, rxnFormKey, metID, stoich=0):
        """ Add a substrate to the reaction.
        
        This fucntion associates a metabolite with a reaction form key.
        Defining a new metabolite to a reaction form key will overwrite any current 
        metabolite association, and will print a warning.
        If a stoichiometry is defined, its value will be used as-is when creating
        the ODE model.
        
        Args:
            rxnFormKey (str): Reaction form key.
            metID (str): Metabolite ID.
            stoich (int): Reactant stoichiometry.
        Returns:
            none
        """
        
        if rxnFormKey in self.__substrates.keys():
            if self.__substrates[rxnFormKey] == metID:
                print(metID + " is already bound to key " + rxnFormKey)
                return 2
            
            print("WARNING: reaction form key \"" + rxnFormKey + \
                "\" was previously assigned to metabolite \"" + \
                    self.__substrates[rxnFormKey] + "\".")
            print("Removing previous assignment and binding key to metabolite \"" \
                + metID + "\".")
            
            self.rmSubstrate(rxnFormKey)
        
        self.__substrates[rxnFormKey] = metID
        
        if stoich != 0:
            self.__stoich[metID] = stoich
    
    def addProduct(self, rxnFormKey, metID, stoich=0):
        
        if rxnFormKey in self.__products.keys():
            if self.__products[rxnFormKey] == metID:
                print(metID + " is already bound to key " + rxnFormKey)
                return 2
                
            print("WARNING: reaction form key \"" + rxnFormKey + \
                "\" was previously assigned to metabolite \"" + \
                    self.__products[rxnFormKey] + "\".")
            print("Removing previous assignment and binding key to metabolite \"" \
                + metID + "\".")
            
            self.rmProduct(rxnFormKey)
        
        self.__products[rxnFormKey] = metID
        
        if stoich != 0:
            self.__stoich[metID] = stoich
    
    def addParameter(self, rxnFormKey, value, verbose = 0):
        
        if rxnFormKey in self.__parameters.keys():
            if self.__parameters[rxnFormKey] == value:
                if verbose > 1:
                    print("Key " + rxnFormKey + " was already assigned the value "\
                        + str(value) + " in reaction " + self.__ID)
                return 2
                
            if verbose > 0:
                print("WARNING: reaction form key \"" + rxnFormKey + \
                    "\" was previously assigned to value \"" + \
                        str(self.__parameters[rxnFormKey]) + "\".")
                print("Overwriting assignment to value \"" + str(value) + "\".")
        
        self.__parameters[rxnFormKey] = value
    
    ## Adds an equivalent FBA reaction to the ODE reaction
    # 
    # The multiplier can be used to change the reaction's
    # direction, in order to match the default direction
    # in the FBA model.
    # 
    def addFBAEquiv(self, fbaEquiv, multiplier = 1):
        
        if isinstance(fbaEquiv,str):
            self.__fbaEquivRxnsSet.add( (fbaEquiv, multiplier) )
            return(0)
            
        elif isinstance(fbaEquiv,set):
            self.__fbaEquivRxnsSet.update(fbaEquiv)
            
            return(0)
        else :
            return(1)
    
    def getFBAEquiv(self):
        return self.__fbaEquivRxnsSet
    
    def cleanFBAEquiv(self):
        self.__fbaEquivRxnsSet = set()
    
    def rmSubstrate(self, rxnFormKey):
        
        if rxnFormKey not in self.__substrates.keys():
            print("ERROR: Key \"" + rxnFormKey + \
                "\" is not associated with a substrate.")
            return 1
        
        returnMetID = self.__substrates[rxnFormKey]
        
        del self.__substrates[rxnFormKey]
        
        return returnMetID
    
    def rmProduct(self, rxnFormKey):
        
        if rxnFormKey not in self.__products.keys():
            print("ERROR: Key \"" + rxnFormKey + \
                "\" is not associated with a product.")
            return 1
        
        returnMetID = self.__products[rxnFormKey]
        
        del self.__products[rxnFormKey]
        
        return returnMetID
    
    def rmParameter(self, rxnFormKey):
        
        if rxnFormKey not in self.__parameters.keys():
            print("ERROR: Key \"" + rxnFormKey + \
                "\" is not associated with a parameter.")
            return 1
        
        del self.__parameters[rxnFormKey]
        
        return 0
    
    def addDependent(self, rxnIndx):
        self.__dependentRxnsIndxs.add(rxnIndx)
    
    def getDependentRxns(self):
        return self.__dependentRxnsIndxs
    
    def hasMetabolite(self, metID):
        return (metID in self.__substrates.values()) or (metID in self.__products.values())
    

class Parameter():
    
    def __init__(self, newName, newVal, newUnit, newFormKey):
        self.name = newName
        self.val = newVal
        self.unit = newUnit
        self.formKey = newFormKey
        self.lb = 0
        self.ub = 0
    
    def __init__(self, newName, newVal, newUnit, newFormKey, parLB, parUB):
        self.name = newName
        self.val = newVal
        self.unit = newUnit
        self.formKey = newFormKey
        self.lb = parLB
        self.ub = parUB

class OptDat():
    
    def __init__(self):
        self.type = ""
        self.indx = -1
        self.field = ""
        self.val = 0
        self.lb = 0
        self.ub = 0
    
    def __init__(self, parDat):
        self.type = ""
        self.indx = -1
        self.field = parDat.formKey
        self.val = parDat.val
        self.lb = parDat.lb
        self.ub = parDat.ub
    

class MetabolicModel():
    
    
    def __init__(self):
        
        ## Relates IDs to indexes
        self.__metDict = {}
        ## Ordered list of metabolite objects
        self.__metList = []
        
        ## Relates rxn IDs to indexes
        self.__rxnDict = {}
        ## Ordered list of reaction objects
        self.__rxnList = []
        
        self.numMainRxns = 0
        
        ## Growth rate of cell
        self.__gr = 0
        
        ## Unit conversion for metabolite concentrations to be used 
        ## with FBA models.
        self.__unitConvFactor = 1
        
        self.__explicitParam = collections.OrderedDict()
        self.__globalParam = collections.OrderedDict()
        self.__rxnParam = collections.OrderedDict()
        
        self.__formDict = collections.OrderedDict()
        
        self.__hybrid_ODE_FBA = False
        self.__cobraModel = 0
        self.__fbaEquivRxns = set()
        self.__cobraModelSolution = 0
        
        self.__optSpaceList = []
        
        self.zeroOrder = RateForm("$K")
        
        self.firstOrder = RateForm("$K*$Sub1")
        self.firstOrder.setGradComp("$Sub1","$K")
        
        self.secondOrder = RateForm("$K*$Sub1*$Sub2")
        self.secondOrder.setGradComp("$Sub1","$K*$Sub2")
        self.secondOrder.setGradComp("$Sub2","$K*$Sub1")
        
        self.thirdOrder = RateForm("$K*$Sub1*$Sub2*$Sub3")
        self.thirdOrder.setGradComp("$Sub1","$K*$Sub2*$Sub3")
        self.thirdOrder.setGradComp("$Sub2","$K*$Sub1*$Sub3")
        self.thirdOrder.setGradComp("$Sub3","$K*$Sub1*$Sub2")
        
        self.massActionRev = RateForm("$K1*($Sub1 - $Prod1/$Keq)")
        self.massActionRev.setGradComp("$Sub1","$K1")
        self.massActionRev.setGradComp("$Prod1","-$K1/$Keq")
        
        self.firstOrderMM = RateForm("$Vmax*$Sub1/($Km + $Sub1)")
        self.firstOrderMM.setGradComp("$Sub1","$Vmax/($Km + $Sub1) - $Vmax/(($Km + $Sub1)**2)")
    
        self.secondOrderMM = RateForm("$Vmax*$Sub1*$Sub2/($Ki1*$Km2 + $Km2*$Sub1\
+ $Km1*$Sub2 + $Sub1*$Sub2)")
        self.secondOrderMM.setGradComp("$Sub1", "$Vmax*$Sub2/($Ki1*$Km2 + $Km2*$Sub1\
+ $Km1*$Sub2 + $Sub1*$Sub2) - ($Vmax*$Sub1*$Sub2/(($Ki1*$Km2 + $Km2*$Sub1\
+ $Km1*$Sub2 + $Sub1*$Sub2)**2))*($Km2 + $Sub2)")
        self.secondOrderMM.setGradComp("$Sub2", "$Vmax*$Sub1/($Ki1*$Km2 + $Km2*$Sub1\
+ $Km1*$Sub2 + $Sub1*$Sub2) - ($Vmax*$Sub1*$Sub2/(($Ki1*$Km2 + $Km2*$Sub1\
+ $Km1*$Sub2 + $Sub1*$Sub2)**2))*($Km1 + $Sub1)")
        
        self.QpHdep = RateForm("(1 + 2*(10**($pHm1 -$pK1)))/(1 +\
10**($pH - $pK1) + 10**(2*$pHm1 -$pH -$pK1))")
        
        self.ZpHdep = RateForm("1+ $Proton/$Kh1 + $Kh2/$Proton")
        
        self.firstOrderMMRev = RateForm("$Vmax*$pHdep*($Sub1 - $Prod1/$Keq) / \
($Km1 + $Sub1 + $Km1*$Prod1/$Kmp1)")
        self.firstOrderMMRev.setGradComp("$Sub1","$Vmax*$pHdep/($Km1 + $Sub1 + \
$Km1*$Prod1/$Kmp1) + $Vmax*$pHdep*($Prod1/$Keq - $Sub1)/(($Km1 + $Sub1 + $Km1*$Prod1/$Kmp1)**2)")
        self.firstOrderMMRev.setGradComp("$Prod1","(-$Vmax*$pHdep/$Keq)/($Km1 + $Sub1 + \
$Km1*$Prod1/$Kmp1) + $Vmax*$pHdep*($Km1/$Kmp1)*($Prod1/$Keq - $Sub1)/(($Km1 + \
$Sub1 + $Km1*$Prod1/$Kmp1)**2)")
        
        self.firstOrderMMRevInh = RateForm("$Vmax*$pHdep*($Sub1 - $Prod1/$Keq) /\
($Km1 + $Sub1 + $Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1)")
        self.firstOrderMMRevInh.setGradComp("$Sub1","$Vmax*$pHdep/($Km1 + $Sub1 + \
$Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1) + $Vmax*$pHdep*($Prod1/$Keq - $Sub1)/(($Km1\
 + $Sub1 + $Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1)**2)")
        self.firstOrderMMRevInh.setGradComp("$Prod1","(-$Vmax*$pHdep/$Keq)/($Km1 + $Sub1 + \
$Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1) + $Vmax*$pHdep*($Km1/$Kmp1)*($Prod1/$Keq\
 - $Sub1)/(($Km1 + $Sub1 + $Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1)**2)")
        self.firstOrderMMRevInh.setGradComp("$Inh1","-$Vmax*$pHdep*($Sub1 - \
$Prod1/$Keq)*($Km1/$KefInh1)/(($Km1 + $Sub1 + $Km1*$Prod1/$Kmp1 + $Km1*$Inh1/$KefInh1)**2)")
        
        self.orderedUniBiRev = RateForm("$Vmax*$pHdep*($Sub1 - $Prod1*$Prod2/$Keq) / \
($Km1 + $Sub1 + $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2)) )")
        self.orderedUniBiRev.setGradComp("$Sub1","$Vmax*$pHdep/($Km1 + $Sub1\
+ $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2)) ) - $Vmax*$pHdep*($Sub1 \
- $Prod1*$Prod2/$Keq)/(($Km1 + $Sub1 + $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)) )**2)")
        self.orderedUniBiRev.setGradComp("$Prod1","(-$Vmax*$pHdep*$Prod2/$Keq)/($Km1 + $Sub1\
+ $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2)) ) - $Vmax*$pHdep*($Sub1 \
- $Prod1*$Prod2/$Keq)*($Km1*(1/$Kmp1 + $Prod2/($Kmp1*$Kmp2)))/(($Km1 + $Sub1 + $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)) )**2)")
        self.orderedUniBiRev.setGradComp("$Prod2","(-$Vmax*$pHdep*$Prod1/$Keq)/($Km1 + $Sub1\
+ $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2)) ) - $Vmax*$pHdep*($Sub1 \
- $Prod1*$Prod2/$Keq)*($Km1*(1/$Kmp2 + $Prod1/($Kmp1*$Kmp2)))/(($Km1 + $Sub1 + $Km1*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)) )**2)")
        
        self.randomBiBiRev = RateForm("$Vmax*$pHdep*($Sub1*$Sub2 - $Prod1*$Prod2/$Keq) \
/ ($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )")
        self.randomBiBiRev.setGradComp("$Sub1","$Vmax*$pHdep*$Sub2 \
/ ($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*($Km2 + $Sub2)/(($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.randomBiBiRev.setGradComp("$Sub2","$Vmax*$pHdep*$Sub1 \
/ ($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*($Km1 + $Sub1)/(($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.randomBiBiRev.setGradComp("$Prod1","($Vmax*$pHdep*$Prod2/$Keq) \
/ ($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*($Kd1*$Km2*(1/$Kmp1 + $Prod2/($Kmp1*$Kmp2)))/(($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.randomBiBiRev.setGradComp("$Prod2","($Vmax*$pHdep*$Prod1/$Keq) \
/ ($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*($Kd1*$Km2*(1/$Kmp2 + $Prod1/($Kmp1*$Kmp2)))/(($Kd1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Kd1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        
        self.orderedBiBiIrev = RateForm("$Vmax*$Sub1*$Sub2 / \
($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2)")
        self.orderedBiBiIrev.setGradComp("$Sub1","$Vmax*$Sub2 / \
($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2) - $Vmax*$Sub1*$Sub2*($Km2 + $Sub2) / \
(($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2)**2)")
        self.orderedBiBiIrev.setGradComp("$Sub2","$Vmax*$Sub1 / \
($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2) - $Vmax*$Sub1*$Sub2*($Km1 + $Sub1) / \
(($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2)**2)")
        
        self.orderedBiBiRev = RateForm("$Vmax*$pHdep*($Sub1*$Sub2 - $Prod1*$Prod2/$Keq) \
/ ($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )")
        self.orderedBiBiRev.setGradComp("$Sub1","$Vmax*$pHdep*$Sub2 \
/ ($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - $Prod1*$Prod2/$Keq)*($Km2 + $Sub2) \
/ (($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.orderedBiBiRev.setGradComp("$Sub2","$Vmax*$pHdep*$Sub1 \
/ ($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - $Prod1*$Prod2/$Keq)*($Km1 + $Sub1) \
/ (($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.orderedBiBiRev.setGradComp("$Prod1","($Vmax*$pHdep*$Prod2/$Keq) \
/ ($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*( $Ki1*$Km2*(1/$Kmp1 + $Prod2/($Kmp1*$Kmp2)) ) \
/ (($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        self.orderedBiBiRev.setGradComp("$Prod2","($Vmax*$pHdep*$Prod1/$Keq) \
/ ($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) ) - $Vmax*$pHdep*($Sub1*$Sub2 - \
$Prod1*$Prod2/$Keq)*( $Ki1*$Km2*(1/$Kmp2 + $Prod1/($Kmp1*$Kmp2)) ) \
/ (($Ki1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2 + $Ki1*$Km2*($Prod1/$Kmp1 \
+ $Prod2/$Kmp2 + ($Prod1*$Prod2)/($Kmp1*$Kmp2) ) )**2)")
        
        self.randomTerBiRev = RateForm("$Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq) / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))")
        self.randomTerBiRev.setGradComp("$Sub1","$Vmax*$pHdep*$Sub2*$Sub3 \
 / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2))) - $Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq)*($Km2*$Km3 + $Km3*$Sub2 + $Km2*$Sub3 + $Sub2*$Sub3) / \
(($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))**2)")
        self.randomTerBiRev.setGradComp("$Sub2","$Vmax*$pHdep*$Sub1*$Sub3 \
 / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2))) - $Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq)*($Km1*$Km3 + $Km3*$Sub1 + $Km1*$Sub3 + $Sub1*$Sub3) / \
(($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))**2)")
        self.randomTerBiRev.setGradComp("$Sub3","$Vmax*$pHdep*$Sub1*$Sub2 \
 / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2))) - $Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq)*($Km1*$Km2 + $Km2*$Sub1 + $Km1*$Sub2 + $Sub1*$Sub2) / \
(($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))**2)")
        self.randomTerBiRev.setGradComp("$Prod1","$Vmax*$pHdep*($Prod2/$Keq) \
 / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2))) - $Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq)*($Km1*$Km2*$Km3*(1/$Kmp1 + $Prod2/($Kmp1*$Kmp2))) \
/ (($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))**2)")
        self.randomTerBiRev.setGradComp("$Prod2","$Vmax*$pHdep*($Prod1/$Keq) \
 / ($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2))) - $Vmax*$pHdep*($Sub1*$Sub2*$Sub3 - \
$Prod1*$Prod2/$Keq)*($Km1*$Km2*$Km3*(1/$Kmp2 + $Prod1/($Kmp1*$Kmp2))) \
/ (($Km1*$Km2*$Km3 + $Km2*$Km3*$Sub1 + $Km1*$Km3*$Sub2 + \
$Km1*$Km2*$Sub3 + $Km3*$Sub1*$Sub2 + $Km2*$Sub1*$Sub3 + $Km1*$Sub2*$Sub3 + \
$Sub1*$Sub2*$Sub3 + $Km1*$Km2*$Km3*($Prod1/$Kmp1 + $Prod2/$Kmp2 + \
($Prod1*$Prod2)/($Kmp1*$Kmp2)))**2)")
        
        
        self.genMWC = RateForm("$n*$Fr*(1 + ($Ft/$Fr)*$Q)/(1 + $Q)")
        
        self.updateAvailableForms()
        
        self.__verbose = 1
    
    ## Sets the verbosity level for functions that build a model.
    # 
    # All functions that add/remove reactions, metabolites and parameters
    # can send messages about its status, warnings for possibly unwanted 
    # or accidental actions, and errors that will make it impossible for the
    # model to be simulated.
    # verbosity = 0 (super quiet) only shows errors.
    # verbosity = 1 (default) shows status and warnings, plus what verbosity 0 shows.
    # verbosity = 2 (very talkative) show all that level 1 show, plus extended status messages.
    # @param self The object pointer
    # @param verbosity The new verbosity level.
    def setVerbosity(self, verbosity):
        self.__verbose = verbosity
    
    ## Returns the verbosity level.
    def getVerbosity(self):
        return self.__verbose
    
    ## Returns the verbosity level.
    def getNumMainRxns(self):
        return self.numMainRxns
    
    ## Returns the list of metabolite objects
    #
    # The function may return the entire list, if no index is not given, or 
    # a single object if a valid index is given.
    # @param self The object pointer.
    # @param index The index to for a given reaction. Default returns the full list
    def getMetList(self, index = -1):
        if (index < 0):
            return self.__metList
        else :
            if (index < len(self.__metList)):
                return self.__metList[index]
            else:
                print("ERROR: Index out of bounds!")
    
    ## Returns the dictionary of metabolites and indexes
    #
    # The function returns the entire dictionary relating metabolite IDs with
    # their respective index in the model's metabolite list.
    # @param self The object pointer.
    def getMetDict(self):
        return self.__metDict
    
    ## Returns the list of reaction objects
    #
    # The function may return the entire list, if no index is not given, or 
    # a single object if a valid index is given.
    # @param self The object pointer.
    # @param index The index to for a given reaction. Default returns the full list
    def getRxnList(self, index = -1):
        if (index < 0):
            return self.__rxnList
        else :
            if (index < len(self.__rxnList)):
                return self.__rxnList[index]
            else:
                print("ERROR: Index out of bounds!")
    
    ## Returns the dictionary of reaction and indexes
    #
    # The function returns the entire dictionary relating metabolite IDs with
    # their respective index in the model's metabolite list.
    # @param self The object pointer.
    def getRxnDict(self):
        return self.__rxnDict
    
    ## Returns a reaction object, given an identifier.
    #
    # The identifier may be a string containing the ID or an int containing the 
    # index.
    # @param self The object pointer.
    # @param rxnIdentifier The identifier, either a string or an int.
    def getReaction(self, rxnIdentifier):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        if(rxnID == ""):
            print("--> Could not find reaction.")
            return 1
        else:
            return self.__rxnList[rxnIndex]
        
    ## Returns a metabolite object, given an identifier.
    #
    # The identifier may be a string containing the ID or an int containing the 
    # index.
    # @param self The object pointer.
    # @param rxnIdentifier The identifier, either a string or an int.
    def getMetabolite(self, metIdentifier):
        
        metID, metIndex = self.parseMetIdentifier(metIdentifier)
        
        if(metID == ""):
            print("--> Could not find metabolite.")
            return 1
        else:
            return self.__metList[metIndex]
    
    ## Returns the list of available reaction forms names.
    #
    # The function updates the available reaction forms before returning the list.
    # @param self The object pointer.
    def getAvailableForms(self, fullDict=False):
        
        self.updateAvailableForms()
        
        if fullDict:
            return self.__formDict
        else:
            return(self.__formDict.keys())
    
    ## Prints the list of available reaction forms names.
    #
    # The function updates the available reaction forms before printing the list.
    # @param self The object pointer.
    def printAvailableForms(self):
        
        print(", ".join(self.getAvailableForms()))
    
    ## Updates the reaction form defined in the current model object.
    #
    # The function parses the object looking for attributes of the class RateForm,
    # and creates a dictionary relating the form name and the form object.
    # @param self The object pointer.
    def updateAvailableForms(self):
        
        for attrName,attrVal in self.__dict__.items():
            if isinstance(attrVal,RateForm):
                if attrName not in self.__formDict.keys():
                    self.__formDict[attrName] = attrVal
    
    ## Returns the explicit parameters dictionary.
    # 
    # The parameters in this dictionary are intended to be written explicitely
    # in the final call function.
    # @param self The object pointer.
    def getExplicitParameters(self):
        return self.__explicitParam
    
    ## Returns all parameter dictionaries.
    # 
    # 1. Explicit are intended to be written explicitely in the final call function.
    # 2. Global parameters are searched for common values used in the same key in 
    #     multiple reactions.
    # 3. Reaction parameters are only accessible to their reactions.
    # @param self The object pointer.
    def getParameters(self):
        return self.__explicitParam, self.__globalParam, self.__rxnParam
    
    def getOptSpace(self):
        """ Returns the optimization space and bounds for this model.
        
        The function will return a list of OptData with the values and bounds
        for parameters that should be optimized. The parameters can be kinetic
        properties of reactions or initial values of metabolite concentrations.
        
        Args:
            
        Returns:
            List of OptData
        
        """
        
        self.__optSpaceList = []
        
        currParVals = []
        
        if self.__verbose > 1:
            print("Explicit Parameters to be optimized:")
        
        for parName, parDat in self.__explicitParam.items():
            if parDat.lb != parDat.ub:
                self.__optSpaceList.append(OptDat(parDat))
                self.__optSpaceList[-1].type = "explPar"
                self.__optSpaceList[-1].indx = parName
                
                currParVals.append(parDat.val)
                
                if self.__verbose > 1:
                    print(self.__optSpaceList[-1].type,
                          self.__optSpaceList[-1].field,
                          self.__optSpaceList[-1].lb, 
                          self.__optSpaceList[-1].ub, 
                          self.__optSpaceList[-1].indx,
                          parName, parDat.val)
        
        if self.__verbose > 1:
            print("Global Parameters to be optimized:")
        
        for paramKey, parDat in self.__globalParam.items():
            if parDat.lb != parDat.ub:
                self.__optSpaceList.append(OptDat(parDat))
                self.__optSpaceList[-1].type = "glblPar"
                self.__optSpaceList[-1].indx = paramKey
                
                currParVals.append(parDat.val)
                
                if self.__verbose > 1:
                    print(self.__optSpaceList[-1].type,
                          self.__optSpaceList[-1].field,
                          parDat.lb, parDat.ub, paramKey, parDat.val)
        
        if self.__verbose > 1:
            print("Reaction-Specific Parameters to be optimized:")
        
        for rxnID, parDict in self.__rxnParam.items():
            for formKey, parDat in parDict.items():
                if parDat.lb != parDat.ub:
                    self.__optSpaceList.append(OptDat(parDat))
                    self.__optSpaceList[-1].type = "rxnPar"
                    self.__optSpaceList[-1].indx = rxnID
                    
                    currParVals.append(parDat.val)
                    
                    if self.__verbose > 1:
                        print(self.__optSpaceList[-1].type,
                          self.__optSpaceList[-1].field,
                          parDat.lb, parDat.ub, rxnID, parDat.val)
        
        #for met in self.__metList:
        
        if self.__verbose > 0:
            print("Number of dimensions for optimization:",len(self.__optSpaceList))
        
        return self.__optSpaceList, currParVals
    
    def applyOpt(self, x):
        """ Applies new parameter values to reactions. Used in parameter optimization.
        
        Stores all new values in their respective parameter objects. Then, iterates
        over reactions to apply the new parameter values.
        
        Args:
            x (list): The list with new values for model parameters.
        Returns:
            none
        
        """
        
        if not self.__optSpaceList:
            print("No optimization space was defined!")
            return -1
        
        indx = 0
        for optItem in self.__optSpaceList:
            
            #print(indx, x[indx], optItem.type,optItem.indx,
                      #optItem.field,optItem.val,
                      #optItem.lb, optItem.ub)
            
            if optItem.type == "explPar":
                self.__explicitParam[optItem.indx].val = x[indx]
            
            elif optItem.type == "glblPar":
                self.__globalParam[optItem.indx].val = x[indx]
            
            elif optItem.type == "rxnPar":
                self.__rxnParam[optItem.indx][optItem.field].val = x[indx]
            
            indx += 1
        
        
        for rxn in self.__rxnList:
            
            rxnID = rxn.getID()
            
            for formKey in rxn.getKeys():
                
                foundKey = False
                
                # If there is a parameter for this specific reaction, 
                # apply this to the unbound key.
                if rxnID in self.__rxnParam.keys():
                    
                    if formKey in self.__rxnParam[rxnID].keys():
                        rxn.addParameter(formKey,self.__rxnParam[rxnID][formKey].val,
                                         self.__verbose)
                        foundKey = True
                
                # If no specific parameter was found but there is a global 
                # parameter available, apply it now.
                if not foundKey and formKey in self.__globalParam.keys():
                    rxn.addParameter(formKey,self.__globalParam[formKey].val, self.__verbose)
            
    
    def parseRxnIdentifier(self, rxnIdentifier):
        """ Returns a list with reaction ID and reaction index, given an identifier.
        
            The identifier may be a string containing the ID or an int containing the 
            index. The function checks if the ID exists and if the index is within the
            bounds for available reactions, then returns both values.
            
            Args:
                rxnIdentifier (str or int) The identifier, either a string or an int.
            Returns:
                ID,indx
                
                Returns the index of the reaction and the string representing its ID.
        """
        
        if isinstance(rxnIdentifier,str):
            rxnID = rxnIdentifier
            
            # Checks if the metabolite has been defined
            if rxnID in self.__rxnDict.keys():
                rxnIndex = self.__rxnDict[rxnID]
            else :
                print("ERROR: Reaction \"" + rxnID + "\" was not previously defined!")
                return ["", -1]
            
        if isinstance(rxnIdentifier,int):
            rxnIndex = rxnIdentifier
            
            # Sanity Checks
            if rxnIndex < len(self.__rxnList):
                rxnID = self.__rxnList[rxnIndex].getID()
            else:
                print("ERROR: Reaction index " + str(rxnIndex) +  \
                    " is out of bounds. " + str(len(self.__rxnList)) + \
                    " reactions were defined so far.")
                return ["", -1]
        
        return [rxnID, rxnIndex]
    
    ## Returns a list with metabolite ID and metabolite index, given an identifier.
    #
    # The identifier may be a string containing the ID or an int containing the 
    # index. The function checks if the ID exists and if the index is within the
    # bounds for available metabolites, then returns both values.
    # @param self The object pointer.
    # @param metIdentifier The identifier, either a string or an int.
    def parseMetIdentifier(self, metIdentifier):
        if isinstance(metIdentifier,str):
            metID = metIdentifier
            
            # Checks if the metabolite has been defined
            if metID in self.__metDict.keys():
                metIndex = self.__metDict[metID]
            else :
                #print("WARNING: Metabolite \"" + metID + \
                    #"\" was not previously defined. Defining it now!")
                #metIndex = self.addMetabolite(metID, metName)
                print("ERROR: Metabolite \"" + metID + "\" was not previously defined!")
                return ["", -1]
                
        if isinstance(metIdentifier,int):
            metIndex = metIdentifier
            
            # Sanity Checks
            if metIndex < len(self.__metList):
                metID = self.__metList[metIndex].getID()
            else:
                print("ERROR: Metabolite index " + str(metIndex) +  \
                    " is out of bounds. " + str(len(self.__metList)) + \
                    " metabolites were defined so far.")
                return ["", -1]
        
        return [metID, metIndex]
    
    ## Adds a metabolite to the model and returns its index.
    #
    # @param self The object pointer.
    # @param metID The metabolite ID (string).
    # @param metName The name of the new metabolite.
    def addMetabolite(self, metID, metName = "", initVal = 0, fbaMetID = "", metMode = ""):
        
        if (metID in self.__metDict.keys()):
            print("Error: metabolite \"" + metID + "\" was already defined in model!")
            return -1
        
        self.__metList.append(Metabolite(metID,metName,initVal,fbaMetID,metMode))
        metIndex = len(self.__metList) -1
        self.__metDict[metID] = metIndex
        
        if self.getVerbosity() > 1:
            finalStr = "(" +  str(metIndex) + ") Created metabolite \"" + metID + "\"; "
            if fbaMetID:
                finalStr += "FBA ID: " + fbaMetID + "; "
            if metName:
                finalStr += "Name \"" + metName + "\"; "
            finalStr += "Initial value: " + str(initVal)
            print(finalStr)
        
        return metIndex
    
    def addRateForm(self, rateName, rateObj):
        
        if not isinstance(rateObj, RateForm):
            print("Error: object is not a DESolver RateForm!")
            return 1
        
        if rateName in self.__formDict.keys():
            print("Error: rate form already defined in this model!")
            return 1
        
        setattr(self, rateName, rateObj)
        
        self.__formDict[rateName] = rateObj
        
        
    def addReaction(self, rxnID, rxnRateForm, rxnName = ""):
        
        if (rxnID in self.__rxnDict.keys() ):
            print("Error: reaction \"" + rxnID + "\" was already defined in this model!")
            return 1
        
        elif rxnRateForm in self.__formDict.keys():
            
            self.__rxnList.append(Reaction(rxnID,rxnName))
            rxnIndex = len(self.__rxnList) -1
            self.__rxnDict[rxnID] = rxnIndex
            self.__rxnList[rxnIndex].setRateForm(self.__formDict[rxnRateForm], rxnRateForm)
            
            if self.getVerbosity() > 1:
                finalStr = "(" +  str(rxnIndex) + ") Created reaction \"" + rxnID + "\"; "
                if rxnName:
                    finalStr += "Name \"" + rxnName + "\"; "
                finalStr += "Rate form: " + rxnRateForm
                print(finalStr)
            
            return rxnIndex
            
        else:
            
            print("ERROR: Reaction form \"" + rxnRateForm + "\" is not defined in this model.")
            return -1
    
    def addSubstrate(self, rxnIdentifier, rxnFormKey, metIdentifier, stoich=0):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        metID, metIndex = self.parseMetIdentifier(metIdentifier)
        
        if(metID == "" or rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        retCode = self.__rxnList[rxnIndex].addSubstrate(rxnFormKey, metID, stoich)
        
        # Code 2 indicates no changes were made
        if retCode != 2:
            if rxnIndex in self.__metList[metIndex].getReactions():
                if self.getVerbosity() > 0:
                    print("WARNING: Metabolite \"" + metID + \
                    "\" has already been associated with reaction \"" + rxnID + "\".")
            self.__metList[metIndex].addReaction(rxnIndex)
            
        if not rxnFormKey in self.__rxnList[rxnIndex].getKeys():
            if self.getVerbosity() > 0:
                print("WARNING: Key \"" + rxnFormKey + "\" is not defined in the \""\
                + self.__rxnList[rxnIndex].getRateFormName() + \
                  "\" rate form for reaction \"" + rxnID + "\".")
    
    def addProduct(self, rxnIdentifier, rxnFormKey, metIdentifier, stoich=0):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        metID, metIndex = self.parseMetIdentifier(metIdentifier)
        
        if(metID == "" or rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        retCode = self.__rxnList[rxnIndex].addProduct(rxnFormKey, metID, stoich)
        
        # Code 2 indicates no changes were made
        if retCode != 2:
            if rxnIndex in self.__metList[metIndex].getReactions():
                if self.getVerbosity() > 0:
                    print("WARNING: Metabolite \"" + metID + \
                    "\" has already been associated with reaction \"" + rxnID + "\".")
            self.__metList[metIndex].addReaction(rxnIndex)
            
        if not rxnFormKey in self.__rxnList[rxnIndex].getKeys():
            if self.getVerbosity() > 0:
                print("WARNING: Key \"" + rxnFormKey + "\" is not defined in the \""\
                + self.__rxnList[rxnIndex].getRateFormName() + \
                  "\" rate form for reaction \"" + rxnID + "\".")
    
    def addParameter(self, rxnIdentifier, rxnFormKey, value, lb=0, ub=0, unit="", parName=""):
        
        
        if rxnIdentifier == "Explicit":
            self.__explicitParam[parName] = Parameter(parName, value,\
                unit, rxnFormKey, lb, ub)
            return 0
        
        elif rxnIdentifier == "Global":
            if rxnFormKey in self.__globalParam.keys():
                print("ERROR: Identical FromKeys are not allowed among global parameters!")
                print("--> Key \"" + rxnFormKey + "\" was defined more than once with no associated reaction.")
            
            self.__globalParam[rxnFormKey] = Parameter(parName, value, \
                unit, rxnFormKey, lb, ub)
            if self.getVerbosity() > 1:
                print("Added " + rxnFormKey + " to global parameters dictionary")
            return 0
        
        # For a reaction-specific parameter, check the reaction identifier 
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        if(rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        # If the parameter is added directly to a reaction object, we keep the
        # data for the model's parameter list.
        
        if not rxnID in self.__rxnParam.keys():
            self.__rxnParam[rxnID] = dict()
        self.__rxnParam[rxnID][rxnFormKey] = Parameter(parName, value,\
                        unit, rxnFormKey, lb, ub)
        
        # Add the parameter to the reaction object
        retCode = self.__rxnList[rxnIndex].addParameter(rxnFormKey, value)
        
        # Code 2 indicates no changes were made
        if retCode != 2:
            if not rxnFormKey in self.__rxnList[rxnIndex].getKeys():
                if self.getVerbosity() > 0:
                    print("WARNING: Key \"" + rxnFormKey + "\" is not defined in the \""\
                    + self.__rxnList[rxnIndex].getRateFormName() + \
                      "\" rate form for reaction \"" + rxnID + "\".")
    
    def rmMetabolite(self, metIdentifier):
        
        metID, metIndex = self.parseMetIdentifier(metIdentifier)
        
        if(metID == ""):
            print("--> No action will be performed!")
            return 1
        
        rxnSet = self.__metList[metIndex].getReactions()
        
        if (len(rxnSet) > 0):
            if self.getVerbosity() > 0:
                print("WARNING: Metabolite associated with model reactions.")
                print("--> No action will be performed!")
            return rxnSet
        else:
            del self.__metDict[metID]
            del self.__metList[metIndex]
            
    
    def rmSubstrate(self, rxnIdentifier, rxnFormKey):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        if(rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        metID = self.__rxnList[rxnIndex].rmSubstrate(rxnFormKey)
        
        if metID == 1:
            
            print("--> No action will be performed!")
            return 1
            
        else:
            
            metIndex = self.__metDict[metID]
            
            if not self.__rxnList[rxnIndex].hasMetabolite(metID):
                self.__metList[metIndex].rmReaction(rxnIndex)
        
    def rmProduct(self, rxnIdentifier, rxnFormKey):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        if(rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        metID = self.__rxnList[rxnIndex].rmProduct(rxnFormKey)
        
        if metID == 1:
            
            print("--> No action will be performed!")
            return 1
            
        else:
            
            metIndex = self.__metDict[metID]
            
            if not self.__rxnList[rxnIndex].hasMetabolite(metID):
                self.__metList[metIndex].rmReaction(rxnIndex)
    
    def rmParameter(self, rxnIdentifier, rxnFormKey):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(rxnIdentifier)
        
        if(rxnID == ""):
            print("--> No action will be performed!")
            return 1
        
        returnVal = self.__rxnList[rxnIndex].rmParameter(rxnFormKey)
        
        if returnVal != 0:
            
            print("--> No action will be performed!")
            return 1
    
    def addDependent(self, parent, dependent):
        
        rxnID, rxnIndex = self.parseRxnIdentifier(parent)
        rxnIDdep, rxnIndexdep = self.parseRxnIdentifier(dependent)
        
        if (rxnIndex < 0):
            print("ERROR! Could not find parent reaction:",parent)
            sys.exit(2)
        if (rxnIndexdep < 0):
            print("ERROR! Could not find dependent reaction:",dependent)
            sys.exit(2)    
        
        self.__rxnList[rxnIndex].addDependent(rxnIndexdep)
        
        self.__rxnList[rxnIndex].addDepKeysVals(self.__rxnList[rxnIndexdep].getKeysVals())
        
    
    def getInitVals(self):
        
        initVals = []
        
        for met in self.__metList:
            
            initVals.append(float(met.getInitVal()))
        
        return initVals
    
    def setFBAModel(self, fbaModel, fileName = ""):
        if isinstance(fbaModel,str):
            try:
                with open(fbaModel, "rb") as infile:
                    self.__cobraModel = load(infile)
            except:
                self.__cobraModel = 0
                print("ERROR: Could not load file \"" + 
                          fbaModel + "\" as a pickled cobra model.")
                return(1)
        else :
            self.__cobraModel = fbaModel
            
        if isinstance(self.__cobraModel,cobra.core.Model):
            
            print("COBRA model loaded successfully: " + str(fbaModel) )
            
            self.__hybrid_ODE_FBA = True
            
            self.compileFbaEquivRxns()
            
            if fileName != "":
                self.parseFbaIgnoreRxnList(fileName)
            
            self.parseConnRxns()
            
            # Optimize FBA model with updated constraints
            self.runFBA()
            
            # Place GR on ODE model for dilution term.
            self.setGR(self.__cobraModelSolution.f)
            
            # Updates connecting fluxes for ODE model
            self.calcConnFlux()
            
            return(0)
            
        else:
            self.__cobraModel = 0
            
            print("ERROR: file \"" + 
                          str(fbaModel) + "\" is not a pickled cobra model.")
            
            return(2)
    
    def getFBAModel(self):
        return self.__cobraModel
    
    def setGR(self, newGR):
        print("Setting Growth Rate of ODE model to:", newGR)
        self.__gr = newGR
    
    def getGR(self):
        return self.__gr
    
    def setUnitConvFactor(self, newfactor):
        self.__unitConvFactor = newfactor
    
    def getUnitConvFactor(self):
        return self.__unitConvFactor
    
    def getFBASolution(self):
        return self.__cobraModelSolution
    
    ## Compiles a set with all FBA equivalent reaction IDs
    # 
    # This function goes through all Reaction Objects and gets 
    # their equivalent FBA reaction IDs.
    # Those will not be used as connecting reactions, even though
    # they operate in metabolites traked by the ODE model.
    def compileFbaEquivRxns(self):
        for rxn in self.__rxnList:
            for fbaid in [ setItem[0] for setItem in rxn.getFBAEquiv()]:
                self.__fbaEquivRxns.add(fbaid)
    
    def parseFbaIgnoreRxnList(self, fileName):
        """ Parse a CSV reaction list with one reaction ID per line
        
        The reaction IDs read from the file will compose a list of reactions 
        from the FBA model which will NOT be considered as connecting reactions.
        It should cover all reactions covered by the ODE model.
        
        Args:
            fileName (str): The file name of the CSV file with reaction IDs.
        Returns:
            Zero (0) if all goes fine.
        
        """
        
        with open(fileName) as csvfile:
            # Detects the CSV format automatically
            #fileDialect = csv.Sniffer().sniff(csvfile.readline())
            #csvfile.seek(0)
            
            # Defines the dictionarry creator.
            #parDictReader = csv.DictReader(csvfile, dialect=fileDialect)
            parDictReader = csv.DictReader(csvfile, delimiter=",")
            
            for row in parDictReader:
                rxnID = row["ReactionID"]
                
                self.__fbaEquivRxns.add(rxnID)
                
        return(0)
    
    ## Parse a CSV parameter file and stores parameters in the reaction objects.
    #
    # The function reads in all parameters from a CSV file and applies reaction
    # specific parameters first, then global parameters (parameters that have no
    # value in the ReactionID column). The latter are applied to all reactions
    # which have unbound keys.
    # 
    # @param self The object pointer.
    # @param parameterFileName The file name of the CSV file with parameters.
    def parseParameters(self, parameterFileName):
        
        # The function builds two dictionaries, one key-ed by reaction name, 
        # containing another dictionary (keyed by FormKey) of objects which 
        # hold names, values and units of that parameter. The other dictionary is 
        # keyed by rxn FormKeys, and contains a similar dictionary as the previous 
        # one, which holds objects with parameter names, units and values.
        
        if self.__verbose > 0:
            print("Parsing parameter file " + parameterFileName + " ...")
        
        with open(parameterFileName) as csvfile:
            # Detects the CSV format automatically
            fileDialect = csv.Sniffer().sniff(csvfile.readline())
            csvfile.seek(0)
            
            # Defines the dictionarry creator.
            parDictReader = csv.DictReader(csvfile, dialect=fileDialect)
            #parDictReader = csv.DictReader(csvfile, delimiter=",")
            
            for row in parDictReader:
                
                
                rxnID = row["ReactionID"]
                paramKey = row["FormKey"]
                parVal = row["parVal"]
                parName = row["parName"]
                parUnit = row["unit"]
                
                parLB = 0
                parUB = 0
                if ("optimizeMin" in row.keys()) and ("optimizeMax" in row.keys()):
                    parLB = row["optimizeMin"]
                    parUB = row["optimizeMax"]
                    
                    if parLB == "":
                        parLB = 0
                    
                    if parUB == "":
                        parUB = 0
                
                if rxnID == "Explicit":
                    self.__explicitParam[parName] = Parameter(parName, parVal,\
                        parUnit, paramKey, parLB, parUB)
                
                elif rxnID == "Global":
                    if paramKey in self.__globalParam.keys():
                        print("ERROR: Identical FromKeys are not allowed among global parameters!")
                        print("--> Key \"" + paramKey + "\" was defined more than once with no associated reaction.")
                    
                    self.__globalParam[paramKey] = Parameter(parName, parVal, \
                        parUnit, paramKey, parLB, parUB)
                    if self.getVerbosity() > 1:
                        print("Added " + paramKey + " to global parameters dictionary")
                else :
                    if not rxnID in self.__rxnParam.keys():
                        self.__rxnParam[rxnID] = dict()
                    
                    if paramKey in self.__rxnParam[rxnID].keys():
                        if self.getVerbosity() > 0:
                            print("WARNING: Parameter file has a new value for key \""\
                                + paramKey + "\" in reaction \"" + rxnID + "\".")
                            print("--> Only the last value was kept!")
                    self.__rxnParam[rxnID][paramKey] = Parameter(parName, parVal,\
                        parUnit, paramKey, parLB, parUB)
                    
                    
        if self.getVerbosity() > 0:
            print("All parameters were loaded from file \"" + parameterFileName + "\".")
        
        self.parameterizeRxns()
    
    def parameterizeRxns(self):
        """ Applies the parameters to reactions, to form final rate forms.
        
        For each reaction, checks all unbound keys.
        Checking only unbound keys assures that parameters passed manually 
        will not be overwritten.
        
        Args:
            
        Returns:
            none
        
        """
        
        for rxn in self.__rxnList:
            
            rxnID = rxn.getID()
            if self.getVerbosity() > 1:
                print("Checking parameters for reaction \"" + rxnID + "\"")
            
            for unbdKey in rxn.getUnboundKeys():
                
                foundKey = False
                
                # If there is a parameter for this specific reaction, 
                # apply this to the unbound key.
                if rxnID in self.__rxnParam.keys():
                    if unbdKey in self.__rxnParam[rxnID].keys():
                        rxn.addParameter(unbdKey,self.__rxnParam[rxnID][unbdKey].val)
                        foundKey = True
                
                # If no specific parameter was found but there is a global 
                # parameter available, apply it now.
                if not foundKey and unbdKey in self.__globalParam.keys():
                    rxn.addParameter(unbdKey,self.__globalParam[unbdKey].val)
            
            # Checks for unused parameters
            if rxnID in self.__rxnParam.keys():
                unusedKeys = set(self.__rxnParam[rxnID].keys()) - set(rxn.getKeys())
                if len(unusedKeys) > 0:
                    if self.getVerbosity() > 0:
                        print("WARNING: Keys defined but not used in reaction form for \""\
                        + rxnID + "\": " + str(unusedKeys) )
                    for unusedKey in unusedKeys:
                        rxn.addParameter(unusedKey,self.__rxnParam[rxnID][unusedKey].val)
    
    def writeParameters(self, fileName = "parameters_OptResults.csv"):
        """ Writes a file with the parameters given to the model.
        
        The function is intended to output the results of a parameter 
        optimization run.
        
        Args:
            fileName (str): Name of file for parameter output.
        Returns:
        
        """
        
        tmpList = []
        
        for parName, parDat in self.__explicitParam.items():
            if parDat.lb != parDat.ub:
                tmpList.append(parDat)
        
        for paramKey, parDat in self.__globalParam.items():
            if parDat.lb != parDat.ub:
                tmpList.append(parDat)
        
        for rxnID, parDict in self.__rxnParam.items():
            for formKey, parDat in parDict.items():
                if parDat.lb != parDat.ub:
                    tmpList.append(parDat)
        
        for parName, parDat in self.__explicitParam.items():
            if parDat.lb == parDat.ub:
                tmpList.append(parDat)
        
        for paramKey, parDat in self.__globalParam.items():
            if parDat.lb == parDat.ub:
                tmpList.append(parDat)
        
        for rxnID, parDict in self.__rxnParam.items():
            for formKey, parDat in parDict.items():
                if parDat.lb == parDat.ub:
                    tmpList.append(parDat)
        
        outfile = open(fileName,"w")
        
        outfile.write("parName,parVal,unit,FormKey,ReactionID,optimizeMin,optimizeMax\n")
        
        for par in tmpList:
            outfile.write(par.name + ",")
            outfile.write(par.val + ",")
            outfile.write(par.unit + ",")
            outfile.write(par.formKey + ",")
            outfile.write(par.lb + ",")
            outfile.write(par.ub + "\n")
        
        outfile.close()
    
    def parseMetabolites(self, metabolitesFileName):
        """Parse a CSV metabolite file and creates the Metabolite objects.
        #
        # The function reads in all information from a CSV file and creates the 
        # Metabolite objects with names and initial concentrations.
        # 
        # @param self The object pointer.
        # @param metabolitesFileName The file name of the CSV file with metabolites.
        """
        
        if self.__verbose > 0:
            print("Parsing metabolites file " + metabolitesFileName + " ...")
        
        with open(metabolitesFileName) as csvfile:
            # Detects the CSV format automatically
            fileDialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            
            # Defines the dictionarry creator.
            metDictReader = csv.DictReader(csvfile, dialect=fileDialect)
            
            for row in metDictReader:
                
                metID = row["metID"]
                
                try:
                    fbaMetID = row["fbaMetID"]
                except:
                    fbaMetID = ""
                
                metName = row["metName"]
                
                try:
                    initVal = float(row["initVal"])
                except:
                    if self.getVerbosity() > 0:
                        print("WARNING: The initial value for metabolite \"" + 
                          metID + "\" cannot be converted to a float!")
                    initVal = row["initVal"]
                
                try:
                    metMode = row["mode"]
                except:
                    metMode = ""
                
                if self.addMetabolite(metID, metName = metName, 
                                      initVal = initVal, 
                                      fbaMetID = fbaMetID, 
                                      metMode=metMode) < 0:
                    print("--> Skipping metabolite!")
                    
    
    def parseConnFlux(self, connectFluxFileName):
        """Parse a CSV flux file and load connecting fluxes.
        #
        # The function reads in all information from a CSV file and loads the 
        # connecting fluxes from an FBA model into this ODE model, checking if
        # the reffered metabolites exist, and finally setting to ZERO the 
        # connecting flux of metabolites not found in the CSV file.
        # 
        # @param self The object pointer.
        # @param connectFluxFileName The file name of the CSV file with fluxes.
        """
        
        print("Parsing connecting fluxes file " + connectFluxFileName + " ...")
        
        with open(connectFluxFileName) as csvfile:
            # Detects the CSV format automatically
            fileDialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            
            # Defines the dictionarry creator.
            metDictReader = csv.DictReader(csvfile, dialect=fileDialect)
            
            for row in metDictReader:
                
                metID = row["metID"]
                
                try:
                    conFlux = float(row["flux"])
                    
                except:
                    print("ERROR: The connecting flux value for metabolite \"" + 
                      metID + "\" cannot be converted to a float!")
                    
                    return(1)
                
                if metID in self.__metDict.keys():
                    metIndex = self.__metDict[metID]
                    
                    self.__metList[metIndex].setConnFlux(conFlux)
        
        return(0)
    
    ## Lists FBA reactions which connect both models.
    #
    # This funciton checks which FBA reactions operate on metabolites
    # tracked by the ODE model.
    # 
    # @param self Self pointer
    def parseConnRxns(self):
        
        print("Parsing connected reactions between models.")
        
        for odemet in self.__metList:
            # Makes sure the connecting reactions list is empty,
            # in case this function is ran repeatedly.
            odemet.cleanConnections()
        
        for rxnIndx in range(len(self.__cobraModel.reactions)):
            
            rxn = self.__cobraModel.reactions[rxnIndx]
            
            if (rxn.id in self.__fbaEquivRxns):
                continue
            
            for odemet in self.__metList:
                
                if odemet.getFBAID() in [ met.id for met in rxn.metabolites ] :
                    
                    odemet.addFBAConnection(rxnIndx, rxn.get_coefficient(odemet.getFBAID()) )
        
        #print("The following connecting reactions were found:")
        #for odemet in self.__metList:
            #print("ODE Metabolite: " + odemet.getID() )
            
            #for rxnTuple in odemet.getConnRxns():
                #print(str(rxnTuple[0]) + " (" + self.__cobraModel.reactions[rxnTuple[0]].id + "): " + str(rxnTuple[1]))
        
    
    def calcConnFlux(self):
        
        if self.getVerbosity() > 0:
            print("Calculating connecting fluxes...")
        
        if self.getVerbosity() > 1:
            for odemet in self.__metList:
                odemet.calcConnFlux(self.__cobraModelSolution, self.__cobraModel)
        else:
            for odemet in self.__metList:
                odemet.calcConnFlux(self.__cobraModelSolution)
                
    def updateFBAConn(self, rxnFluxDict):
        
        fbaOverlapFlux = dict()
        
        # Place ODE fluxes as constraints on FBA model
        
        # We make sure all equivalent reaction fluxes are added before
        # applying lower and upper bounds.
        for oderxn in self.__rxnList:
            for setItem in oderxn.getFBAEquiv():
                if setItem[1] in fbaOverlapFlux.keys():
                    fbaOverlapFlux[setItem[0]] += rxnFluxDict[oderxn.getID()]*setItem[1]
                else:
                    fbaOverlapFlux[setItem[0]] = rxnFluxDict[oderxn.getID()]*setItem[1]
        
        for key, val in fbaOverlapFlux.items():
            self.__cobraModel.reactions.get_by_id(key).lower_bound = val*0.9
            self.__cobraModel.reactions.get_by_id(key).upper_bound = val*1.1
            print(key,":",val)
        
        with open("test_fluxBounds_inDesolver.csv","w") as outfile:
            for rxn in self.__cobraModel.reactions:
                outfile.write(rxn.id + "," + str(rxn.lower_bound) + "," + str(rxn.upper_bound) + "\n")
        
        # Optimize FBA model with updated constraints
        self.runFBA()
        
        # Place GR on ODE model for dilution term.
        self.setGR(self.__cobraModelSolution.f)
        
        # Updates connecting fluxes for ODE model
        self.calcConnFlux()
        
        return(0)
    
    def runFBA(self):
        
        self.__cobraModelSolution = cobra.flux_analysis.optimize_minimal_flux(self.__cobraModel)
    
    def __getDepKeyVals(self,rxnIndx):
        
        retDict = dict()
        
        # Gets the keys and values for the first reaction. 
        retDict.update(self.__rxnList[rxnIndx].getKeysVals())
        
        # Gets keys and values for all dependent reactions, recursively.
        for depID in self.__rxnList[rxnIndx].getDependentRxns():
            
            retDict.update( self.__getDepKeyVals(depID) )
        
        return retDict
    
    def prepModel(self):
        
        for rxn in self.__rxnList:
            
            # Sets the numbers used to identify each reaction's rate.
            # This is only done for reactions whose result has not been
            # explicitly set by the user during model constroction.
            if rxn.getResult() == 0 and rxn.getNumber() == 0:
                self.numMainRxns += 1
                rxn.setNumber(self.numMainRxns)
            
            # Builds recursively a complete list of all keys-values 
            # for all dependent reactions.
            for depID in rxn.getDependentRxns():
                
                rxn.addDepKeysVals(self.__getDepKeyVals(depID))
        
        return 0
    
    def checkModel(self):
        
        if self.getVerbosity():
            print("\nChecking model...")
        
        ## Create an object to hold errors/problems with the error msg
        ## that originated them, maybe with an associated obj (Rxn, met, param?).
        #if self.__modelStatus != 0:
            #print("ERRORS occurred during model building!")
            #return self.__modelStatus
        
        for rxn in self.__rxnList:
            
            if rxn.getResult() != 0 and rxn.getNumber() == 0:
                # If it is an auxiliary reaction, it cannot have 
                # products or substrates assigned to it, since it would lead
                # to code being generated for it mistakenly.
                
                if (len(rxn.getProducts()) > 0 or len(rxn.getSubstrates()) > 0) and rxn.getCheckRxn():
                    print("ERROR: Auxiliary Reaction cannot have Products/Reactants, \
only parameters. Offending Reaction:", rxn.getID())
                    return 1
            
            
            undbKeySet = rxn.getUnboundKeys()
            
            if len(undbKeySet) > 0:
                print("ERROR: Unbound keys found in reaction \"" + \
                    rxn.getName() + "\" (" + rxn.getID() + ") :" + \
                    str(undbKeySet) )
                return 1
            
        return 0
    
    def checkModelJac(self):
        
        if self.getVerbosity():
            print("\nChecking model Jacobian...")
        
        ## Gets all free metabolite IDs
        metIDList = [met.getID() for met in self.__metList]
        
        ## Gets all dependent metabolite IDs
        depMetIDSet = reduce(set.union,[met.getDependMets() for met in self.__metList])
        
        ## Maps dependent metabolites to their respective parent (free metabolite).
        depMetMapping = dict()
        for met in self.__metList:
            for depMetID in met.getDependMets():
                if depMetID in depMetMapping.keys():
                    print("ERROR! Dependent Metabolite cannot be dependent on more than one free metabolite.")
                    sys.exit(2)
                depMetMapping[depMetID] = met.getID()
        
        # Checks if the reaction has gradient component for all metabolites (free
        # or dependent) defined in its rate forms and in its dependent reactions.
        for rxn in self.__rxnList:
            
            #print("\n",rxn.getID())
            
            #if rxn.getID().startswith("PpsA"):
                #if rxn.getID() == "PpsA":
                    #print("\n ---> WARNING: no jacobian terms for PpsA!!\n")
                #continue
            #if rxn.getID() == "PpsA_ZADPi":
                #print("rxn.getKeysVals():",rxn.getKeysVals())
                #print(" rxn.getDependentKeysVals():", rxn.getDependentKeysVals())
                #print("rxn.getDependentRxns():",rxn.getDependentRxns())
                #print("rxn.getGrad():")
                #print(rxn.getGrad())
            
            for key,val in rxn.getDependentKeysVals().items():
                if key in rxn.getKeysVals().keys():
                    if val != rxn.getKeysVals()[key]:
                        print("ERROR! Keys repeated in dependent reactions CANNOT have different values!")
                        print("Reaction:",rxn.getID(),"; Key:",key,"Value:",val,"; Value in dependent rxn:",rxn.getKeysVals()[key])
                        sys.exit(2)
            
            metSet = set()
            depMetSet = set()
            
            depRxnMetSet = set()
            depRxnDepMetSet = set()
            
            # First we build sets of values associated with keys for the reaction
            # and dependent reactions.
            for possibleID in [ rxn.getKeysVals()[key] for key in rxn.getKeys()]:
                if possibleID in metIDList:
                    metSet.add(possibleID)
                    
                if possibleID in depMetIDSet:
                    depMetSet.add((possibleID,depMetMapping[possibleID]))
            
            for possibleID in rxn.getDependentKeysVals().values():
                if possibleID in metIDList:
                    depRxnMetSet.add(possibleID)
                if possibleID in depMetIDSet:
                    depRxnDepMetSet.add((possibleID,depMetMapping[possibleID]))
            
            depMetSetIter = depMetSet.copy()
            depRxnDepMetSetIter = depRxnDepMetSet.copy()
            
            # Next, we get the set of metabolites for which we have gradient 
            # components in this reaction.
            gradMetSet = set()
            for metID in [gradDat[0] for gradDat in rxn.getGrad()]:
                if metID in gradMetSet:
                    print("ERROR! More than one gradient component was defined \
for the same metabolite!\n---> If one metabolite has many roles (substrate, \
allosteric activator, inhibitor), create a derived rate for for this reactions \
with a specialized gradient component!")
                    print("Reaction:",rxn.getID(), "; met:", metID)
                else:
                    gradMetSet.add(metID)
            
            # Now we compare  the sets and check if any key was associated with
            # a free or dependent metabolite for which no gradient component
            # was defined:
            
            # 1) For main reactions (NOT dependent reactions), gradients can be defined
            # with respect to dependent metablites. The chain rule will be applied
            # automatically.
            if rxn.getNumber() != 0:
                for dep,free in depMetMapping.items():
                    if dep in gradMetSet:
                        if free in gradMetSet:
                            print("ERROR! Gradient components were defined \
for BOTH Free and Dependent Metabolites in a main reaction!\n---> Only one should be defined.")
                            print("Reaction:",rxn.getID(), "; met (dependent, free):", dep,free)
                            sys.exit(2)
                        if len(rxn.getDependentRxns()) > 0:
                            print("ERROR! Gradient component(s) defined \
for Dependent Metabolite in a main reaction WITH dependent reactions!\n---> \
This is only allowed when the main reaction has NO dependent reactions.")
                            print("Reaction:",rxn.getID(), "; met (dependent, free):", dep,free)
                            sys.exit(2)
                        gradMetSet.remove(dep)
                        gradMetSet.add(free)
            
            else:
                for dep,free in depMetMapping.items():
                    if dep in gradMetSet:
                        print("ERROR! Gradient components were defined \
for Dependent Metabolites in a dependent reaction!\n---> Only gradients \
with respect to free metabolites can be defined in dependent reactions!")
                        print("Reaction:",rxn.getID(), "; met (dependent, free):", dep,free)
                        sys.exit(2)
            
            #if rxn.getID() == "MaeB_Q_Er":
                #print("metSet:",metSet)
                #print("depMetSet:",depMetSet)
                #print("depRxnMetSet:",depRxnMetSet)
                #print("depRxnDepMetSet:",depRxnDepMetSet)
                #print("gradMetSet:",gradMetSet)
                
            # 2) For free metabolites, the comparison can be done simply by removing
            # metabolites form the sets we built initially.
            metSet.difference_update(gradMetSet)
            depRxnMetSet.difference_update(gradMetSet)
            
            # 3) For dependent metabolites, we need to check if their parent metabolite
            # has a gradient component.
            for item in depMetSetIter:
                if item[1] in gradMetSet:
                    depMetSet.remove(item)
            
            for item in depRxnDepMetSetIter:
                if item[1] in gradMetSet:
                    depRxnDepMetSet.remove(item)
            
            if len(metSet) > 0 or len(depMetSet):
                print("ERROR! Free and/or Dependent Metabolites are used, but no \
gradient component was defined!")
                print("Reaction:",rxn.getID(),"; Free metabolites:", metSet,
                      "; Dependent metabolites",depMetSet)
                sys.exit(2)
            
            if len(depRxnMetSet) > 0 or len(depRxnDepMetSet):
                print("ERROR! Free and/or Dependent Metabolites are used IN DEPENDENT \
REACTIONS, but no gradient component was defined!")
                print("Reaction:",rxn.getID(),"; Free metabolites:", depRxnMetSet,
                      "; Dependent metabolites",depRxnDepMetSet)
                sys.exit(2)
        
        return 0
    
    def buildCobraModel(self,modelName="ode_model",checkModel=True):
        
        self.prepModel()
        
        if checkModel:
            if self.checkModel() != 0:
                print("--> Call function cannot be built!")
                return 1
        
        if self.__verbose > 1:
            print("Building COBRA model...")
        
        ## Model created from the example code in cobrapy examples.
        
        cobra_model = cobra.Model(modelName)
        
        for rxn in self.__rxnList:
            
            if rxn.getNumber() == 0:
                continue
            
            if self.__verbose > 1:
                print("Adding reaction", rxn.getID(), "to cobra model.")
            
            reaction = cobra.Reaction(rxn.getID() )
            reaction.name = rxn.getName()
            reaction.subsystem = 'Metabolism'
            reaction.lower_bound = 0.  # This is the default
            reaction.upper_bound = 1000.  # This is the default
            reaction.objective_coefficient = 0. # this is the default
            
            
            # We need to create metabolites as well. If we were using an existing model, we
            # could use get_by_id to get the apporpriate Metabolite objects instead.
            
            cMetsDict = dict()
            
            mets = set(rxn.getSubstrates()).union( set(rxn.getProducts()))
            
            for met in mets:
                stoic = rxn.getStoichiometry(met)
                
                metIndx = self.__metDict[met]
                
                cMet = cobra.Metabolite(met, formula="", name= self.__metList[metIndx].getName(), compartment='c')
                
                cMetsDict[cMet] = stoic
                
            # Adding metabolites to a reaction requires using a dictionary of the
            # metabolites and their stoichiometric coefficients. A group of metabolites can
            # be added all at once, or they can be added one at a time.
            
            reaction.add_metabolites(cMetsDict)
            
            cobra_model.add_reaction(reaction)
            
        
        return cobra_model
    

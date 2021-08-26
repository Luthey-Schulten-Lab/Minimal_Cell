"""
Author: David Bianchi

A file to define all of the reactions in the Lipid Module
"""

import numpy as np
import pandas as pd
from collections import defaultdict, OrderedDict

import odecell

from in_out_two import calcCellVolume
#import defSrcSinkRxns
#import defLipCentNucTransMetRxns as defLipCentNucTransMetRxns
import defLipCentNucTransMetRxns_two as defLipCentNucTransMetRxns
# Import the lipid patch
import lipid_patch_Zane_two as lipid_patch
import PPP_patch_Zane as ppp_patch

### Constants
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

countToMiliMol = 1000/(NA*V)

#def calcCellVolume(pmap):
    
#    SurfaceArea = pmap['CellSA']
    
#    cellRadius = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
    
#    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
#     print('Volume',cellVolume)
    
#    return cellVolume

def partTomM(particles,pmap):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's
    r_cell = 200e-9 # 200 nm radius, 400 nm diameter
    V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell
    
    cellVolume = calcCellVolume(pmap)

    conc = (particles*1000.0)/(NA*cellVolume)

    return conc


def reptModel(model):
    """
    Report on the constructed hybrid model - but probably would only want to do after the first time step

    Arguments: 
    model (model obj.): The ODE Kinetic Model

    Returns:

    None
    """

    dictTypes = defaultdict(int)
    typeList = ["Transcription","Translation","Degradation"]

    for rxn in model.getRxnList():
        
        if rxn.getResult():
            # If an explicit result has been set, this is a dependent reaction.
            dictTypes["Dependent reactions"] += 1
            continue
        
        for rxntype in typeList:
            if rxntype in rxn.getID():
                dictTypes[rxntype] += 1

                
    #print( "There are {} ODEs in the model:".format(len(model.getRxnList())) )

    outList = list(dictTypes.items())
    outList.sort(key=lambda x: x[1], reverse=True)
    for key,val in outList:
        print("{:>20} :   {}".format(key,val) )
        return 0
    
    return


# Rather than using the rate forms in the SBML for the metabolism, this allows us to
# use forward and reverse kcats rather than the geometric mean form.

def Enzymatic(subs,prods):
        
    def numerator(subs,prods):
        
        subterm = [ '( $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subNumer = ' * '.join(subterm)
        
        prodterm = [ '( $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodNumer = ' * '.join(prodterm)
        
        numerator = '( ' + '$kcatF * ' + subNumer + ' - ' + '$kcatR * ' + prodNumer + ' )'
        return numerator
    
    def denominator(subs,prods):
        
        subterm = [ '( 1 + $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subDenom = ' * '.join(subterm)
        
        prodterm = [ '( 1 + $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodDenom = ' * '.join(prodterm)
        
        denominator = '( ' + subDenom + ' + ' + prodDenom + ' - 1 )'
        return denominator
        
    rate = '$onoff * $Enzyme * ( ' + numerator(subs,prods) + ' / ' + denominator(subs,prods) + ' )'
    
    return rate


def defineRxns(model,pmap):
    """
    Define all of the reactions and rateforms needed for the current module to an existing module.

    """

    ### ADD AN ADDITIONAL FILE TO DEFINE METABOLIC RXNS, as well as a patch for some literature values for lipid metabolism
    defLipCentNucTransMetRxns.addReactionsToModel(model,pmap)
    lipid_patch.main(model,pmap)
    ppp_patch.main(model,pmap)
    print("Reactions defined")

    ### ADD ATPase
    #RateForm
#     ATPaseRateLaw = Enzymatic(2,1)

#     RateName = 'ATPase_Rate'

#     model.addRateForm(RateName, odecell.modelbuilder.RateForm(ATPaseRateLaw))
#     model.updateAvailableForms()

#     rxnIndx = model.addReaction('ATPase','ATPase_Rate','ATP synthase')
#     model.addParameter('ATPase','Enzyme',100*countToMiliMol)
#     model.addParameter('ATPase','kcatF',20)

#     model.addParameter('ATPase','kcatR',285)

#     model.addSubstrate('ATPase','Sub1','M_adp_c')
#     model.addParameter('ATPase','KmSub1',0.1)
#     model.addSubstrate('ATPase','Sub2','M_pi_c')
#     model.addParameter('ATPase','KmSub2',4.2)

#     model.addProduct('ATPase','Prod1','M_atp_c')
#     model.addParameter('ATPase','KmProd1',0.6)

#     model.addParameter('ATPase', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch") 
    
    
#     PitRateLaw = Enzymatic(2,3)

#     RateName = 'Piabc_Rate'

#     model.addRateForm(RateName, odecell.modelbuilder.RateForm(PitRateLaw))
#     model.updateAvailableForms()

#     rxnIndx = model.addReaction('PIabc','Piabc_Rate','Pi transport')
#     model.addParameter(rxnIndx,'Enzyme',56*countToMiliMol)

#     model.addParameter('PIabc','kcatF', 50)
#     model.addParameter('PIabc','kcatR',0)

#     pi_e_conc = 140

#     model.addParameter('PIabc','Sub1',pi_e_conc)
#     model.addParameter('PIabc','KmSub1',0.0031)
#     model.addSubstrate('PIabc','Sub2','M_atp_c')
#     model.addParameter('PIabc','KmSub2',0.023)

#     model.addProduct('PIabc','Prod1','M_pi_c',stoich=2)
#     model.addParameter('PIabc','KmProd1',0.02)
#     model.addParameter('PIabc','Prod2','M_pi_c')
#     model.addParameter('PIabc','KmProd2',0.385)
#     model.addProduct('PIabc','Prod3','M_adp_c')
#     model.addParameter('PIabc','KmProd3',0.654)

#     model.addParameter('PIabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")


#   Add the Source and Sink Reactions to the model object

    #defSrcSinkRxns.addSrcSinkRxns(model,pmap)

    #print("source sink defined")


    ### REPORT ON THE MODEL HERE OR SEPARATE FUNCTION???
    reptModel(model)

    return model

#def main():
    #__init__

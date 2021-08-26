"""
Author: David Bianchi

A file to add patch changes to the JCVI-Syn3A lipid metabolism constructed from the input .tsv file
"""

# Import the simulation and model building module
import odecell

# Import the Reactions module
#import Rxns

# Constant
volToSA = 1.5874 # Ratio of volume change to surface area change for 200nm radius cell

def addPppParams(model):
    """
    Adds modified parameters from the lipid metabolism to the metabolic model

    Parameters:

    model (model object) - the metabolic model

    Returns:

    None
    """
    #### PDH_E1 Parameters ####
    #model.addParameter('PDH_E1',rxnFormKey='kcatF',value=486)
    #model.addParameter('PDH_E1',rxnFormKey='kcatR',value=1.2763)
    

    #model.addParameter('PDH_E1',rxnFormKey='KmSub1',value=0.1003) #lpl_PdhC
    #model.addParameter('PDH_E1',rxnFormKey='KmSub2',value=0.0537) #pyr
    #model.addParameter('PDH_E1',rxnFormKey='KmProd1',value=0.0997) #acdhlpl_PdhC
    #model.addParameter('PDH_E1',rxnFormKey='KmProd2',value=0.0997) #co2    

    ### GAPDP Parameters ####
    model.addParameter('GAPDP','KmSub2',0.385) # nadp
    model.addParameter('GAPDP','KmProd2',0.202) # nadph
    model.addParameter('GAPDP','kcatF',2.8)
    model.addParameter('GAPDP','kcatR',0)

    ### FMETTRS Parameters ###
    model.addParameter('FMETTRS','kcatF',0.45)

    ### MTHFC Parameters ###
    model.addParameter('MTHFC','kcatF',185)

    #### GHMT2 Paramters ####
    model.addParameter('GHMT2','kcatF',0.0)
    model.addParameter('GHMT2','kcatR',0.0)
    #model.addParameter('PDH_E2',rxnFormKey='kcatF',value=169.4488)
    #model.addParameter('PDH_E2',rxnFormKey='kcatR',value=0.562)

    #model.addParameter('PDH_E2',rxnFormKey='KmSub1',value=0.1003) # acdhlpl_PdhC
    #model.addParameter('PDH_E2',rxnFormKey='KmSub2',value=0.1003) # coa
    #model.addParameter('PDH_E2',rxnFormKey='KmProd1',value=0.0111) # accoa
    #model.addParameter('PDH_E2',rxnFormKey='KmProd2',value=0.0997) # dhlpl_PdhC

    #### TKT1 Parameters ####
    model.addParameter('TKT1',rxnFormKey='kcatF',value=20.58)
    model.addParameter('TKT1',rxnFormKey='kcatR',value=0.8)
    
    model.addParameter('TKT1',rxnFormKey='KmSub1',value=0.743) #g3p
    model.addParameter('TKT1',rxnFormKey='KmSub2',value=3.7298) #s7p
    model.addParameter('TKT1',rxnFormKey='KmProd1',value=0.4717) #r5p
    model.addParameter('TKT1',rxnFormKey='KmProd2',value=0.134) #xu5p
    
    #### TKT2 Parameters ####
    model.addParameter('TKT2',rxnFormKey='kcatF',value=26.87)
    model.addParameter('TKT2',rxnFormKey='kcatR',value=1.4)
    
    model.addParameter('TKT2',rxnFormKey='KmSub1',value=0.25) #f6p
    model.addParameter('TKT2',rxnFormKey='KmSub2',value=0.743) #g3p
    model.addParameter('TKT2',rxnFormKey='KmProd1',value=0.0227) #e4p
    model.addParameter('TKT2',rxnFormKey='KmProd2',value=0.134) #xu5p
    
    #### TALA Parameters ####
    model.addParameter('TALA',rxnFormKey='kcatF',value=22.3)
    model.addParameter('TALA',rxnFormKey='kcatR',value=0.54)
    
    model.addParameter('TALA',rxnFormKey='KmSub1',value=0.0401) #e4p
    model.addParameter('TALA',rxnFormKey='KmSub2',value=0.6688) #f6p
    model.addParameter('TALA',rxnFormKey='KmProd1',value=1.9) #g3p
    model.addParameter('TALA',rxnFormKey='KmProd2',value=0.285) #s7p
    
    #### Speed up DGSN Pathway ####
    model.addParameter('DGSNK',rxnFormKey='kcatF',value=2.25)

    #### Speed up DADN pathway ####
    model.addParameter('PUNP2',rxnFormKey='kcatF',value=13.3)

    #### Speed up FBA rxn ####
    model.addParameter('FBA',rxnFormKey='kcatF',value=64.5)
    #model.addParameter('FBA',rxnFormKey='KmSub1',value=0.17)

    model.addParameter('RNDR2',rxnFormKey='KmSub1',value=0.24)
#    model.addParameter('RNDR3',rxnFormKey='KmSub1',value=0.31)
#   #### RPE Parameters ####
    #model.addParameter('RPE',rxnFormKey='kcatF',value=1.0)
    #model.addParameter('RPE',rxnFormKey='kcatR',value=1.0)
    
    #model.addParameter('RPE',rxnFormKey='KmSub1',value=1.0)
    #model.addParameter('RPE',rxnFormKey='KmProd1',value=1.0)    
    
#     #### RPI Parameters ####
    model.addParameter('RPI',rxnFormKey='kcatF',value=10.0)
    model.addParameter('RPI',rxnFormKey='kcatR',value=1.0)
    
    #model.addParameter('RPI',rxnFormKey='KmSub1',value=1.0)
    #model.addParameter('RPI',rxnFormKey='KmProd1',value=1.0)


#     print('Updated PPP Parameters')

    return



def main(model,pmap):
    """
    Add all lipid related parameters

    Parameters:

        model (model object) - the metabolic model

        pmap (map) - The particle map (LM)

    Returns:
        
        None
    """

    addPppParams(model)

#     addTransportParams(model,pmap)

    #translationSources(model)

    #addLipidMetabs(model)

    return

if __name__ == '__main__':
    # execute only if run as a script
    main()

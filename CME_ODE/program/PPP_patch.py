"""
Author: David Bianchi

A file to add patch changes to the JCVI-Syn3A lipid metabolism constructed from the input .tsv file
"""

# Import the simulation and model building module
import odecell

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

    model.addParameter('RNDR2',rxnFormKey='KmSub1',value=0.24)
    
#     #### RPI Parameters ####
    model.addParameter('RPI',rxnFormKey='kcatF',value=10.0)
    model.addParameter('RPI',rxnFormKey='kcatR',value=1.0)

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

    return

if __name__ == '__main__':
    # execute only if run as a script
    main()

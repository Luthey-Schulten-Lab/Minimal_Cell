"""
Author: David Bianchi, Zane Thornburg

A file to add patch changes to the JCVI-Syn3A lipid metabolism constructed from the input .tsv file
"""

# Import the simulation and model building module
import odecell

# Import the Reactions module
import Rxns_ODE as Rxns

# Constant
volToSA = 1.5874 # Ratio of volume change to surface area change for 200nm radius cell

def addLipidParams(model):
    """
    Adds modified parameters from the lipid metabolism to the metabolic model

    Parameters:

    model (model object) - the metabolic model

    Returns:

    None
    """

    ### ACPPAT: Raise forward kcat, Raise substrate km to literature values
    kcatf_ACPPAT = 20.7 # 1/s kcat propanoyl phosphate - S. enterica EC 2.3.1.222
    km_ap_ACPPAT = 0.0179 # mM - ap - S. enterica (BRENDA)

    model.addParameter('ACPPAT',rxnFormKey='kcatF',parName='kcatF_ACPPAT',unit='1/s',value=kcatf_ACPPAT)
    model.addParameter('ACPPAT',rxnFormKey='KmSub2',parName='km_Sub2_ACPPAT',unit='mM',value=km_ap_ACPPAT)

    ### AGPAT: Fix kcat and km values to literature values for increased flux
    #km_1ag3p_AGPAT = 0.003 # mM S. oleracea (BRENDA)
    kcatf_AGPAT = 144.4 # 1/s, C. moschata (SABIO-RK) - converted from Vmax, Slabas et al

    #model.addParameter('AGPAT',rxnFormKey='KmSub1',parName='km_Sub1_AGPAT',unit='mM',value=km_1ag3p_AGPAT)
    model.addParameter('AGPAT',rxnFormKey='kcatF',parName='kcatF_AGPAT',unit='1/s',value=kcatf_AGPAT)
    model.addParameter('AGPAT',rxnFormKey='kcatR',parName='kcatR_AGPAT',unit='1/s',value=0.0) # Irreversible

    ### PGSA: Fix kcat backwards to zero, this is an irreversible reaction, otherwise CTP issue
    model.addParameter('PGSA',rxnFormKey='kcatR',parName='kcatR_PGSA',unit='1/s',value=0.0) # Irreversible Reaction 

    # Try with irreversible PGPP
    model.addParameter('PGPP',rxnFormKey='kcatR',parName='kcatR_PGPP',unit='1/s',value=0.0)

    ### UDPGALM: galactomutase rate selection from literature
    kcatf_UDPGALM = 27.0 # 1/s E. coli (BRENDA)
    kcatr_UDPGALM = 1.5 # 1/s E. coli (BRENDA)
    kmProd1_UDPGALM = 0.45 # mM Aspergillus (BRENDA)


    model.addParameter('UDPGALM',rxnFormKey='KmProd1',parName='km_Prod1_UDPGALM',unit='mM',value=kmProd1_UDPGALM)
    model.addParameter('UDPGALM',rxnFormKey='kcatR',parName='kcatR_UDPGALM',unit='1/s',value=kcatr_UDPGALM) 
    model.addParameter('UDPGALM',rxnFormKey='kcatF',parName='kcatF_UDPGALM',unit='1/s',value=kcatf_UDPGALM)

    ### DAGGALT: Glycolipid synthesis rate selection from literature
    kcatr_DAGGALT = 0.0 # !/s Irreversible reaction
    #kcatf_DAGGALT = 24.86 # 1/s C. tepidum (BRENDA) #2.09 # 1/s Spinacea olera (BRENDA)

    model.addParameter('DAGGALT',rxnFormKey='kcatR',parName='kcatR_UDPGALM',unit='1/s',value=kcatr_DAGGALT) 
    #model.addParameter('DAGGALT',rxnFormKey='kcatF',parName='kcat1_DAGGALT',unit='1/s',value=kcatf_DAGGALT)


    return

def translationSources(model):
    """
    Source reactions for translated proteins, namely apoACP

    Parameters:

        model (model object): The metabolic model

    Returns:
        
        None

    """

     ### Add non-ATP dependent Fatty acid incorporation
    Kcat_trans_apoACP = 9.99E-04 # mM/s - apoACP _0621 translation flux

    model.addReaction("apoACP_SOURCE","FAFBA","apo-ACP Translation Source")
    model.addProduct("apoACP_SOURCE","Prod1","M_apoACP_c")
    model.addParameter("apoACP_SOURCE","Kcat_R",Kcat_trans_apoACP)

def FAFBA():
    """
    Active Transport rate form.
    
    Example: Bring in Fatty Acid at FBA Flux
                
    Heuristic Rules: Ki = 1/10 Normal Intracellular Concentration (NIC)
                    Ka = 1/25 Normal Intracellular Concentration (NIC)
    
    """
    rateform='($Kcat_R)'#*($Sub1/($Ka_ATP+$Sub1))'
    return rateform

def COAabcR():
    """
    Activate transport rate form for Coenzyme A transport

    Example: coa_e + atp -> coa_c + adp + pi

    Kcat, Km values from Santos et al (2018) - Cobalamine transport

    Enzyme value: Unknown, EcfA - 274, EcfT - 25 (copy number, proteomics), EcfS - Assumed to be 10
    """
    rateform='($KcatR)*($Enzyme)*($Subex/($Subex+$Km1))'
    return rateform

def ActTrans():
    """
    Active Transport rate form.
    
    Example: (Kcat_R_COAabc*enzCOAabc*M_coa_e/(Km_R_COAabc+M_coa_e))*(M_atp_c/(Ka_ATP+M_atp_c))
                *(KiAbc/(KiAbc+M_coa_c))
                
    Heuristic Rules: Ki = 10 Normal Intracellular Concentration (NIC)
                    Ka = 1/25 Normal Intracellular Concentration (NIC)
    
    """

    rateform='($Kcat_R*$Enz_R)'#($Kcat_R*$Enz_R*$Sub1/($Km_R+$Sub1))'#*($Sub2/$Ka_ATP+$Sub2)*($Ki_R/($Ki_R+$Prod1))'
    return rateform

def addTransportParams(model,pmap):
    """
    Adds the transport reaction parameters to the metabolic model object

    Parameters:

        model (model object): the metabolic model

        pmap (map): the particle map (LM)

    Returns:
        
        None
    """

    ### Add the rate forms that match FBA fluxes and account for ATP dependence
    model.FAFBA = odecell.modelbuilder.RateForm(FAFBA())
    model.updateAvailableForms()

    ### Add the rate forms that match FBA fluxes and account for ATP dependence
    model.COAabcR = odecell.modelbuilder.RateForm(COAabcR())
    model.updateAvailableForms()

    Enz_COAabc = Rxns.partTomM(min(pmap['P_641'],pmap['P_642'],pmap['P_643'],pmap['P_836']),pmap) ## But this has to be all COAabc genes (See M. Melo)

    Kcat_R_COAabc = 2.0 # 1/s, Santos et al (2018), CbrT params
    Km_R_COAabc = 2.1e-6 # mM, Santos et al (2018), CbrT params

    M_coa_e = 0.0013 # 1.3 uM (Tables 1&2 def med v3 Table 2)

    ### Define the COAabc reaction
    rxnIndx=model.addReaction("COAabc","COAabcR","Coenzyme A Transport")
    #model.addSubstrate(rxnIndx,"Sub1",'M_coa_e')
    model.addProduct("COAabc","Prod1",'M_coa_c')
    #model.addParameter(rxnIndx,"Prod1",M_coa_c)
    #model.addProduct(rxnIndx,"Prod2","M_atpUsed_c")
    model.addParameter("COAabc","Subex",M_coa_e)
    model.addSubstrate("COAabc","Sub2",'M_atp_c') # Update RL
    model.addProduct("COAabc","Prod2",'M_adp_c') # Update RL
    model.addParameter("COAabc","Prod3",'M_pi_c') # Update RL
    #model.addParameter(rxnIndx,"Kcat_R",7.0e-6)
    model.addParameter("COAabc","KcatR",Kcat_R_COAabc)
    model.addParameter("COAabc","Km1",Km_R_COAabc)
    model.addParameter("COAabc","Enzyme",Enz_COAabc)
    #model.addParameter(rxnIndx,"Ki_R",Ki_R_COAabc)
    #model.addParameter(rxnIndx,"Ka_ATP",Ka_R_COAabc)

    ### Add Sphingomyelin uptake
    model.addMetabolite('M_sm_c','Sphingomyelin',initVal=Rxns.partTomM(pmap['M_sm_c'],pmap))
    Kcat_FBA_SM = 1.10e-3*(volToSA/2.0)#*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA

    rxnIndx=model.addReaction("SMuptake","FAFBA","Sphingomyelin Uptake")
    model.addProduct("SMuptake","Prod1","M_sm_c")
    model.addParameter("SMuptake","Kcat_R",Kcat_FBA_SM)

    # Add PC to the model
    
    model.addMetabolite('M_pc_c','Phosphatidylcholine',initVal=Rxns.partTomM(pmap['M_pc_c'],pmap))
    Kcat_FBA_PC = 1.42e-4*(volToSA/2.0)#*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA
    rxnIndx=model.addReaction("PCuptake","FAFBA","Phosphatidylcholine Uptake")
    model.addProduct("PCuptake","Prod1","M_pc_c")
    model.addProduct("PCuptake","Prod2","M_adp_c")
    model.addProduct("PCuptake","Prod3","M_pi_c")
    model.addSubstrate("PCuptake","Sub1","M_atp_c")
    model.addParameter("PCuptake","Kcat_R",Kcat_FBA_PC)

    ### Add Triacylglycerol uptake
    model.addMetabolite('M_tag_c','triacylglycerol',initVal=Rxns.partTomM(pmap['M_tag_c'],pmap))
    Kcat_FBA_TAG = 0.0003095*(volToSA/2.0)#*volToSA # mM/s, Uptake Flux Triacylglycerol from FBA

    rxnIndx = model.addReaction("TAGt","FAFBA","TAG uptake")
    model.addProduct("TAGt","Prod1","M_tag_c")
    model.addParameter("TAGt","Kcat_R",Kcat_FBA_TAG)


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

    addLipidParams(model)

    addTransportParams(model,pmap)

    return

if __name__ == '__main__':
    # execute only if run as a script
    main()

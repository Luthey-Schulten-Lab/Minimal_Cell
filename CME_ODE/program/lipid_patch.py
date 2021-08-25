"""
Author: David Bianchi

A file to add patch changes to the JCVI-Syn3A lipid metabolism constructed from the input .tsv file
"""

# Import the simulation and model building module
import odecell

# Import the Reactions module
import Rxns_two as Rxns

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

    ### GLYK: Lower glycerol kinase kcatfwd than from Parameter Balancing
    #model.addParameter('GLYK',rxnFormKey='kcat1',parName='kcat1_GLYK',unit='1/s',value=1670.0) # BRENDA: Candida Mycoderma

    #model.addMetabolite('M_atpUsed_GLYK_c','ATP used by lipid rxns',initVal=0.0)
    #model.addMetabolite('M_adpMade_GLYK_c','ADP made by GLYK',initVal=0.0)
    ### Necessary or pre-defined in lipid/central
    #model.addProduct("GLYK","Prod3","M_atpUsed_GLYK_c")
    #model.addProduct("GLYK","Prod4","M_adpMade_GLYK_c")

    #model.addMetabolite('M_adpMade_FAKr_c','ADP made by FAKr',initVal=0.0)
    #model.addProduct('FAKr','Prod3','M_adpMade_FAKr_c')
    #kcatf_FAKr = 38.0 # 1/s kcat fa - S. aureus, Rock et al 2016
    #model.addParameter('FAKr',rxnFormKey='kcatF',parName='kcatF_FAKr',unit='1/s',value=kcatf_FAKr)
    #model.addMetabolite('M_coaUsed_

    ### Add Phosphoesterase Enzyme Level - GAP FILL in the metabolic module
    #enz_PAPA = 86*5e-5 # 86 from proteomics gene: JCVISYN3A_
    #model.addParameter('PAPA',rxnFormKey='Enzyme',parName='enzLvl_PAPA',unit='mM',value=enz_PAPA) 

    ### ACPPAT: Raise forward kcat, Raise substrate km to literature values
    kcatf_ACPPAT = 20.7 # 1/s kcat propanoyl phosphate - S. enterica EC 2.3.1.222
    km_ap_ACPPAT = 0.0179 # mM - ap - S. enterica (BRENDA)

    model.addParameter('ACPPAT',rxnFormKey='kcatF',parName='kcatF_ACPPAT',unit='1/s',value=kcatf_ACPPAT)
    model.addParameter('ACPPAT',rxnFormKey='KmSub2',parName='km_Sub2_ACPPAT',unit='mM',value=km_ap_ACPPAT)

    ### DASYN: Fix DASYN enzyme level (JCVISYN3A_0304), was 1 from proteomics 
    #enzLvl_DASYN = 40*5e-5 # 40 proteomics
    #model.addParameter('DASYN',rxnFormKey='Enzyme',parName='enzLvl_DASYN',unit='mM',value=enzLvl_DASYN)

    ### APG3PAT: Fix AGPAT enzyme level (JCVISYN3A_0117), was 1 from proteomics
    #enzLvl_APG3PAT = 100*5e-5 # mM _0117
    #model.addParameter('APG3PAT',rxnFormKey='Enzyme',parName='enzLvl_APG3PAT',unit='mM',value=enzLvl_APG3PAT)

    ### AGPAT: Fix kcat and km values to literature values for increased flux
    #km_1ag3p_AGPAT = 0.003 # mM S. oleracea (BRENDA)
    kcatf_AGPAT = 144.4 # 1/s, C. moschata (SABIO-RK) - converted from Vmax, Slabas et al

    #model.addParameter('AGPAT',rxnFormKey='KmSub1',parName='km_Sub1_AGPAT',unit='mM',value=km_1ag3p_AGPAT)
    model.addParameter('AGPAT',rxnFormKey='kcatF',parName='kcatF_AGPAT',unit='1/s',value=kcatf_AGPAT)
    model.addParameter('AGPAT',rxnFormKey='kcatR',parName='kcatR_AGPAT',unit='1/s',value=0.0) # Irreversible

    ### PGSA: Fix kcat backwards to zero, this is an irreversible reaction, otherwise CTP issue
    model.addParameter('PGSA',rxnFormKey='kcatR',parName='kcatR_PGSA',unit='1/s',value=0.0) # Irreversible Reaction 

    ### CLPNS: Fix kcat forward to the lower available kcat value
    #CLPNS_kcat_eco = 3.6 # 1/s, E. coli (BRENDA), Hiraoka et al J. Biochem. (1991)
    #CLPNS_kcat_r_hal = 24.79 # 1/s, From Haldane Relationship (Main Text Leibermeister 2010, Eq. 15)
    # Try with irreversible PGPP
    model.addParameter('PGPP',rxnFormKey='kcatR',parName='kcatR_PGPP',unit='1/s',value=0.0)
    #model.addParameter('CLPNS',rxnFormKey='kcatF',parName='kcatF_CLPNS',unit='1/s',value=CLPNS_kcat_eco)
    ### PGMT: phosphoglucomutase reverse reaction Km from literature
    #km_g1p_PGMT = 1.1 # mM, glucose 1-phosphate - E. coli (BRENDA) - Fujimoto et al Biochim. Biophys. Acta (1965)
    #kcat_fwd_PGMT = 6.0 # 1/s - fitting #0.668 # 1/s, M. smegmatis (BRENDA) - Fujimoto et al Biochim. Biophys. Acta (1965)
    #kcat_rev_PGMT = 2.5 # 1/s - Vmaxs spreadsheet
    #print(kcat_fwd_PGMT)
    #model.addParameter('PGMT',rxnFormKey='kcatR',parName='kcatR_PGMT',unit='1/s',value=kcat_rev_PGMT)
    #model.addParameter('PGMT',rxnFormKey='KmProd1',parName='km_Prod1_PGMT',unit='mM',value=km_g1p_PGMT)

    ### BPNT: Bisphosphate nucleotidase -> AMP production
    #kcat_fwd_BPNT = 5.4 # 1/s, BRENDA M. tuberculosis - Hatzios et al Biochemistry (2008) 
    #model.addParameter('BPNT',rxnFormKey='kcatF',parName='kcatF_BPNT',unit='1/s',value=kcat_fwd_BPNT) 
    #model.addMetabolite('M_ampMade_BPNT_c','AMP made by BPNT',initVal=0.0)
    #model.addMetabolite('M_piMade_BPNT_c','PI made by BPNT',initVal=0.0)
    #model.addProduct("BPNT","Prod3","M_ampMade_BPNT_c")
    #model.addProduct("BPNT","Prod4","M_piMade_BPNT_c")


    ### UDPG4E: epimerase rate selection from literature
    #kmSub1_UDPG4E = 0.76 # mM # Camplyobacter (BRENDA)  #29.9 # mM E. coli (BRENDA)
    #kmProd1_UDPG4E = 0.18 # mM # E. coli (BRENDA) #0.1181 # mM M. Maripaludis (Queen's U. 2018 - Sharma et al Thesis)
    #kcat2_UDPG4E = 12.9 # 1/s Aeromonas hydrophila (BRENDA)
    #kcat1_UDPG4E = 74.0 # 1/s Thermatoga (BRENDA) - EC 5.1.3.2

    #model.addParameter('UDPG4E',rxnFormKey='KmSub1',parName='km_Sub1_UDPG4E',unit='mM',value=kmSub1_UDPG4E) # E. coli (BRENDA)
    #model.addParameter('UDPG4E',rxnFormKey='KmProd1',parName='km_Prod1_UDPG4E',unit='mM',value=kmProd1_UDPG4E)
    #model.addParameter('UDPG4E',rxnFormKey='kcatR',parName='kcatF_UDPG4E',unit='1/s',value=kcat2_UDPG4E) 
    #model.addParameter('UDPG4E',rxnFormKey='kcatF',parName='kcatR_UDPG4E',unit='1/s',value=kcat1_UDPG4E) # Thermatoga (BRENDA) EC 5.1.3.2

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

    ### Add non-ATP dependent Fatty acid incorporation
    #Kcat_FBA_FA = 0.007877588767506535*(volToSA/2.0)#*volToSA#/VolToSA # mM/s - FBA Flux

    #rxnIndx=model.addReaction("FAt","FAFBA","Fatty Acid Transport")
    #model.addProduct("FAt","Prod1","M_fa_c")
    #model.addParameter("FAt","Kcat_R",Kcat_FBA_FA)

    ### Add the rate forms that match FBA fluxes and account for ATP dependence
    model.COAabcR = odecell.modelbuilder.RateForm(COAabcR())
    model.updateAvailableForms()

    Enz_COAabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0641_c'],pmap['M_PTN_JCVISYN3A_0642_c'],pmap['M_PTN_JCVISYN3A_0643_c'],pmap['M_PTN_JCVISYN3A_0836_c']),pmap) ## But this has to be all COAabc genes (See M. Melo)

    Kcat_R_COAabc = 2.0 # 1/s, Santos et al (2018), CbrT params
    Km_R_COAabc = 2.1e-6 # mM, Santos et al (2018), CbrT params
    #Ka_R_COAabc = 3.6/25.0 # mM, Heuristic Rule Atlas et al (2008), NIC: 3.63 mM Rabinowitz et al 2018
    #Ki_R_COAabc = 10*0.8 # mM, Heuristic Rule Atlas et al (2008), NIC: 0.8 mM Rabinowitz et al 2018

    M_coa_e = 0.0013 # 1.3 uM (Tables 1&2 def med v3 Table 2)
    #M_coa_e = 0.0029 # mM (C5Mod-CMRL)

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

    ### Add ATP dependent COAabc uptake
    #Kcat_FBA_COAabc = 5.182980193627005e-06 # FBA Flux mM/s
    #Kcat_FBA_COAabc = 6.48e-6
    #model.addMetabolite('M_atpUsed_COAabc_c','ATP used by COA transport',initVal=0.0)


    ### Kcat, Km and Enzyme go here
    #Kcat_R_COAabc = 21.1 #7.5e-4 #21.1 # 1/s (21.1 from Cobalamine conversion - doesn't require volume Shuler-Active.ipynb)\n",
    #M_coa_e = 0.0013 # mM (K. wise defined media)
    #Km_R_COAabc = 2.1e-6 # mM 2.1e-9 M or 2.1e-6 mM (Santos et al - 2018)\n",
    #Enz_COAabc = 10*(5e-5) # mM, corresponds to ~10 enzymes\n",
    #"Ki_R_COAabc = 2*0.8*10.0 #0.12*10.0 # mM - 10 NIC CoA (From pnuemonia:0.12 mM,From Rabinowitz et al 2016: 1.37 mM) \n",


    #rxnIndx=model.addReaction("COAabc","COAabc","Coenzyme A Transport")
    #model.addProduct("COAabc","Prod1","M_coa_c")
    #model.addProduct("COAabc","Prod2","M_atpUsed_COAabc_c")
    #model.addProduct("COAabc","Prod3","M_pi_c")
    #model.addProduct("COAabc","Prod4","M_adp_c")
    #model.addSubstrate("COAabc","Sub1","M_atp_c")
    #model.addParameter("COAabc","Kcat_R",Kcat_R_COAabc)
    #model.addParameter("COAabc","Km1",Km_R_COAabc)
    #model.addParameter("COAabc","Subex",M_coa_e)
    #model.addParameter("COAabc","Enzyme",Enz_COAabc)

    ### Add Sphingomyelin uptake
    model.addMetabolite('M_sm_c','Sphingomyelin',initVal=Rxns.partTomM(pmap['M_sm_c'],pmap))
    Kcat_FBA_SM = 1.10e-3*(volToSA/2.0)*1.25#*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA

    rxnIndx=model.addReaction("SMuptake","FAFBA","Sphingomyelin Uptake")
    model.addProduct("SMuptake","Prod1","M_sm_c")
    model.addParameter("SMuptake","Kcat_R",Kcat_FBA_SM)

    # Add PC to the model
    
    model.addMetabolite('M_pc_c','Phosphatidylcholine',initVal=Rxns.partTomM(pmap['M_pc_c'],pmap))
    Kcat_FBA_PC = 1.42e-4*(volToSA/2.0)*1.25#*volToSA#/VolToSA # mM/s, Uptake flux from lipidomics adjusted FBA
    rxnIndx=model.addReaction("PCuptake","FAFBA","Phosphatidylcholine Uptake")
    model.addProduct("PCuptake","Prod1","M_pc_c")
    model.addProduct("PCuptake","Prod2","M_adp_c")
    model.addProduct("PCuptake","Prod3","M_pi_c")
    model.addSubstrate("PCuptake","Sub1","M_atp_c")
    model.addParameter("PCuptake","Kcat_R",Kcat_FBA_PC)


    ### Add Cholesterol uptake
    #model.addMetabolite('M_chsterol_c','cholesterol',initVal=0.0)
    #Kcat_FBA_Chol = 0.003604 # mM/s, Uptake flux Cholesterol from lipidomics
    
    #rxnIndx = model.addReaction("CHOLt","FAFBA","Cholesterol Uptake")
    #model.addProduct(rxnIndx,"Prod1","M_chsterol_c")
    #model.addParameter(rxnIndx,"Kcat_R",Kcat_FBA_Chol)


    ### Add Triacylglycerol uptake
    model.addMetabolite('M_tag_c','triacylglycerol',initVal=Rxns.partTomM(pmap['M_tag_c'],pmap))
    Kcat_FBA_TAG = 0.0003095*(volToSA/2.0)*1.25#*volToSA # mM/s, Uptake Flux Triacylglycerol from FBA

    rxnIndx = model.addReaction("TAGt","FAFBA","TAG uptake")
    model.addProduct("TAGt","Prod1","M_tag_c")
    model.addParameter("TAGt","Kcat_R",Kcat_FBA_TAG)

    ### FOR ADDING GLYCt
    # P_R_GLYCt = 4.50E-10 #m/s #Check SA/V
    #Kcat_FBA_Glyc = 0.00595 # mM/s FBA Flux

    #rxnIndx = model.addReaction("GLYCt","FAFBA","Glycerol Uptake")
    #model.addProduct(rxnIndx,"Prod1","M_glyc_c")
    #model.addParameter(rxnIndx,"Kcat_R",Kcat_FBA_Glyc)

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

    #translationSources(model)

    #addLipidMetabs(model)

    return

if __name__ == '__main__':
    # execute only if run as a script
    main()

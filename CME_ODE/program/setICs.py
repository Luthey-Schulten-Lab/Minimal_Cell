"""
Sets the initial conditions for metabolite values in the Lattice Microbes CME particle map.
    Initialize specialized enzyme concentrations that have erroneous Mass Spec. Values (1 etc.)

Reads all initial concentration conditions from the !Quantity table of metabolic module .tsv files and 
supplements with additional values

Author: David Bianchi
Date: 8/24/20
"""

import numpy as np

# SBtab classes source code
from sbtab import SBtab

# SBtab validator
from sbtab import validatorSBtab


NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

def mMtoPart(mM):
    particle = (mM/1000)*NA*V
    return particle


def ValidateandSeparate(fn):
    """
    Validate tsv read in file

    Parameters:
        
        fn (string): filename containing concentration parameters

    Returns: 
        
        ComDF_sub (pd.Dataframe): Pandas dataframe containing compound information
    """

        # open a file and read it
    with open(fn,"r") as infile:
        file_content = infile.read()

    # create an SBtab Document Object Sd
    Sd = SBtab.SBtabDocument('BalancedModel', file_content, fn)

    validator = validatorSBtab.ValidateDocument(Sd)
    # print("warnings:",validator.validate_document())
    warnings = validator.validate_document()
    for warn in warnings:
        if warn[1]:
            print(warn[1][0])


    # Read SBTab tables that will be used in the kinetic model
    SBtabQnt = Sd.get_sbtab_by_name("Quantity")
    SBtabCom = Sd.get_sbtab_by_name("Compound")

    # Transform SBtab tables into Pandas dataframes
    # RxnDF = SBtabRxn.to_data_frame()
    QntDF_sub = SBtabQnt.to_data_frame()
    ComDF_sub = SBtabCom.to_data_frame()
    
    return ComDF_sub
    

def initNucConcs(sim):
    """
    Initialize nucleotide concentrations

    Parameters:

        sim (lm:cme:simulation): The lattice microbes CME simulation object

    Returns:

        None
    """

    file_name = "../model_data/Nucleotide_Kinetic_Parameters.tsv"

    ComDF_nuc = ValidateandSeparate(file_name)

    nuclMets = ["cytd_c","cmp_c","cdp_c","ctp_c","dcyt_c","dcmp_c","dcdp_c","dctp_c",
           "gsn_c","gua_c","gmp_c","gdp_c","gtp_c","dgsn_c","dgmp_c","dgdp_c","dgtp_c",
           "ura_c","uri_c","ump_c","udp_c","utp_c","duri_c","dump_c","dudp_c","dutp_c",
           "adn_c","ade_c","dad_2_c","damp_c","dadp_c","datp_c","thymd_c","dtmp_c","dtdp_c","dttp_c",
           "glu__L_c","gln__L_c","ppi_c"]#,"trdrd_c"]


# concentrationList = []

    for met in nuclMets:
        metID = str("M_" + met)
        conc = ComDF_nuc.loc[ ComDF_nuc["Compound"] == metID, "InitialConcentration" ].values[0]
        parts = int(round(mMtoPart(float(conc))))
        sim.defineSpecies([metID])
        sim.addParticles(metID,count=parts)

    return None


def initCentConcs(sim):
    """
    Initialize central metabolism concentrations

    Parameters:

        sim (lm:cme:simulation): The lattice microbes CME simulation object

    Returns:

        None
    """

    file_name = "../model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"

    ComDF_cent = ValidateandSeparate(file_name)

    centralMets = ["g6p_c","g1p_c","f6p_c","man6p_c","gam6p_c","acgam6p_c","acmanap_c","pi_c","amp_c","adp_c","atp_c","acmana_c","thfglu3_c","10fthfglu3_c","nh3_c","pep_c","pyr_c","xu5p__D_c","ru5p__D_c","r5p_c","prpp_c","e4p_c","r1p_c","fdp_c", "s7p_c","dhap_c","g3p_c","2dr5p_c","2dr1p_c","nad_c","nadh_c","13dpg_c","nadp_c","nadph_c","3pg_c","2pg_c","lac__L_c","co2_c","accoa_c","actp_c","ac_c","acald_c","o2_c",
"nac_c","nicrnt_c","dnad_c","ribflv_c","fmn_c","fad_c","pydx5p_c","thmpp_c","5fthf_c","5fthfglu3_c","methfglu3_c","mlthfglu3_c","sprm_c"]
# concentrationList = []

    for met in centralMets:
        metID = "M_" + met
        print(metID)
        conc = ComDF_cent.loc[ ComDF_cent["Compound"] == metID, "InitialConcentration" ].values[0]
        parts = int(round(mMtoPart(float(conc))))
        print(metID,parts)
        sim.defineSpecies([metID])
        sim.addParticles(metID,count=parts)

#     sim.defineSpecies(["M_lpl_PdhC_c"])
#     sim.addParticles("M_lpl_PdhC_c",int(round(mMtoPart(0.009024644220757/2))))

#     sim.defineSpecies(["M_dhlpl_PdhC_c"])
#     sim.addParticles("M_dhlpl_PdhC_c",int(round(mMtoPart(0.009024644220757/2))))

#     sim.defineSpecies(["M_acdhlpl_PdhC_c"])
#     sim.addParticles("M_acdhlpl_PdhC_c",int(round(mMtoPart(0.0))))    

    return None


def initLipConcs(sim):
    """ 
    Initialize central metabolism concentrations

    Parameters:

        sim (lm:cme:simulation): The lattice microbes CME simulation object

    eturns:

        None
    """

    file_name = "../model_data/lipid_NoH2O_balanced_model.tsv"

    ComDF_lip = ValidateandSeparate(file_name)

    lipMets = ["udpg_c","udpgal_c","ap_c","pa_c","glyc_c","glyc3p_c","coa_c","pap_c","1ag3p_c","cdpdag_c","pg3p_c","pg_c","clpn_c","12dgr_c","galfur12dgr_c",
          "udpgalfur_c","fa_c"]#,"ACP_c","ACP_R_c"] # "fa_c" #"chsterol_c", ACP_R_c, ACP_c

    for met in lipMets:
        metID = "M_" + met
        conc = ComDF_lip.loc[ ComDF_lip["Compound"] == metID, "InitialConcentration" ].values[0]
        
        sim.defineSpecies([metID])
        
        parts = int(round(mMtoPart(float(conc))))

        sim.addParticles(metID,count=parts)

    lipDict_additional ={"M_tag_c":1.95,
             "M_pc_c":1.20,
             "M_sm_c":9.22,
             "M_chsterol_c":23.29,
             'M_atpUsed_GLYK_c':0.0,
             'M_adpMade_GLYK_c':0.0,
             'M_adpMade_FAKr_c':0.0,
             'M_ampMade_BPNT_c':0.0,
             'M_piMade_BPNT_c':0.0}

    for key, val in lipDict_additional.items():
        sim.defineSpecies([key])

    for key, val in lipDict_additional.items():
        val = int(round(mMtoPart(val))) 
        sim.addParticles(key,count=val)
    

    return None

def initAAConcs(sim):
    """ 
    Initialize central metabolism concentrations

    Parameters:

        sim (lm:cme:simulation): The lattice microbes CME simulation object

    Returns:

        None
    """

    file_name = "../model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"

    ComDF_aa = ValidateandSeparate(file_name)

    aaMetIDs = ["M_ala__L_c", "M_arg__L_c",
    "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_gly_c",
    "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
    "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_amet_c"]


    for metID in aaMetIDs:
        conc = ComDF_aa.loc[ ComDF_aa["Compound"] == metID, "InitialConcentration" ].values[0]
        parts = int(round(mMtoPart(float(conc))))
        sim.defineSpecies([metID])
        sim.addParticles(metID,count=parts)

    return None


def initTRNAConcs(sim):
    """
    Initialize the tRNA concentration values

    Follow what is in there for the tRNA concentrations already
    
    Parameters:
        sim (cme:lm:sim) - the lattice microbes CME simulation

    Returns:
        None
    """

#     charged_trna = ["M_alatrna_c", "M_argtrna_c",
#         "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c",
#         "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c",
#         "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c"]

#     uncharged_trna = ["M_trnaala_c", "M_trnaarg_c",
#         "M_trnaasn_c", "M_trnaasp_c", "M_trnacys_c", "M_trnaglu_c", "M_trnagln_c", "M_trnagly_c",
#         "M_trnahis_c", "M_trnaile_c", "M_trnaleu_c", "M_trnalys_c", "M_trnamet_c", "M_trnaphe_c",
#         "M_trnapro_c", "M_trnaser_c", "M_trnathr_c", "M_trnatrp_c", "M_trnatyr_c", "M_trnaval_c"]

#     for met in charged_trna:
#         metPartCounts = int(round(7000/20*0.8))
#         sim.defineSpecies([met])
#         sim.addParticles(met,count=metPartCounts)

#     for met in uncharged_trna:
#         metPartCounts = int(round(3750/20*0.2))
#         sim.defineSpecies([met])
#         sim.addParticles(met,count=metPartCounts)

    sim.defineSpecies(['M_fmettrna_c'])
    sim.addParticles('M_fmettrna_c',10)

    sim.defineSpecies(['M_glutrnagln_c'])
    sim.addParticles('M_glutrnagln_c',10)    

    return None

def initTransportConcs(sim):
    """
    Initialize the transport metabolite values from the notebook
    """

    transport_Dict = {"M_ac_e":0.01,"M_pyr_e":0.0,"M_lac__L_e":0.0,"M_h_e":0.0,
    "M_cytd_e":0.0,"M_ura_e":0.0,"M_adn_e":0.0,"M_dad_2_e":0.0,"M_gsn_e":0.0,
    "M_dgsn_e":0.0,"M_dcyt_e":0.0,"M_duri_e":0.0,"M_thymd_e":0.0,"M_uri_e":0.0,#"M_co2_c":0.0234,"M_o2_c":0.1,
    "M_glc__D_e":40,"M_coa_e":0.0013,"M_glyc_e":5.0,"M_fa_e":1.4,
    'M_chsterol_e':0.005,'M_nac_e':0.0016,'M_arg__L_e':3.32,'M_ala__L_e':2.81,
    'M_asp__L_e':2.25,'M_asn__L_e':0.0,'M_cys__L_e':16.5,'M_glu__L_e':5.1,
    'M_his__L_e':0.95,'M_ile__L_e':1.53,'M_leu__L_e':4.58,'M_lys__L_e':3.83,
    'M_met__L_e':1.01,'M_phe__L_e':1.52,'M_pro__L_e':3.48,'M_ser__L_e':2.38,
    'M_thr__L_e':2.52,'M_trp__L_e':0.49,'M_tyr__L_e':2.21,'M_val__L_e':2.14,
    'M_gly_e':6.67}
#      "ptsi":0.01750*0.05, "ptsi_P":0.01750*0.95,"ptsh":0.01438*0.05,
#     "ptsh_P":0.01438*0.95, "crr":0.01557*0.05,"crr_P":0.01557*0.95,"ptsg":0.04121*0.15,
#     "ptsg_P":0.04121*0.85,}

    for key, val in transport_Dict.items():
        sim.defineSpecies([key])

    for key, val in transport_Dict.items():
        val = int(round(mMtoPart(val)))
        sim.addParticles(key,count=val)

    return None

def initEnzConcs(sim):
    """
    Initialize specialized enzyme concentrations that have erroneous Mass Spec. Values (1 etc.)
    """
    
    # AG3PAT Original MS proteomics value was 1, set to 100.0
    #pmap['M_PTN_JCVISYN3A_0117']=100.0
    #sim.defineSpecies(['M_PTN_JCVISYN3A_0117_c'])
    #sim.addParticles('M_PTN_JCVISYN3A_0117_c',100)

    # DASYN Orignal MS proteomics value was 1, set to default 10.0 
    #pmap['M_PTN_JCVISNY3A_0304']=10.0
    #sim.defineSpecies(['M_PTN_JCVISYN3A_0304_c'])
    #sim.addParticles('M_PTN_JCVISYN3A_0304_c',40)    

    # ACPS concentration starts at 1 (MS Proteomics), set to 50
    #sim.defineSpecies(['M_PTN_JCVISYN3A_0513_c'])
    #sim.addParticles('M_PTN_JCVISYN3A_0513_c',50)

    return
    
def __main__(sim):
    """
    Add the concentrations for each subsystem

    Parameters:

        sim (lm:cme:simulation): The lattice microbes CME simulation object

    Returns:

        None
    """


    initLipConcs(sim)
    initNucConcs(sim)
    initCentConcs(sim)
    initAAConcs(sim)
    initTRNAConcs(sim)
    initTransportConcs(sim)
    #initEnzConcs(sim)    

    sim.defineSpecies(['CellSA'])
    #sim.addParticles('CellSA',int(round(304311+198189))) # Case where fa surface area accounted for (40% prot.)
    #sim.addParticles('CellSA',int(round(266605+235895))) # Case where fa surface area not accounted for (47% prot.)
    #sim.addParticles('CellSA',int(round(235984+266515))) # Case with lower cholesterol SA (53% prot)
    sim.addParticles('CellSA',int(round(231875+270956))) # Rescale FA, 54% protein
    #sim.addParticles('CellSA',int(round(4*np.pi*(200**2))))
    

    sim.defineSpecies(['M_k_c'])
    sim.addParticles('M_k_c',int(mMtoPart(10)))
    
    sim.defineSpecies(['M_ca2_c'])
    sim.addParticles('M_ca2_c',int(mMtoPart(1.41)))
    
    sim.defineSpecies(['M_mg2_c'])
    sim.addParticles('M_mg2_c',int(mMtoPart(2.36)))
    
    sim.defineSpecies(['M_glc__D_c'])
    sim.addParticles('M_glc__D_c',int(mMtoPart(1.0)))
    
    sim.defineSpecies(['CellSA_Lip'])
    #sim.addParticles('CellSA_Lip',int(round(304311))) # Case where fa surface area accounted for (40% prot.)
    #sim.addParticles('CellSA_Lip',int(round(266605))) # Case where fa surface area not accounted for (47% prot.)
    #sim.addParticles('CellSA_Lip',int(round(235984))) # Case with lower cholesterol SA (53% prot)
   # sim.addParticles('CellSA_Lip',int(round(262988))) # 11/10 rescale, 53% Lipid
    sim.addParticles('CellSA_Lip',int(round(231875))) # Rescale FA, 46% Lipid

    # Add a cell surface area covered by protein species, equal to amount covered by prot from SA estimate
    sim.defineSpecies(['CellSA_Prot'])   
    #sim.addParticles('CellSA_Prot',int(round(198189))) # Case where fa surface area accounted for (40% prot.)
    #sim.addParticles('CellSA_Prot',int(round(235895))) # Case where fa surface area not accounted for (47% prot.)
    #sim.addParticles('CellSA_Prot',int(round(266515))) # Case with lower cholesterol SA (53%)
    #sim.addParticles('CellSA_Prot',int(round(239512))) # 11/10 rescale, 47% protein
    sim.addParticles('CellSA_Prot',int(round(270956))) # Rescale FA, 54% protein 
    sim.defineSpecies(['CellV'])

    sim.addParticles('CellV',int(round(335))) # Cell Volume x 10^19 L


    return None

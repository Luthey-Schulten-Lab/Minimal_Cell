"""
A file defining Lipid, Central, Nucleotide and Transport Module Metabolic Reactions for JCVI Syn3A

Author: David Bianchi
"""

import Rxns as Rxns
import importlib
import odecell
import odecell.modelbuilder as modelbuilder
importlib.reload(odecell)

import libsbml
import csv
import pandas as pd
import numpy as np
import re

from collections import defaultdict, OrderedDict

import sys

# SBtab classes source code
from sbtab import SBtab

# SBtab validator
from sbtab import validatorSBtab

# Converter SBtab -> SBML
from sbtab import sbtab2sbml

# Converter SBML -> SBtab
from sbtab import sbml2sbtab


### CONSTANTS

### Constants
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

countToMiliMol = 1000/(NA*V)


# The reconstruction matches reactions with gene-protein-reactions (GPR) that use MMSYN1* IDs.
reconstPD = pd.read_excel("../model_data/reconstruction.xlsx", sheet_name='Reactions') # _DB or no

# The annotation matches MMSYN1* IDs with JCVISYN3* IDs (or "locus tags"). # _DB or no
annotatPD = pd.read_excel("../model_data/FBA/Syn3A_annotation_compilation.xlsx",
                         sheet_name="Syn3A_annotation_compilation_condensed")

# The genome data matches "locus tags" with AOE* protein IDs.
# It provides both the gene sequence, needed for transcription reactions in the ODE model,
# and the protein sequence, needed for translation reactions in the model.
# This is the NCBI Gene Bank-formated file Locus "CP016816" (https://www.ncbi.nlm.nih.gov/nuccore/CP016816.1)
# genomeFile3A = '../model_data/syn3A.gb'
# genome3A = next(SeqIO.parse(genomeFile3A, "gb"))

# The proteomics matches AOE IDs with quantitative proteomics data.
proteomPD = pd.read_excel("../model_data/proteomics.xlsx", sheet_name="Proteomics", skiprows=[0] )

# The manual GPR conversion matches proteins with known abundances (AOE IDs), with reactions that did not 
#   have associated genes in the Syn3A reconstruction.
manGPRPD = pd.read_csv("../model_data/manual_GPR_conversion.csv", header=None, names=["MM","AOE"])


sbmlFile = "../model_data/iMB155_NoH2O.xml"

docSBML = libsbml.readSBMLFromFile(sbmlFile)
modelSBML = docSBML.getModel()

speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
speciesNamesLower = [x.lower() for x in speciesNames]
speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]


# # Read global parameters CSV file
# globalParsDF = pd.read_csv("../model_data/GlobalParameters.csv")


# Default metabolite concentration (mM)
defaultMetConcentration = 0.1

# Default Protein concentration (mM): set to 1 micro-Mol 
defaultPtnConcentration = defaultMetConcentration/100


#######################################

##### Import Central Metabolism  ######

#######################################


file_name = "../model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"

# open a file and read it
with open(file_name,"r") as infile:
    file_content = infile.read()

# create an SBtab Document Object Sd
Sd = SBtab.SBtabDocument('BalancedModel', file_content, file_name)

validator = validatorSBtab.ValidateDocument(Sd)
# print("warnings:",validator.validate_document())
warnings = validator.validate_document()
for warn in warnings:
    if warn[1]:
        print(warn[1][0])
        
        
# Read SBTab tables that will be used in the kinetic model
SBtabRxnCent = Sd.get_sbtab_by_name("Reaction")
SBtabParCent = Sd.get_sbtab_by_name("Parameter")
SBtabQntCent = Sd.get_sbtab_by_name("Quantity")
SBtabComCent = Sd.get_sbtab_by_name("Compound")

# Transform SBtab tables into Pandas dataframes
RxnDF = SBtabRxnCent.to_data_frame()
ParDF = SBtabParCent.to_data_frame()
QntDF = SBtabQntCent.to_data_frame()
CentQntDF = SBtabQntCent.to_data_frame()
ComDF = SBtabComCent.to_data_frame()

# Initialize dictionary for "quantity" values
QntDFDict = dict()


# Reaction IDs for central metabolism reactions
cntrMetRxn = ["PGI","PFK","FBA","TPI","GAPD","PGK","PGM","ENO","PYK",
              "LDH_L","PDH_acald","PDH_E3","PTAr","ACKr","NOX","TALA","TKT1",
              "TKT2","RPE","RPI","PRPPS","PPM","PPM2","DRPA","GAPDP"] #"PDH_E1","PDH_E2" #"PDH_E1","PDH_E2", "GAPDP", 
              # Not in the metabolic model: "MAN6PI",,"AMANK","AMANPEr","AGDC","G6PDA",

cntrMetRxnID = [ "R_" + rxn for rxn in cntrMetRxn ]

RxnDF_CentralMet = RxnDF.loc[ RxnDF.Reaction.isin(cntrMetRxnID) ]


selectArray = np.ndarray(QntDF.shape[0],dtype=np.bool)
selectArray.fill(False)

# Turns out it isn't very easy to select rows based on a list of substrings that may appear in 
#  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
for rxnIDR in cntrMetRxnID:
    selectArray = np.logical_or(selectArray, QntDF['Quantity'].str.contains(rxnIDR).values)

QntDF = QntDF.loc[selectArray]


# Read in the parameters in the parameter dataframe for reactions in Central metabolism.
# This gives us the forward and reverse kcats to use later in kinetic rates.

selectArray = np.ndarray(ParDF.shape[0],dtype=np.bool)
selectArray.fill(False)

for rxnIDR in cntrMetRxnID:
    selectArray = np.logical_or(selectArray, ParDF['Reaction:SBML:reaction:id'].str.contains(rxnIDR).values)

ParDF = ParDF.loc[selectArray]

centralMets = ["g6p_c","f6p_c","man6p_c","gam6p_c","acgam6p_c","acmanap_c","pi_c","amp_c","adp_c","atp_c","acmana_c",
              "nh3_c","pep_c","pyr_c","xu5p__D_c","ru5p__D_c","r5p_c","prpp_c","e4p_c","r1p_c","fdp_c",
              "s7p_c","dhap_c","g3p_c","2dr5p_c","2dr1p_c","nad_c","nadh_c","13dpg_c","nadp_c","nadph_c",
              "3pg_c","2pg_c","lac__L_c","co2_c","acdhlpl_PdhC_c","coa_c","accoa_c","actp_c","ac_c","dhlpl_PdhC_c",
              "lpl_PdhC_c","acald_c","o2_c","co2_c","10fthfglu3_c","thfglu3_c"]

concentrationList = []

metIDcheck = []

for met in centralMets:
    metID = "M_" + met
    conc = ComDF.loc[ ComDF["Compound"] == metID, "InitialConcentration" ].values[0]
    concentrationList.append([metID,conc])
    metIDcheck.append(metID)


# Read through the parameter dataframe and pull out forward and reverse kcat values.
# Then create a dataframe of just the kcat values.

Cent_kcats = []

for rxnID in cntrMetRxnID:
    
    kcatFID = "kcatF_" + rxnID

    kcatF = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[3]
    kcatFunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[3]
    
    kF = [kcatFID,kcatF,kcatFunits]

    kcatRID = "kcatR_" + rxnID
    
    kcatR = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[4]
    kcatRunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[4]
    
    kR = [kcatRID,kcatR,kcatRunits]
    
    Cent_kcats.append(kF)
    Cent_kcats.append(kR)
    
CentKcatDF = pd.DataFrame( Cent_kcats, columns = ['Parameter','Value','Units'])


Cent_Kms = []

for rxnID in cntrMetRxnID:
    
    reaction = RxnDF.loc[ RxnDF["Reaction"] == rxnID, "ReactionFormula" ].values[0]

    reaction = reaction.replace(' ','')

    substrates = reaction.split("<=>")[0]

    if '2.0' in substrates:
        substrates = substrates.replace('2.0','')
    if '3.0' in substrates:
        substrates = substrates.replace('3.0','')
    if '4.0' in substrates:
        substrates = substrates.replace('4.0','')
    
    if '+' in substrates:
        
        plusNum = substrates.count('+')
        
        for i in range(plusNum+1):
            
            substrate = substrates.split('+')[i]


            kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID
    
            km_value = CentQntDF.loc[ CentQntDF["Quantity"] == kmID, "Value" ].values[0]
        
            Cent_Kms.append([kmID,km_value,"mM"])
        
    elif '+' not in substrates:
        
        kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID
    
        km_value = CentQntDF.loc[ CentQntDF["Quantity"] == kmID, "Value" ].values[0]

        Cent_Kms.append([kmID,km_value,"mM"])
    
    products = reaction.split("<=>")[1]

    if '2.0' in products:
        products = products.replace('2.0','')
    if '3.0' in products:
        products = products.replace('3.0','')
    if '4.0' in products:
        products = products.replace('4.0','')

    if '+' in products:
        
        plusNum = products.count('+')
        
        for i in range(plusNum+1):
            
            product = products.split('+')[i]
            
            kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID
    
            km_value = CentQntDF.loc[ CentQntDF["Quantity"] == kmID, "Value" ].values[0]
        
            Cent_Kms.append([kmID,km_value,"mM"])
            
    elif '+' not in products:
        
        kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID
    
        km_value = CentQntDF.loc[ CentQntDF["Quantity"] == kmID, "Value" ].values[0]

        Cent_Kms.append([kmID,km_value,"mM"])
            
CentKmDF = pd.DataFrame(Cent_Kms, columns = ["Parameter","Value","Units"])

############################################

######## Amino Acid Metabolism

############################################

aaMetRxn = ["FMETTRS"]#,"GLNTRAT","GLNTRAT2","METTRS","ILETRS",
            #"VALTRS","LEUTRS","CYSTRS","GLUTRS","GLUTRS_Gln","ARGTRS","TYRTRS",
            #"TRPTRS","SERTRS","THRTRS","PROTRS","ASPTRS","ASNTRS","LYSTRS","HISTRS",
            #"PHETRS","ALATRS","GLYTRS"] # MAT, AHCi
             #"ARGt2r","ASPt2r","CYSt2r","GLUt2r","GLYt2r","ISOt2r",
             #"ALAt2r","ASNt2r","LEUt2r","GLNt2r","HISt2r","LYSt2r","PROt2r","PHEt2r","THRt2r",
             #"TRPt2r","TYRt2r","VALt2r","SERt2r","METt2r"]

aaMetRxnID = [ "R_" + rxn for rxn in aaMetRxn ]

RxnDF_aaMet = RxnDF.loc[ RxnDF.Reaction.isin(aaMetRxnID) ]

RxnDF = pd.concat( [RxnDF_CentralMet, RxnDF_aaMet] )

# Transform SBtab tables into Pandas dataframes
RxnDF = SBtabRxnCent.to_data_frame()
ParDF = SBtabParCent.to_data_frame()
QntDF = SBtabQntCent.to_data_frame()
CentQntDF = SBtabQntCent.to_data_frame()
ComDF = SBtabComCent.to_data_frame()

selectArray = np.ndarray(QntDF.shape[0],dtype=np.bool)

selectArray.fill(False)

# Turns out it isn't very easy to select rows based on a list of substrings that may appear in 
#  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
for rxnIDR in aaMetRxnID:
    selectArray = np.logical_or(selectArray, QntDF['Quantity'].str.contains(rxnIDR).values)

QntDF = QntDF.loc[selectArray]

# Read in the parameters in the parameter dataframe for reactions in Central metabolism.
# This gives us the forward and reverse kcats to use later in kinetic rates.

selectArray = np.ndarray(ParDF.shape[0],dtype=np.bool)
selectArray.fill(False)

for rxnIDR in aaMetRxnID:
    selectArray = np.logical_or(selectArray, ParDF['Reaction:SBML:reaction:id'].str.contains(rxnIDR).values)

ParDF = ParDF.loc[selectArray]

# Read through the parameter dataframe and pull out forward and reverse kcat values.
# Then create a dataframe of just the kcat values.

kcats = []

for rxnID in aaMetRxnID:

    kcatFID = "kcatF_" + rxnID

    kcatF = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[3]
    kcatFunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[3]

    kF = [kcatFID,kcatF,kcatFunits]

    kcatRID = "kcatR_" + rxnID

    kcatR = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[4]
    kcatRunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[4]

    kR = [kcatRID,kcatR,kcatRunits]


    kcats.append(kF)
    kcats.append(kR)

AAKcatDF = pd.DataFrame( kcats, columns = ['Parameter','Value','Units'])

Aa_Kms = []

for rxnID in aaMetRxnID:

    reaction = RxnDF.loc[ RxnDF["Reaction"] == rxnID, "ReactionFormula" ].values[0]

    reaction = reaction.replace(' ','')

    substrates = reaction.split("<=>")[0]

    if '2.0' in substrates:
        substrates = substrates.replace('2.0','')
    if '3.0' in substrates:
        substrates = substrates.replace('3.0','')
    if '4.0' in substrates:
        substrates = substrates.replace('4.0','')

    if '+' in substrates:

        plusNum = substrates.count('+')

        for i in range(plusNum+1):

            substrate = substrates.split('+')[i]


            kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID

            km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]

            Aa_Kms.append([kmID,km_value,"mM"])

    elif '+' not in substrates:

        kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID

        km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]

        Aa_Kms.append([kmID,km_value,"mM"])

    products = reaction.split("<=>")[1]

    if '2.0' in products:
        products = products.replace('2.0','')
    if '3.0' in products:
        products = products.replace('3.0','')
    if '4.0' in products:
        products = products.replace('4.0','')

    if '+' in products:

        plusNum = products.count('+')

        for i in range(plusNum+1):

            product = products.split('+')[i]

            kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID

            km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]

            Aa_Kms.append([kmID,km_value,"mM"])

    elif '+' not in products:

        kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID

        km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]

        Aa_Kms.append([kmID,km_value,"mM"])

AaKmDF = pd.DataFrame(Aa_Kms, columns = ["Parameter","Value","Units"])


aaMetIDs = ["M_ala__L_c", "M_arg__L_c",
    "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_glu__L_c", "M_gln__L_c", "M_gly_c",
    "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
    "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_alatrna_c",    "M_argtrna_c",
    "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c",
    "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c",
    "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c",
    "M_trnaala_c", "M_trnaarg_c",
    "M_trnaasn_c", "M_trnaasp_c", "M_trnacys_c", "M_trnaglu_c", "M_trnagln_c", "M_trnagly_c",
    "M_trnahis_c", "M_trnaile_c", "M_trnaleu_c", "M_trnalys_c", "M_trnamet_c", "M_trnaphe_c",
    "M_trnapro_c", "M_trnaser_c", "M_trnathr_c", "M_trnatrp_c", "M_trnatyr_c", "M_trnaval_c",
    "M_glutrnagln_c", "M_amet_c","M_alatrna_c", "M_argtrna_c",
    "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c",
    "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c",
    "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c"]

aaMetDict = {}
for met in aaMetIDs:
    metID = met
    conc = ComDF.loc[ ComDF["Compound"] == metID, "InitialConcentration" ].values[0]
    aaMetDict[met]=conc
    concentrationList.append([metID,conc])
    metIDcheck.append(metID)

print(aaMetDict)


##########################################

####### Import Nucleotide Metabolism #####

##########################################

file_name = "../model_data/Nucleotide_Kinetic_Parameters.tsv"

# open a file and read it
with open(file_name,"r") as infile:
    file_content = infile.read()

# create an SBtab Document Object Sd
Sd = SBtab.SBtabDocument('BalancedModel', file_content, file_name)

validator = validatorSBtab.ValidateDocument(Sd)
# print("warnings:",validator.validate_document())
warnings = validator.validate_document()
for warn in warnings:
    if warn[1]:
        print(warn[1][0])
        
        
# Read SBTab tables that will be used in the kinetic model
SBtabRxn = Sd.get_sbtab_by_name("Reaction")
SBtabPar = Sd.get_sbtab_by_name("Parameter")
SBtabQnt = Sd.get_sbtab_by_name("Quantity")
SBtabCom = Sd.get_sbtab_by_name("Compound")

# Transform SBtab tables into Pandas dataframes
# RxnDF = SBtabRxn.to_data_frame()
ParDF_nuc = SBtabPar.to_data_frame()
QntDF_nuc = SBtabQnt.to_data_frame()
ComDF_nuc = SBtabCom.to_data_frame()


SBtabRxn = Sd.get_sbtab_by_name("Reaction")
SBtabCom = Sd.get_sbtab_by_name("Compound")
SBtabQnt = Sd.get_sbtab_by_name("Quantity")

ComDF = pd.concat( [ComDF, SBtabCom.to_data_frame()], sort = True  ) 
QntDF_Nuclt = SBtabQnt.to_data_frame()

RxnDF_Nuclt = SBtabRxn.to_data_frame()

ComDF.loc[ ComDF["Compound"] == 'ptsg' ]


nuclMetRxn = list(pd.read_csv("../model_data/nucleo_rxns_list.txt", header=None).loc[:,0].values)
nuclMetRxnID = [ "R_" + x for x in nuclMetRxn]

nuclTurnOff = ['NTD9','DCDPMP','DCTPMP','DCTPDP','CTPDP'] # NTD1, NTD5, NTD6, NTD8
for rxn in nuclTurnOff:
    rxnID = 'R_' + rxn
    nuclMetRxnID.remove(rxnID)


indx = 0
for rxnID in nuclMetRxnID:
    if (rxnID == 'R_PYK'):
        del nuclMetRxnID[indx]
    else:
        indx = indx + 1
indx = 0        
for rxnID in nuclMetRxnID:
    if (rxnID == 'R_PGK'):
        del nuclMetRxnID[indx]
    else:
        indx = indx + 1


RxnDF_Nuclt = RxnDF_Nuclt.loc[ RxnDF_Nuclt.Reaction.isin(nuclMetRxnID) ]


RxnDF = pd.concat( [RxnDF_CentralMet, RxnDF_aaMet, RxnDF_Nuclt] )


# Select only the parameters related to reactions in this section of the model.

selectArray = np.ndarray(QntDF_Nuclt.shape[0],dtype=np.bool)
selectArray.fill(False)

# Turns out it isn't very easy to select rows based on a list of substrings that may appear in 
#  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
for rxnIDR in nuclMetRxnID:
    selectArray = np.logical_or(selectArray, QntDF_Nuclt['Quantity'].str.contains(rxnIDR).values)

QntDF_Nuclt = QntDF_Nuclt.loc[selectArray]

QntDF = pd.concat( [QntDF, QntDF_Nuclt], sort = False  )



nuclMets = ["cytd_c","cmp_c","cdp_c","ctp_c","dcyt_c","dcmp_c","dcdp_c","dctp_c",
           "gsn_c","gua_c","gmp_c","gdp_c","gtp_c","dgsn_c","dgmp_c","dgdp_c","dgtp_c",
           "ura_c","uri_c","ump_c","udp_c","utp_c","duri_c","dump_c","dudp_c","dutp_c",
           "adn_c","ade_c","dad_2_c","damp_c","dadp_c","datp_c","thymd_c","dtmp_c","dtdp_c","dttp_c",
           "glu__L_c","gln__L_c","ppi_c","trdrd_c"]


# concentrationList = []

for met in nuclMets:
    metID = "M_" + met
    conc = ComDF_nuc.loc[ ComDF_nuc["Compound"] == metID, "InitialConcentration" ].values[0]
    concentrationList.append([metID,conc])
    metIDcheck.append(metID)
    #print(metID,conc)


concDF = pd.DataFrame( concentrationList, columns = ['Metabolite','InitialConcentration'] )


# Read through the parameter dataframe and pull out forward and reverse kcat values.
# Then create a dataframe of just the kcat values.

Nuc_kcats = []

for rxnID in nuclMetRxnID:
    
    kcatFID = "kcatF_" + rxnID

    kcatF = ParDF_nuc.loc[ ParDF_nuc["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[3]
    kcatFunits = ParDF_nuc.loc[ ParDF_nuc["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[3]
    
    kF = [kcatFID,kcatF,kcatFunits]

    kcatRID = "kcatR_" + rxnID
    
    kcatR = ParDF_nuc.loc[ ParDF_nuc["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[4]
    kcatRunits = ParDF_nuc.loc[ ParDF_nuc["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[4]
    
    kR = [kcatRID,kcatR,kcatRunits]
    
    Nuc_kcats.append(kF)
    Nuc_kcats.append(kR)
    
NucKcatDF = pd.DataFrame( Nuc_kcats, columns = ['Parameter','Value','Units'])

pykRate = 3204 # Rate kcat 1/s

pgkRate = 220 # Rate kcat 1/s
editRxnsDict = {'R_PYK':pykRate*1,'R_PYK2':pykRate*.11,'R_PYK3':pykRate*.21,
                'R_PYK4':pykRate*0.17,'R_PYK5':pykRate*.18,'R_PYK6':pykRate*0.26,
                'R_PYK7':pykRate*0.13,'R_PYK8':pykRate*.05,'R_PYK9':pykRate*.01,
                'R_PGK':pgkRate*1,'R_PGK2':pgkRate*0.17,'R_PGK3':pgkRate*0.64,
                'R_PGK4':pgkRate*0.16} #scaling values from Pollack et al. OMICS 2002 
tmp = list(NucKcatDF['Value'].values)
tmpRxns = ['R_PYK2','R_PYK3','R_PYK4','R_PYK5','R_PYK6','R_PYK7','R_PYK8','R_PYK9',
          'R_PGK2','R_PGK3','R_PGK4']
for rxn in tmpRxns:
    tmp[list(NucKcatDF['Parameter'].values).index('kcatF_'+rxn)] = editRxnsDict[rxn]
NucKcatDF['Value'] = tmp

KcatDF = pd.concat( [CentKcatDF, NucKcatDF], sort = False )

Nuc_Kms = []

for rxnID in nuclMetRxnID:
    
    reaction = RxnDF.loc[ RxnDF["Reaction"] == rxnID, "ReactionFormula" ].values[0]


    reaction = reaction.replace(' ','')

    substrates = reaction.split("<=>")[0]

    if '2.0' in substrates:
        substrates = substrates.replace('2.0','')
    if '3.0' in substrates:
        substrates = substrates.replace('3.0','')
    if '4.0' in substrates:
        substrates = substrates.replace('4.0','')
    
    if '+' in substrates:
        
        plusNum = substrates.count('+')
        
        for i in range(plusNum+1):
            
            substrate = substrates.split('+')[i]


            kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID

    
            km_value = QntDF_nuc.loc[ QntDF_nuc["Quantity"] == kmID, "Value" ].values[0]

        
            Nuc_Kms.append([kmID,km_value,"mM"])
        
    elif '+' not in substrates:
        
        kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID

    
        km_value = QntDF_nuc.loc[ QntDF_nuc["Quantity"] == kmID, "Value" ].values[0]


        Nuc_Kms.append([kmID,km_value,"mM"])
    
    products = reaction.split("<=>")[1]


    if '2.0' in products:
        products = products.replace('2.0','')
    if '3.0' in products:
        products = products.replace('3.0','')
    if '4.0' in products:
        products = products.replace('4.0','')

    if '+' in products:
        
        plusNum = products.count('+')
        
        for i in range(plusNum+1):
            
            product = products.split('+')[i]
            
            kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID
    
            km_value = QntDF_nuc.loc[ QntDF_nuc["Quantity"] == kmID, "Value" ].values[0]
        
            Nuc_Kms.append([kmID,km_value,"mM"])
            
    elif '+' not in products:
        
        kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID

    
        km_value = QntDF_nuc.loc[ QntDF_nuc["Quantity"] == kmID, "Value" ].values[0]


        Nuc_Kms.append([kmID,km_value,"mM"])
            
NucKmDF = pd.DataFrame(Nuc_Kms, columns = ["Parameter","Value","Units"])



KmDF = pd.concat( [CentKmDF,NucKmDF], sort = False )


##########################################

####### Import Lipid Metabolism #####

##########################################
file_name = "../model_data/lipid_NoH2O_balanced_model.tsv"

# open a file and read it
with open(file_name,"r") as infile:
    file_content = infile.read()

# create an SBtab Document Object Sd
Sd = SBtab.SBtabDocument('BalancedModel', file_content, file_name)

validator = validatorSBtab.ValidateDocument(Sd)
warnings = validator.validate_document()
for warn in warnings:
    if warn[1]:
        print(warn[1][0])

# Read SBTab tables that will be used in the kinetic model
SBtabRxn = Sd.get_sbtab_by_name("Reaction")
SBtabPar = Sd.get_sbtab_by_name("Parameter")
SBtabQnt = Sd.get_sbtab_by_name("Quantity")
SBtabCom = Sd.get_sbtab_by_name("Compound")

# Transform SBtab tables into Pandas dataframes
# RxnDF = SBtabRxn.to_data_frame()
ParDF_lip = SBtabPar.to_data_frame()
QntDF_lip = SBtabQnt.to_data_frame()
ComDF_lip = SBtabCom.to_data_frame()

# Initialize dictionary for "quantity" values
# QntDFDict = dict()

SBtabRxn = Sd.get_sbtab_by_name("Reaction")
SBtabCom = Sd.get_sbtab_by_name("Compound")
SBtabQnt = Sd.get_sbtab_by_name("Quantity")

ComDF = pd.concat( [ComDF, SBtabCom.to_data_frame()], sort = True  )
QntDF_Lipid = SBtabQnt.to_data_frame()
RxnDF_Lipid = SBtabRxn.to_data_frame()

lipMetRxn = ['GLYK','ACPS','BPNT','FAKr','ACPPAT','APG3PAT','AGPAT','DASYN','PGSA','PGPP','CLPNS','PGMT','PAPA','GALU','UDPG4E','UDPGALM','DAGGALT','PSSYN','DAGPST']

lipTurnOff = ['DAGPST','PSSYN']
for rxn in lipTurnOff:
    #rxnID = 'R_' + rxn
    lipMetRxn.remove(rxn)

lipMetRxnID = [ "R_" + x for x in lipMetRxn]
indx = 0
for rxnID in lipMetRxnID:
    ### For now use Central Metabolism for PGMT
    if (rxnID == 'R_PGK'):
        del lipMetRxnID[indx]
         
    else:
        indx = indx + 1
indx = 0
for rxnID in lipMetRxnID:
    if (rxnID == 'R_PYK'):
        del lipMetRxnID[indx]

    else:
        indx = indx + 1

RxnDF_Lipid = RxnDF_Lipid.loc[ RxnDF_Lipid.Reaction.isin(lipMetRxnID) ]

RxnDF = pd.concat( [RxnDF_CentralMet, RxnDF_aaMet,RxnDF_Nuclt, RxnDF_Lipid] )

selectArray = np.ndarray(QntDF_Lipid.shape[0],dtype=np.bool)
selectArray.fill(False)

# Turns out it isn't very easy to select rows based on a list of substrings that may appear in 
#  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
for rxnIDR in lipMetRxnID:
    selectArray = np.logical_or(selectArray, QntDF_Lipid['Quantity'].str.contains(rxnIDR).values)

QntDF_Lipid = QntDF_Lipid.loc[selectArray]

QntDF = pd.concat( [QntDF, QntDF_Lipid], sort = False  )

lipMets = ["udpgal_c","fa_c","ap_c","pa_c","glyc_c","glyc3p_c","apoACP_c","coa_c","pap_c",
          "ACP_c","ACP_R_c","1ag3p_c","cdpdag_c","pg3p_c","pg_c","clpn_c","12dgr_c","galfur12dgr_c",
          "udpgalfur_c","glyc_c","chsterol_c"]

for met in lipMets:
    metID = "M_" + met
    #print(metID)
    conc = ComDF_lip.loc[ ComDF_lip["Compound"] == metID, "InitialConcentration" ].values[0]
    concentrationList.append([metID,conc])
    metIDcheck.append(metID)

concDF = pd.DataFrame( concentrationList, columns = ['Metabolite','InitialConcentration'] )

# Read through the parameter dataframe and pull out forward and reverse kcat values.
# Then create a dataframe of just the kcat values.

Lip_kcats = []

for rxnID in lipMetRxnID:
    
    kcatFID = "kcatF_" + rxnID

    kcatF = ParDF_lip.loc[ ParDF_lip["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[3]
    kcatFunits = ParDF_lip.loc[ ParDF_lip["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[3]
    
    kF = [kcatFID,kcatF,kcatFunits]

    kcatRID = "kcatR_" + rxnID
    
    kcatR = ParDF_lip.loc[ ParDF_lip["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[4]
    kcatRunits = ParDF_lip.loc[ ParDF_lip["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[4]
    
    kR = [kcatRID,kcatR,kcatRunits]
    
    Lip_kcats.append(kF)
    Lip_kcats.append(kR)
    
LipKcatDF = pd.DataFrame( Lip_kcats, columns = ['Parameter','Value','Units'])

Lip_Kms = []

for rxnID in lipMetRxnID:
    
    reaction = RxnDF.loc[ RxnDF["Reaction"] == rxnID, "ReactionFormula" ].values[0]


    reaction = reaction.replace(' ','')

    substrates = reaction.split("<=>")[0]

    if '2.0' in substrates:
        substrates = substrates.replace('2.0','')
    if '3.0' in substrates:
        substrates = substrates.replace('3.0','')
    if '4.0' in substrates:
        substrates = substrates.replace('4.0','')
    if '88.0' in substrates:
        substrates = substrates.replace('88.0','')
    
    if '+' in substrates:
        
        plusNum = substrates.count('+')
        
        for i in range(plusNum+1):
            
            substrate = substrates.split('+')[i]



            kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID

    
            km_value = QntDF_lip.loc[ QntDF_lip["Quantity"] == kmID, "Value" ].values[0]

        
            Lip_Kms.append([kmID,km_value,"mM"])
        
    elif '+' not in substrates:
        
        kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID
        
    
        km_value = QntDF_lip.loc[ QntDF_lip["Quantity"] == kmID, "Value" ].values[0]


        Lip_Kms.append([kmID,km_value,"mM"])
    
    products = reaction.split("<=>")[1]

    if '2.0' in products:
        products = products.replace('2.0','')
    if '3.0' in products:
        products = products.replace('3.0','')
    if '4.0' in products:
        products = products.replace('4.0','')
    if '87.0' in products:
        products = products.replace('87.0','')

    if '+' in products:
        
        plusNum = products.count('+')
        
        for i in range(plusNum+1):
            
            product = products.split('+')[i]
            
            kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID
    
            print(kmID)
    
            km_value = QntDF_lip.loc[ QntDF_lip["Quantity"] == kmID, "Value" ].values[0]
        
            Lip_Kms.append([kmID,km_value,"mM"])
            
    elif '+' not in products:
        
        kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID

    
        km_value = QntDF_lip.loc[ QntDF_lip["Quantity"] == kmID, "Value" ].values[0]

        Lip_Kms.append([kmID,km_value,"mM"])
            
LipKmDF = pd.DataFrame(Lip_Kms, columns = ["Parameter","Value","Units"])


RxnDF = SBtabRxnCent.to_data_frame()
ParDF = SBtabParCent.to_data_frame()
QntDF = SBtabQntCent.to_data_frame()
CentQntDF = SBtabQntCent.to_data_frame()
ComDF = SBtabComCent.to_data_frame()

cofactMetRxn = ["NCTPPRT","NNATr","NADS","NADK","RBFK","FMNAT","5FTHFPGS","FTHFCL","MTHFC","GHMT","MTHFD"]#,"NADHK","GHMT2"]

cofactMetRxnID = [ "R_" + rxn for rxn in cofactMetRxn ]

RxnDF_cofactMet = RxnDF.loc[ RxnDF.Reaction.isin(cofactMetRxnID) ]

RxnDF = pd.concat( [RxnDF_CentralMet, RxnDF_aaMet,RxnDF_Nuclt, RxnDF_Lipid, RxnDF_cofactMet] )

selectArray = np.ndarray(QntDF.shape[0],dtype=np.bool)

selectArray.fill(False)

# Turns out it isn't very easy to select rows based on a list of substrings that may appear in 
#  values of a column. We iterate over reaction IDs and create a logical index of rows to keep.
for rxnIDR in cofactMetRxnID:
    selectArray = np.logical_or(selectArray, QntDF['Quantity'].str.contains(rxnIDR).values)

QntDF = QntDF.loc[selectArray]

# Read in the parameters in the parameter dataframe for reactions in Central metabolism.
# This gives us the forward and reverse kcats to use later in kinetic rates.

selectArray = np.ndarray(ParDF.shape[0],dtype=np.bool)
selectArray.fill(False)

for rxnIDR in cofactMetRxnID:
    selectArray = np.logical_or(selectArray, ParDF['Reaction:SBML:reaction:id'].str.contains(rxnIDR).values)

ParDF = ParDF.loc[selectArray]

# Read through the parameter dataframe and pull out forward and reverse kcat values.
# Then create a dataframe of just the kcat values.

kcats = []

for rxnID in cofactMetRxnID:

    kcatFID = "kcatF_" + rxnID

    kcatF = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[3]
    kcatFunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[3]

    kF = [kcatFID,kcatF,kcatFunits]

    kcatRID = "kcatR_" + rxnID

    kcatR = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Mode" ].values[4]
    kcatRunits = ParDF.loc[ ParDF["Reaction:SBML:reaction:id"] == rxnID, "Unit" ].values[4]

    kR = [kcatRID,kcatR,kcatRunits]


    kcats.append(kF)
    kcats.append(kR)

cofactKcatDF = pd.DataFrame( kcats, columns = ['Parameter','Value','Units'])

KcatDF = pd.concat( [CentKcatDF, NucKcatDF, AAKcatDF, LipKcatDF, cofactKcatDF], sort = False )

cofact_Kms = []

for rxnID in cofactMetRxnID:

    reaction = RxnDF.loc[ RxnDF["Reaction"] == rxnID, "ReactionFormula" ].values[0]


    reaction = reaction.replace(' ','')

    substrates = reaction.split("<=>")[0]


    if '2.0' in substrates:
        substrates = substrates.replace('2.0','')
    if '3.0' in substrates:
        substrates = substrates.replace('3.0','')
    if '4.0' in substrates:
        substrates = substrates.replace('4.0','')

    if '+' in substrates:

        plusNum = substrates.count('+')

        for i in range(plusNum+1):

            substrate = substrates.split('+')[i]


            kmID = 'kmc_' + rxnID + '_' + substrate + '_' + rxnID


            km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]


            cofact_Kms.append([kmID,km_value,"mM"])

    elif '+' not in substrates:

        kmID = 'kmc_' + rxnID + '_' + substrates + '_' + rxnID


        km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]


        cofact_Kms.append([kmID,km_value,"mM"])

    products = reaction.split("<=>")[1]


    if '2.0' in products:
        products = products.replace('2.0','')
    if '3.0' in products:
        products = products.replace('3.0','')
    if '4.0' in products:
        products = products.replace('4.0','')

    if '+' in products:

        plusNum = products.count('+')

        for i in range(plusNum+1):

            product = products.split('+')[i]

            kmID = 'kmc_' + rxnID + '_' + product + '_' + rxnID

            km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]

            cofact_Kms.append([kmID,km_value,"mM"])

    elif '+' not in products:

        kmID = 'kmc_' + rxnID + '_' + products + '_' + rxnID


        km_value = QntDF.loc[ QntDF["Quantity"] == kmID, "Value" ].values[0]


        cofact_Kms.append([kmID,km_value,"mM"])

cofactKmDF = pd.DataFrame(cofact_Kms, columns = ["Parameter","Value","Units"])


cofactMetIDs = ["M_nac_c","M_nicrnt_c","M_dnad_c","M_nh3_c","M_ribflv_c","M_fmn_c","M_fad_c","M_sprm_c","M_thmpp_c",
               "M_5fthf_c","M_5fthfglu3_c","M_methfglu3_c","M_10fthfglu3_c","M_thfglu3_c","M_mlthfglu3_c"]

cofactMetDict = {}
for met in cofactMetIDs:
    metID = met
    conc = ComDF.loc[ ComDF["Compound"] == metID, "InitialConcentration" ].values[0]
    cofactMetDict[met]=conc
    concentrationList.append([metID,conc])
    metIDcheck.append(metID)

concDF = pd.DataFrame( concentrationList, columns = ['Metabolite','InitialConcentration'] )

KmDF = pd.concat( [CentKmDF,NucKmDF,AaKmDF,LipKmDF,cofactKmDF], sort = False )

#####################################################

######### Add Transport Reactions ###################

#####################################################

# file_name = "../model_data/glucose_transport_format.tsv"
# file_name = "../model_data/transport_NoH2O.tsv"
file_name = "../model_data/transport_NoH2O_Zane-TB-DB.tsv" #.tsv"


# open a file and read it
with open(file_name,"r") as infile:
    file_content = infile.read()

# create an SBtab Document Object Sd
Sd = SBtab.SBtabDocument('BalancedModel', file_content, file_name)

validator = validatorSBtab.ValidateDocument(Sd)
warnings = validator.validate_document()
for warn in warnings:
    if warn[1]:
        continue
        

for tab in Sd.sbtabs:
    continue
    # print("Table name:",tab.table_name, "; Table type:",tab.table_type)
    
    
SBtabRxn = Sd.get_sbtab_by_name("Reaction")
SBtabCom = Sd.get_sbtab_by_name("Compound")
SBtabQnt = Sd.get_sbtab_by_name("Quantity")

ComDF = pd.concat( [ComDF, SBtabCom.to_data_frame()], sort = True )
QntDF = pd.concat( [QntDF, SBtabQnt.to_data_frame()], sort = True )

RxnDF_Transport = SBtabRxn.to_data_frame()


## WARNING ## 
# We are missing an ATPase reaction, so there will be no explicit tracking of proton concentration.
# Therefore, we are removing the proton transport reaction.

RxnDF_Transport = RxnDF_Transport.loc[ RxnDF_Transport.Reaction != "R_Ht" ]


# Concatonate the Central and transport reaction dataframes.

RxnDF = pd.concat( [RxnDF_CentralMet, RxnDF_Nuclt, RxnDF_Lipid, RxnDF_aaMet, RxnDF_Transport, RxnDF_cofactMet] )
print(RxnDF)


# Read quantities in from the quantities dataframe.
for items in QntDF.itertuples(index=False):
    if items[0]:
        if len(items[0]) == 1:
            continue
            #print(items)
        if items[0] not in QntDFDict.keys():
            QntDFDict[ items[0] ] = float(items.Value)
    
AvailQnts = list(QntDFDict.keys())
AvailQnts


# sbmlFile = "../model_data/iMB155.xml"
sbmlFile = "../model_data/iMB155_NoH2O.xml"

docSBML = libsbml.readSBMLFromFile(sbmlFile)
modelSBML = docSBML.getModel()

speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
speciesNamesLower = [x.lower() for x in speciesNames]
speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]


# Adapt model to specialized transport rate forms for lipid and glucose transport.

# List of reactions that need to be removed from the model in 
# order to be re-defined. The Reaction Formulas need to be redefined.

transpRxnsList = ["FAt",  "CHOLt", "GLCpts"]

for rxnName in transpRxnsList:
    modelSBML.removeReaction( rxnNamesSBML.index(rxnName)  )
    
    
# List of reactions that need to be created in the SBML model.
transpRxnsList = ["GLCpts0","GLCpts1","GLCpts2","GLCpts3","GLCpts4","GLYCt","FAt","CHOLt"]

newRxnsDF = RxnDF_Transport.loc[ RxnDF_Transport.Name.isin(transpRxnsList) ]
newRxnsDF


rxnCounter = 1

for items in newRxnsDF.itertuples():
    print("[{}/{}]".format(rxnCounter,len(newRxnsDF)), items.Name, "(", items.Reaction, ") ; ",items.ReactionFormula)
    print()
    
    # Create new reaction
    rxnTmp = modelSBML.createReaction()
    rxnTmp.setId("R_" + items.Name)
    rxnTmp.setName(items.Name)
    rxnTmp.setReversible(True)
    
    # Extract reactants and products from the reaction formula
    substratesList = items.ReactionFormula.split("<=>")[0]
    productsList = items.ReactionFormula.split("<=>")[1]
    
    substratesList = substratesList.split("+")
    substratesList = [ met.strip() for met in substratesList]
    
    productsList = productsList.split("+")
    productsList = [ met.strip() for met in productsList]
    
    print(substratesList)
    print(productsList)
    
    # Add substrates to the reaction object
    for spcTmpID in substratesList:
        
        # If the species does not already exist in the SBML model, create it.
        if spcTmpID not in speciesIDs:
            spcTmpObj = modelSBML.createSpecies()
            spcTmpObj.setId(spcTmpID)
            spcTmpObj.setName(spcTmpID)
            spcTmpObj.setCompartment('c')
            spcTmpObj.setConstant(False)
            spcTmpObj.setInitialAmount(0)
            spcTmpObj.setSubstanceUnits('mmole')
            
        species_tmp = rxnTmp.createReactant()
        species_tmp.setSpecies(spcTmpID)
        species_tmp.setStoichiometry(1)
    
    # Add products to the reaction object
    for spcTmpID in productsList:
        
        # If the species does not already exist in the SBML model, create it.
        if spcTmpID not in speciesIDs:
            spcTmpObj = modelSBML.createSpecies()
            spcTmpObj.setId(spcTmpID)
            spcTmpObj.setName(spcTmpID)
            spcTmpObj.setCompartment('c')
            spcTmpObj.setConstant(False)
            spcTmpObj.setInitialAmount(0)
            spcTmpObj.setSubstanceUnits('mmole')
    
        species_tmp = rxnTmp.createProduct()
        species_tmp.setSpecies(spcTmpID)
        species_tmp.setStoichiometry(1)
    
    print("------------------")
    rxnCounter += 1
    
    
# Reset variables with new species and reactions
speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
speciesNamesLower = [x.lower() for x in speciesNames]
speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]


transportRxns = ['GLCpts0','GLCpts1','GLCpts2','GLCpts3','GLCpts4','ACt',
                 'PYRt2r','L_LACt2r','GLYCt','FAt','CHOLt'] #'COAabc']
#                  'NACabc','SPRMabc'] 
#'CYTDabc', 'DCYTabc', 'DURIabc', 'THMDabc', 'ADNabc', 'DADNabc','GSNabc','DGSNabc','URIabc',
# URAt2 #COAabc #NACabc leave this out for now, def in patch
#'ARGt2r','ASPt2pr','GLYt2r','ISOt2r','ALAt2r','ASNt2r','ASNt2r','LEUt2r','HISt2r','LYSt2r','PROt2r','PHEt2r','THRt2r','TRPt2r','TYRt2r','VALt2r','SERt2r','METt2r','CYSt2r','GLUt2pr'



def getSpecIDs(rxnName):
    
    returnList = []
    
    rxnObj = modelSBML.getReaction( rxnNamesSBML.index(rxnName) )
    
    # Use model SBML to get IDs, names, and stoichiometries for reactants
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfReactants() ]
    spcStoich = [ -1*float(x.getStoichiometry()) for x in rxnObj.getListOfReactants() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
    returnList.append( (spcNames, specIDs, spcStoich) )
    
    # Now do the same for products
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfProducts() ]
    spcStoich = [ float(x.getStoichiometry()) for x in rxnObj.getListOfProducts() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
    returnList.append( (spcNames, specIDs, spcStoich) )
    
    return returnList


# Add metabolites to ODE model.

constMetList = []

def addSpcToODE(spcID, concDF, model, pmap, constantVals, explicitParsDict):
    
    # If the species was not defined in reaction input files, return error code
    if spcID not in list(concDF.Metabolite):
        print("ERROR: species ID \"{0}\" not found!".format(spcID) )
        return 1
    
    # If species has already been added to model, return
    if spcID in list(model.getMetDict().keys()):
        return 0
    
    spcName = concDF.loc[ concDF.Metabolite == spcID, "Metabolite" ].values[0]
    
    
    if spcID not in constantVals:
        
        # If it is an external concentration, not previously set, add it as a constant
        if spcID.endswith("_e"):
            constantVals.append(spcID)
            spcConc = concDF.loc[ concDF.Metabolite == spcID, "InitialConcentration" ].values[0]
            model.addParameter("Explicit", spcID, spcConc, unit="mM", 
                   parName="External concentration ({})".format(spcID) )
            
        else:
            spcConc = Rxns.partTomM(pmap[spcID],pmap)
            model.addMetabolite(spcID, spcName, spcConc)
            
    else:
        # For escher maps
        constMetList.append(spcID)
        spcConc = concDF.loc[ concDF.Metabolite == spcID, "InitialConcentration" ].values[0]
        explicitParsDict[spcID] = spcConc
    
    return 0


def addExtSpcToODE(spcID, ComDF, concDF, model, pmap, constantVals, explicitParsDict):
    
    # If the species was not defined in reaction input files, return error code
    if spcID not in list(ComDF.Compound):
        print("ERROR: species ID \"{0}\" not found!".format(spcID) )
        return 1
    
    # If species has already been added to model, return
    if spcID in list(model.getMetDict().keys()):
        return 0
    
    spcName = ComDF.loc[ ComDF.Compound == spcID, "Name" ].values[0]
    
    if spcID not in constantVals:
        
        # If it is an external concentration, not previously set, add it as a constant
        if spcID.endswith("_e"):
            constantVals.append(spcID)
            spcConc = ComDF.loc[ ComDF.Compound == spcID, "InitialConcentration" ].values[0]
            model.addParameter("Explicit", spcID, spcConc, unit="mM", 
                   parName="External concentration ({})".format(spcID) )
            
        else:

            spcConc = Rxns.partTomM(pmap[spcID],pmap)
            model.addMetabolite(spcID, spcName, spcConc)

            
    else:
        # For escher maps
        constMetList.append(spcID)
        spcConc = Rxns.partTomM(pmap[spcID],pmap)
        explicitParsDict[spcID] = spcConc
    
    return 0


# List of reaction IDs in *metabolic* reconstruction
reconstRxnIDList = list(reconstPD["Reaction ID"])

# for indx,row in reconstPD.iterrows():
def createGeneExpression(rxnID, pmap):
    
    
    # Creates protein species and auxiliary functions for effective protein levels.
    # Creates necessary trnscription and translation reactions, when info is available.
    
    if rxnID not in reconstRxnIDList:
        print("WARNING!! Reaction {} not found in reconstruction!".format(rxnID))
        return defaultPtnConcentration
    
    gprStr = str(reconstPD.loc[ reconstPD["Reaction ID"] == rxnID ]["GPR rule"].values[0])
    
    
    # Checks if GPR tule is available.
    if (gprStr == "nan") :
#         print("\n\tReaction", rxnID, "has no GPR!!\n")
        # Keeps default value of 0.001 mM enzyme concentration.
        return rxnID, defaultPtnConcentration
    
    # Get list of genes in MM* code.
    genes = re.findall(r'\bMM[A-Z0-9_]{1,}\b' , gprStr )
    
    # Loop over genes associated with the reaction
    unkIter = 1
    rxnPtns = []
    for mmcode in genes:
        
        # Checks if a translation to JCVISYN2* code is available
        try:
            # See if this is a gap fill
            if mmcode in manGPRPD.MM.values:
                jcvi2ID = "JCVIman_" + mmcode

            else:
                jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
        except:
            print('gpr sadness', mmcode, rxnID)
        
        locusNum = mmcode.split('SYN1_')[1]
        jcvi3AID = 'JCVISYN3A_' + locusNum
        
        
        newMetID = "M_PTN_" + jcvi3AID + "_c"
        
        # Keep list of protein species associated with the current Rxn.
        rxnPtns.append(newMetID)
                
    
    enzLvlVar = "---"
    # Check if this reaction has multiple genes.
    if len(rxnPtns) == 1:
        # If there is just one protein, this is used directly in the reaction's rate law
        enzLvlVar = rxnPtns[0]
        enzConc = Rxns.partTomM(pmap[rxnPtns[0]],pmap)
        
    else:
        # IF there are multiple proteins and a logical rule, we need an extra variable
        # which adds the proteins' concentrations (OR rule) or get the minimal concentration (AND rule).
        
        # To create an "effective" enzyme level, the model will calculate, at every step, the "sum" or "min" 
        #   of multiple protein concentrations.
        enzLvlVar = rxnID+"_effecEnz"
        
        
        if " or " in gprStr.lower():
            enzCnts = 0
            for i in range(len(rxnPtns)):
                cnt = pmap[rxnPtns[i]]
                enzCnts = enzCnts + cnt
            enzConc = Rxns.partTomM(enzCnts,pmap)
        else:
            enzCnts = []
            for i in range(len(rxnPtns)):
                cnt = pmap[rxnPtns[i]]
                enzCnts.append(cnt)
            enzCnt = min(enzCnts)
            enzConc = Rxns.partTomM(enzCnt,pmap)
    return enzLvlVar, enzConc


def addGlobalParams(model):
# Parse global parameters to model
    model.parseParameters("../model_data/GlobalParameters_Zane-TB-DB.csv")

    explPar, globPar, rxnPar = model.getParameters()

    # Store explicit and global parameters (fixed values used throughout the model)
    # They may contain metabolites with constant concentration, from the external media
    explicitParsDict = {}
    for name,par in explPar.items():
        explicitParsDict[par.formKey] = par.val

    constantVals = list(explicitParsDict.keys())
    
    return constantVals, explicitParsDict


def addReactionsToModel(model, pmap):

    # To store the rate forms and their names
    rateFormsDict = {}
    knownForms = set()

    ### Diagnostics:
    # parameters with *zero/inf/nan* value:
    errorPar = []

    rxnCounter = 1
    
    constantVals, explicitParsDict = addGlobalParams(model)
    
    for items in RxnDF.itertuples():


        if (items.Name in nuclMetRxn) or (items.Name in cntrMetRxn) or (items.Name in lipMetRxn) or (items.Name in aaMetRxn) or (items.Name in cofactMetRxn):


            # For clarity
            rxnID = items.Name


            EnzymeParam, EnzymeConc = createGeneExpression(rxnID, pmap)

            specLists = getSpecIDs(rxnID)

            reactantsDict = {}

            # Matches the substrates in the rate form to standardized IDs.
            substratesList = []
            substrates = 0
            spcCounter = 1
                
            for i in range(len(specLists[0][1])):

                spc = specLists[0][1][i]
                spcStoich = specLists[0][2][i]

                for j in range(int(abs(spcStoich))):
                    subsStr = "Sub" + str(spcCounter)
                    substratesList.append( (spc, subsStr, spcStoich) )
                    spcCounter +=1


                # we create this lookup dict to standardize parameter names later on
                reactantsDict[spc] = subsStr

                for j in range(int(abs(spcStoich))):
                    substrates = substrates + 1

            # Matches the products in the rate form to standardized IDs.
            productsList = []

            spcCounter = 1
            products = 0

                
            for i in range(len(specLists[1][1])):

                spc = specLists[1][1][i]
                spcStoich = specLists[1][2][i]

                for j in range(int(abs(spcStoich))):
                    prodStr = "Prod" + str(spcCounter)
                    productsList.append( (spc, prodStr, spcStoich) )
                    spcCounter +=1



                reactantsDict[spc] = prodStr


                for j in range(int(abs(spcStoich))):
                    products = products + 1



            rxnName = 'R_' + rxnID

            rateLaw = Rxns.Enzymatic(substrates,products)


            RateName = rxnID + '_Rate'

            model.addRateForm(RateName, odecell.modelbuilder.RateForm(rateLaw))

            rxnIndx = model.addReaction(rxnID, RateName, rxnName="Reaction " + rxnID)

            for spc in reactantsDict.keys():
                if spc in list(model.getMetDict().keys()):
                    # In case we already know about this metabolite.
                    continue
                if spc in metIDcheck:
                    if addSpcToODE(spc, concDF, model, pmap, constantVals, explicitParsDict):
                        # In case the metabolite ID referenced by a rate form is not found
                        #  among defined metabolites.
                        metError = True
                        break

                else:
                    if addExtSpcToODE(spc, ComDF, concDF, model, pmap, constantVals, explicitParsDict):
                        # In case the metabolite ID referenced by a rate form is not found
                        #  among defined metabolites.
                        metError = True
                        break


            rxnMetsAdded = []
                    
            for (modelID, stdID, spcStoich) in substratesList:
                if modelID in constantVals:
                    model.addParameter(rxnIndx, stdID, modelID )
                elif modelID in rxnMetsAdded:
                    model.addParameter(rxnIndx, stdID, modelID )
                else:
                    model.addSubstrate(rxnIndx, stdID, modelID, stoich=spcStoich)
                    rxnMetsAdded.append(modelID)

                km_ID = 'kmc_R_' + rxnID + '_' + modelID + '_R_' + rxnID

                km_Value = KmDF.loc[ KmDF["Parameter"] == km_ID, "Value" ].values[0]

                KM_rxnID = "Km" + stdID.replace('$','')

                model.addParameter(rxnID, KM_rxnID, km_Value)


                    
            for (modelID, stdID, spcStoich) in productsList:
                if modelID in constantVals:
                    model.addParameter(rxnIndx, stdID, modelID )
                elif modelID in rxnMetsAdded:
                    model.addParameter(rxnIndx, stdID, modelID )
                else:
                    model.addProduct(rxnIndx, stdID, modelID, stoich=spcStoich)
                    rxnMetsAdded.append(modelID)

                km_ID = 'kmc_R_' + rxnID + '_' + modelID + '_R_' + rxnID

                km_Value = KmDF.loc[ KmDF["Parameter"] == km_ID, "Value" ].values[0]

                KM_rxnID = "Km" + stdID.replace('$','')

                model.addParameter(rxnID, KM_rxnID, km_Value)


            kcatFID = "kcatF_" + rxnName
            kcatF = KcatDF.loc[ KcatDF["Parameter"] == kcatFID, "Value" ].values[0]

            model.addParameter(rxnIndx, 'kcatF', kcatF)

            kcatRID = "kcatR_" + rxnName
            kcatR = KcatDF.loc[ KcatDF["Parameter"] == kcatRID, "Value" ].values[0]

            model.addParameter(rxnIndx, 'kcatR', kcatR)

                # Add enzyme level parameter.
            if "Enzyme" in model.getReaction(rxnID).getKeys():

                model.addParameter(rxnIndx, "Enzyme", EnzymeConc, unit="mM", parName="Enzyme concentration" )

            if "onoff" in model.getReaction(rxnID).getKeys():
                model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#             print("\n---------------------------\n")
            rxnCounter += 1




        elif items.Name in transportRxns:
            

            kl = re.sub('\s{2,}', ' ', items.KineticLaw)


            # For clarity
            rxnID = items.Name

            specLists = getSpecIDs(rxnID)

            reactantsDict = {}

            # Matches the substrates in the rate form to standardized IDs.
            substratesList = []
            spcCounter = 1
            for i in range(len(specLists[0][1])):

                spc = specLists[0][1][i]
                spcStoich = specLists[0][2][i]

                subsStr = "Sub" + str(spcCounter)
#                 print( "Setting {:<6} to ({}) {}".format(subsStr, spcStoich, spc) )

                substratesList.append( (spc, "$"+subsStr, spcStoich) )

                # we create this lookup dict to standardize parameter names later on
                reactantsDict[spc] = subsStr

                # We use regular expression syntax to assure that no substrings will
                # be matched. For example, when substituting the species M_abc_c for $Sub1, 
                # we would "sub-match" the parameter kmc_R_XYZ_M_abc_c, and create kmc_R_XYZ_Sub1,
                # which would be bad...
                kl = re.sub( r'\b' + str(spc) + r'\b', "$"+subsStr, kl)

                spcCounter +=1

            # Matches the products in the rate form to standardized IDs.
            productsList = []
            spcCounter = 1
            for i in range(len(specLists[1][1])):

                spc = specLists[1][1][i]
                spcStoich = specLists[1][2][i]

                prodStr = "Prod" + str(spcCounter)
#                 print( "Setting {:<6} to ({}) {}".format(prodStr, spcStoich, spc) )

                productsList.append( (spc, "$"+prodStr, spcStoich) )
                reactantsDict[spc] = prodStr

                kl = re.sub( r'\b' + str(spc) + r'\b', "$"+prodStr, kl)

                spcCounter +=1


            # Break main loop
            metError = False

            # Sanity Checks that the metabolites found here are in the ODECell model.
            for spc in reactantsDict.keys():
                if spc in list(model.getMetDict().keys()):
                    # In case we already know about this metabolite.
                    continue
                if spc in metIDcheck:
                    if addSpcToODE(spc, concDF, model, pmap, constantVals, explicitParsDict):
                        # In case the metabolite ID referenced by a rate form is not found
                        #  among defined metabolites.
                        metError = True
                        break

                else:
                    if addExtSpcToODE(spc, ComDF, concDF, model, pmap, constantVals, explicitParsDict):
                        # In case the metabolite ID referenced by a rate form is not found
                        #  among defined metabolites.
                        metError = True
                        break
                        
            if metError:
                break

            paramsDict = {}
            paramCount = defaultdict(int)

            # Substitute parameter values
            for qnt in AvailQnts:
                # Checks that the parameter refers to the reaction
                # Each parameter has the "R_XXX" reaction code in its name.
                # Needed because "R_AB1" parameters would be found in "R_AB12" laws.
                # Check if the quantity is an explicit paramter (a global constant)
                #   shared between rate laws.
                if (items.Reaction not in qnt) or (qnt in constantVals):
                    continue

                # Sanity Checks if the parameter is found in the kinetic law
                if qnt not in kl:
        #             print("WARNING: quantity {} not found reaction {}.".format(qnt, items.Reaction) )
                    continue

                paramType = ""

                # Classify the parameter, store ID, name and value.
                if "hco" in qnt:
                    paramType = "hco"
                    paramCount["hco"] += 1
                    paramsDict[qnt] = ("$hco"+str(paramCount["hco"]), "reaction cooperativity" )
#                     print(  '{:<32} - {:10} - {:25} - {}'.format(
#                         "Found cooperativity constant", paramsDict[qnt][0], qnt, QntDFDict[qnt])  )

                elif "keq" in qnt:
                    paramType = "keq"
                    paramCount["keq"] += 1
                    paramsDict[qnt] = ("$keq"+str(paramCount["keq"]), "equilibrium constant" )
#                     print(  '{:<32} - {:10} - {:25} - {}'.format(
#                         "Found equilibrium constant", paramsDict[qnt][0], qnt, QntDFDict[qnt])  )

                elif "kcr" in qnt:
                    paramType = "kcr"
                    paramCount["kcr"] += 1
                    paramsDict[qnt] = ("$kcat"+str(paramCount["kcr"]), "catalitic rate constant" )
#                     print(  '{:<32} - {:10} - {:25} - {}'.format(
#                         "Found catalitic rate constant", paramsDict[qnt][0], qnt, QntDFDict[qnt])  )

                elif "km" in qnt:
                    paramType = "km"
                    paramCount["km"] += 1

                    # Find species related to the Km value
                    relatedSpcs = set()
                    for spc in [ x[0] for x in substratesList + productsList ]:
                        if spc in qnt:
                            relatedSpcs.add(spc)

                    if not len(relatedSpcs):
                        print("WARNING: no species found for parameter", qnt)
                    else:
                        paramsDict[qnt] = ("$km_" + "_".join( [ reactantsDict[x] for x in relatedSpcs] ) , 
                                           "species constant" )

                else:
                    paramType = "OTHER"
#                     print("\n\nWARNING: parameter {} could not be classified!".format(qnt))
                    paramCount["OTHER"] += 1

                    paramsDict[qnt] = ("$Const"+str(paramCount["OTHER"]), "rate law-specific constant" )


                # Sanity check:
                if QntDFDict[qnt] == 0 or np.isnan(QntDFDict[qnt]) or np.isinf(QntDFDict[qnt]) :
                    errorPar.append(( qnt, paramType, QntDFDict[qnt]))
#                     print("\n\tWARNING!! Parameter {} has value {}!".format(qnt, QntDFDict[qnt]) )
                    if 'keq' in qnt:
                        QntDFDict[qnt] = 10**-7


                # Use regular expression to extract *only* the parameter ID and substitute it
                # with a standardized ID.
                kl = re.sub( r'\b' + str(qnt) + r'\b', paramsDict[qnt][0], kl)

            # Transform "power" operator to Python notation.
            kl = kl.replace("^", "**")

            # Transform "square root" operator to NumPY notation.
            kl = kl.replace("sqrt", "np.sqrt")

            enzLevlPar = ""
            # Create argument for protein level.
            if re.search( "^\( 0.001 / 1.0\)", kl):
                # If the search returns anything but a "None", it found a match.
                # We need the test to create a reaction-specific parameter later on.
                # Creates key for rate law.
                enzLevlPar = "$enzLvl"

                # Substitute hard-coded default value for key in rate law.
                kl = re.sub( "^\( 0.001 / 1.0\)", enzLevlPar, kl)



            # DEBUG -- create an ON/OFF switch for all metabolic reactions:
            kl = "$onoff * (" + kl + " )"



            if kl not in knownForms:
                # If we found a new rate form, report on it and store it for easy search.
                rateFormsDict[kl] = ("RateForm" + str(len(rateFormsDict)), 
                                     "{}s_{}p".format(len(substratesList), len(productsList)) )
                knownForms.add(kl)

                # Now add the rate form to the model.
                model.addRateForm(rateFormsDict[kl][0], odecell.modelbuilder.RateForm(kl))

            if enzLevlPar:
                # In case we are explicitly tracking enzyme levels:
                # Add gene expression for necessary enzyme(s) for the current reaction.
                # This has to run *before*  the actual reaction is added in case a GPR rule 
                #   needs an auxiliary reaction to determine effective enzyme levels.
                enzLvlVar = createGeneExpression(rxnID)

            # Add reaction to model, using the rate form defined above
            rxnIndx = model.addReaction(rxnID, rateFormsDict[kl][0], rxnName="Reaction " + rxnID)

            # Add substrate species
            for (modelID, stdID, spcStoich) in substratesList:
                if modelID in constantVals:
                    model.addParameter(rxnIndx, stdID[1:], modelID )
                else:
                    model.addSubstrate(rxnIndx, stdID[1:], modelID, stoich=spcStoich)

            # Add product species
            for (modelID, stdID, spcStoich) in productsList:
                if modelID in constantVals:
                    model.addParameter(rxnIndx, stdID[1:], modelID )
                else:
                    model.addProduct(rxnIndx, stdID[1:], modelID, stoich=spcStoich)

            # Add reaction-specific parameter values
            for modelID,val in paramsDict.items():
                stdID = val[0][1:]
                model.addParameter(rxnIndx, stdID, QntDFDict[modelID], unit="", parName=val[1] )

            # Add enzyme level parameter.
            if "enzLvl" in model.getReaction(rxnID).getKeys():
                model.addParameter(rxnIndx, "enzLvl", enzLvlVar, unit="mM", parName="Enzyme concentration" )

            # DEBUG
            if "onoff" in model.getReaction(rxnID).getKeys():
                model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

            rxnCounter += 1

        #### PUNP5  --- Explicit Addition

#RateForm
    model.explicitTwoSubTwoProd = odecell.modelbuilder.RateForm('$onoff * $EnzConc * (($Ksub*($Sub1/$KmSub1)*($Sub2/$KmSub2)-$Kprod*($Prod1/$KmProd1)*($Prod2/$KmProd2)) / ((1+$Sub1/$KmSub1)*(1+$Sub2/$KmSub2)+(1+$Prod1/$KmProd1)*(1+$Prod2/$KmProd2)-1))')
    model.updateAvailableForms()

    model.addMetabolite('M_uri_c','uridine',Rxns.partTomM(pmap['M_uri_c'],pmap)) # mM - from pb file  
    

    rxnIndx = model.addReaction('PUNP5','explicitTwoSubTwoProd','Nucleotide Metabolism Reactions R_PUNP5')
  
    enzConc_PUNP5 = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0747_c'],pmap) # mM
    model.addParameter('PUNP5','EnzConc',enzConc_PUNP5)#'M_PTN_MMSYN1_0747_c')
    model.addParameter('PUNP5','Ksub',26.123325)#avg value from the other version (PUNP1-4)
    #BRENDA has a value for Kcat to be 0.09 1/s in Plasmodium falciparum
    model.addParameter('PUNP5','Kprod',109.8933)#avg value from the other version

    model.addSubstrate('PUNP5','Sub1','M_uri_c')#,stoich=-1.0)
    model.addParameter('PUNP5','KmSub1',0.115)#From BRENDA Plasmodium falciparum at pH=7.4 have (0.085 mM)
    model.addSubstrate('PUNP5','Sub2','M_pi_c')#,stoich=-1.0)
    model.addParameter('PUNP5','KmSub2',2.1566)#avg value from the other version (PUNP1-4)

    model.addProduct('PUNP5','Prod1','M_ura_c')#,stoich=1.0)
    model.addParameter('PUNP5','KmProd1',2.8972)#avg value from the other version (PUNP1-4)
    model.addProduct('PUNP5','Prod2','M_r1p_c')#,stoich=1.0)
    model.addParameter('PUNP5','KmProd2',0.0137)#avg value from the other version (PUNP1-4)

    model.addParameter('PUNP5', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")

    rPUNP5=model.getReaction('PUNP5')
    
    ### ADD ATPase
    enzConc_ATPase = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0789_c'],pmap)
    #RateForm
    ATPaseRateLaw = Rxns.Enzymatic(2,1)

    RateName = 'ATPase_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(ATPaseRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('ATPase','ATPase_Rate','ATP synthase')
    model.addParameter('ATPase','Enzyme',enzConc_ATPase)
    model.addParameter('ATPase','kcatF',20)

    model.addParameter('ATPase','kcatR',217/3) #285

    model.addSubstrate('ATPase','Sub1','M_adp_c')
    model.addParameter('ATPase','KmSub1',0.1)
    model.addSubstrate('ATPase','Sub2','M_pi_c')
    model.addParameter('ATPase','KmSub2',4.2)

    model.addProduct('ATPase','Prod1','M_atp_c')
    model.addParameter('ATPase','KmProd1',0.6)

    model.addParameter('ATPase', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch") 
    
    
    enzConc_PIabc = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0427_c'],pmap)
    
    PitRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'Piabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(PitRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('PIabc','Piabc_Rate','Pi transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_PIabc)

    model.addParameter('PIabc','kcatF',25) 
    model.addParameter('PIabc','kcatR',0)

    pi_e_conc = 134.0 # mM, Growth Medium Buffers

    model.addParameter('PIabc','Sub1',pi_e_conc)
    model.addParameter('PIabc','KmSub1',0.0031)
    model.addSubstrate('PIabc','Sub2','M_atp_c')
    model.addParameter('PIabc','KmSub2',0.023)

    model.addProduct('PIabc','Prod1','M_pi_c',stoich=2)
    model.addParameter('PIabc','KmProd1',0.02)
    model.addParameter('PIabc','Prod2','M_pi_c')
    model.addParameter('PIabc','KmProd2',0.385)
    model.addProduct('PIabc','Prod3','M_adp_c')
    model.addParameter('PIabc','KmProd3',0.654)

    model.addParameter('PIabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    enzConc_NUCabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0010_c'],pmap['M_PTN_JCVISYN3A_0009_c']),pmap)
    
    AdnRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'ADNabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(AdnRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('ADNabc','ADNabc_Rate','ADN transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('ADNabc','kcatF', 2.0) #2.0
    model.addParameter('ADNabc','kcatR',0)

    adn_e_conc = 0.15

    model.addParameter('ADNabc','Sub1',adn_e_conc)
    model.addParameter('ADNabc','KmSub1',1.9/10**3)
    model.addSubstrate('ADNabc','Sub2','M_atp_c')
    model.addParameter('ADNabc','KmSub2',0.6)

    model.addProduct('ADNabc','Prod1','M_adn_c')
    model.addParameter('ADNabc','KmProd1',0.02)
    model.addProduct('ADNabc','Prod2','M_pi_c')
    model.addParameter('ADNabc','KmProd2',10)
    model.addProduct('ADNabc','Prod3','M_adp_c')
    model.addParameter('ADNabc','KmProd3',2.8)
    
    model.addParameter('ADNabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    DAdnRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'DADNabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(DAdnRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('DADNabc','DADNabc_Rate','DADN transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('DADNabc','kcatF',0.7) # 0.5
    model.addParameter('DADNabc','kcatR',0)

    dad_2_e_conc = 0.02

    model.addParameter('DADNabc','Sub1',adn_e_conc)
    model.addParameter('DADNabc','KmSub1',1.9/10**3)
    model.addSubstrate('DADNabc','Sub2','M_atp_c')
    model.addParameter('DADNabc','KmSub2',0.6)

    model.addProduct('DADNabc','Prod1','M_dad_2_c')
    model.addParameter('DADNabc','KmProd1',0.02)
    model.addProduct('DADNabc','Prod2','M_pi_c')
    model.addParameter('DADNabc','KmProd2',10)
    model.addProduct('DADNabc','Prod3','M_adp_c')
    model.addParameter('DADNabc','KmProd3',2.8)
    
    model.addParameter('DADNabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    GsnRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'GSNabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(GsnRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('GSNabc','GSNabc_Rate','GSN transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('GSNabc','kcatF', 1.0) #0.75)
    model.addParameter('GSNabc','kcatR',0)

    gsn_e_conc = 0.13

    model.addParameter('GSNabc','Sub1',gsn_e_conc)
    model.addParameter('GSNabc','KmSub1',1.9/10**3)
    model.addSubstrate('GSNabc','Sub2','M_atp_c')
    model.addParameter('GSNabc','KmSub2',0.6)

    model.addProduct('GSNabc','Prod1','M_gsn_c')
    model.addParameter('GSNabc','KmProd1',0.02)
    model.addProduct('GSNabc','Prod2','M_pi_c')
    model.addParameter('GSNabc','KmProd2',10)
    model.addProduct('GSNabc','Prod3','M_adp_c')
    model.addParameter('GSNabc','KmProd3',2.8)
    
    model.addParameter('GSNabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    DGsnRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'DGSNabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(DGsnRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('DGSNabc','DGSNabc_Rate','DGSN transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('DGSNabc','kcatF', 0.7) # 0.5
    model.addParameter('DGSNabc','kcatR',0)

    dgsn_e_conc = 0.02

    model.addParameter('DGSNabc','Sub1',dgsn_e_conc)
    model.addParameter('DGSNabc','KmSub1',1.9/10**3)
    model.addSubstrate('DGSNabc','Sub2','M_atp_c')
    model.addParameter('DGSNabc','KmSub2',0.6)

    model.addProduct('DGSNabc','Prod1','M_dgsn_c')
    model.addParameter('DGSNabc','KmProd1',0.02)
    model.addProduct('DGSNabc','Prod2','M_pi_c')
    model.addParameter('DGSNabc','KmProd2',10)
    model.addProduct('DGSNabc','Prod3','M_adp_c')
    model.addParameter('DGSNabc','KmProd3',2.8)
    
    model.addParameter('DGSNabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    UriRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'URIabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(UriRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('URIabc','URIabc_Rate','URI transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('URIabc','kcatF', 2.0) #2
    model.addParameter('URIabc','kcatR',0)

    uri_e_conc = 0.18

    model.addParameter('URIabc','Sub1',uri_e_conc)
    model.addParameter('URIabc','KmSub1',1.7/10**3)
    model.addSubstrate('URIabc','Sub2','M_atp_c')
    model.addParameter('URIabc','KmSub2',0.6)

    model.addProduct('URIabc','Prod1','M_uri_c')
    model.addParameter('URIabc','KmProd1',0.02)
    model.addProduct('URIabc','Prod2','M_pi_c')
    model.addParameter('URIabc','KmProd2',10)
    model.addProduct('URIabc','Prod3','M_adp_c')
    model.addParameter('URIabc','KmProd3',2.8)
    
    model.addParameter('URIabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    DCytRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'DCYTabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(DCytRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('DCYTabc','DCYTabc_Rate','DCYT transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('DCYTabc','kcatF', 0.5) #0.5
    model.addParameter('DCYTabc','kcatR',0)

    dcyt_e_conc = 0.02

    model.addParameter('DCYTabc','Sub1',dcyt_e_conc)
    model.addParameter('DCYTabc','KmSub1',1.1/10**3)
    model.addSubstrate('DCYTabc','Sub2','M_atp_c')
    model.addParameter('DCYTabc','KmSub2',0.6)

    model.addProduct('DCYTabc','Prod1','M_dcyt_c')
    model.addParameter('DCYTabc','KmProd1',0.02)
    model.addProduct('DCYTabc','Prod2','M_pi_c')
    model.addParameter('DCYTabc','KmProd2',10)
    model.addProduct('DCYTabc','Prod3','M_adp_c')
    model.addParameter('DCYTabc','KmProd3',2.8)
    
    model.addParameter('DCYTabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    ThmdRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'THMDabc_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(ThmdRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('THMDabc','THMDabc_Rate','THMYD transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NUCabc)

    model.addParameter('THMDabc','kcatF',0.75) # 0.75
    model.addParameter('THMDabc','kcatR',0)

    thmyd_e_conc = 0.08 # 0.08 + 0.02 C5Mod-CMRL

    model.addParameter('THMDabc','Sub1',thmyd_e_conc)
    model.addParameter('THMDabc','KmSub1',1.7/10**3)
    model.addSubstrate('THMDabc','Sub2','M_atp_c')
    model.addParameter('THMDabc','KmSub2',0.6)

    model.addProduct('THMDabc','Prod1','M_thymd_c')
    model.addParameter('THMDabc','KmProd1',0.02)
    model.addProduct('THMDabc','Prod2','M_pi_c')
    model.addParameter('THMDabc','KmProd2',10)
    model.addProduct('THMDabc','Prod3','M_adp_c')
    model.addParameter('THMDabc','KmProd3',2.8)
    
    model.addParameter('THMDabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    SprmabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'SPRMabc_Rate'
    
    enzConc_SPRMabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0195_c'],pmap['M_PTN_JCVISYN3A_0196_c'],pmap['M_PTN_JCVISYN3A_0197_c']),pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(SprmabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('SPRMabc','SPRMabc_Rate','SPRM transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_SPRMabc)

    model.addParameter('SPRMabc','kcatF',3)
    model.addParameter('SPRMabc','kcatR',0)

    sprm_e_conc = 0.1

    model.addParameter('SPRMabc','Sub1',sprm_e_conc)
    model.addParameter('SPRMabc','KmSub1',2)
    model.addSubstrate('SPRMabc','Sub2','M_atp_c')
    
    model.addMetabolite('M_sprm_c','sprm',Rxns.partTomM(pmap['M_sprm_c'],pmap))
    model.addProduct('SPRMabc','Prod1','M_sprm_c')
    model.addParameter('SPRMabc','KmProd1',2)
    model.addProduct('SPRMabc','Prod2','M_pi_c')
    model.addProduct('SPRMabc','Prod3','M_adp_c')
    
    model.addParameter('SPRMabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    NacabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'NACabc_Rate'
    
    enzConc_NACabc = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0314_c'],pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(NacabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('NACabc','NACabc_Rate','NAC transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_NACabc)

    model.addParameter('NACabc','kcatF',1.66)
    model.addParameter('NACabc','kcatR',0)

    nac_e_conc = 0.016

    model.addParameter('NACabc','Sub1',nac_e_conc)
    model.addParameter('NACabc','KmSub1',0.0019)
    model.addSubstrate('NACabc','Sub2','M_atp_c')

    model.addProduct('NACabc','Prod1','M_nac_c')
    model.addParameter('NACabc','KmProd1',0.0019)
    model.addProduct('NACabc','Prod2','M_pi_c')
    model.addProduct('NACabc','Prod3','M_adp_c')
    
    model.addParameter('NACabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    RIBFLVabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'RIBFLVabc_Rate'
    
    enzConc_RIBFLVabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0641_c'],pmap['M_PTN_JCVISYN3A_0642_c'],pmap['M_PTN_JCVISYN3A_0643_c'],pmap['M_PTN_JCVISYN3A_0877_c']),pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(RIBFLVabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('RIBFLVabc','RIBFLVabc_Rate','RIBFLV transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_RIBFLVabc)

    model.addParameter('RIBFLVabc','kcatF',1.66)
    model.addParameter('RIBFLVabc','kcatR',0)

    ribflv_e_conc = 0.003

    model.addParameter('RIBFLVabc','Sub1',ribflv_e_conc)
    model.addParameter('RIBFLVabc','KmSub1',0.0019)
    model.addSubstrate('RIBFLVabc','Sub2','M_atp_c')

    model.addProduct('RIBFLVabc','Prod1','M_ribflv_c')
    model.addParameter('RIBFLVabc','KmProd1',0.0019)
    model.addProduct('RIBFLVabc','Prod2','M_pi_c')
    model.addProduct('RIBFLVabc','Prod3','M_adp_c')
    
    model.addParameter('RIBFLVabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    FTHFabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'FTHFabc_Rate'
    
    enzConc_FTHFabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0641_c'],pmap['M_PTN_JCVISYN3A_0642_c'],pmap['M_PTN_JCVISYN3A_0643_c'],pmap['M_PTN_JCVISYN3A_0822_c']),pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(FTHFabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('5FTHFabc','FTHFabc_Rate','FTHF transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_FTHFabc)

    model.addParameter('5FTHFabc','kcatF',1.66)
    model.addParameter('5FTHFabc','kcatR',0)

    fthf_e_conc = 0.05

    model.addParameter('5FTHFabc','Sub1',fthf_e_conc)
    model.addParameter('5FTHFabc','KmSub1',0.0019)
    model.addSubstrate('5FTHFabc','Sub2','M_atp_c')

    model.addProduct('5FTHFabc','Prod1','M_5fthf_c')
    model.addParameter('5FTHFabc','KmProd1',0.0019)
    model.addProduct('5FTHFabc','Prod2','M_pi_c')
    model.addProduct('5FTHFabc','Prod3','M_adp_c')
    
    model.addParameter('5FTHFabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    THMPPabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'THMPPabc_Rate'
    
    enzConc_THMPPabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0706_c'],pmap['M_PTN_JCVISYN3A_0707_c'],pmap['M_PTN_JCVISYN3A_0708_c']),pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(THMPPabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('THMPPabc','THMPPabc_Rate','THMPP transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_THMPPabc)

    model.addParameter('THMPPabc','kcatF',1.66)
    model.addParameter('THMPPabc','kcatR',0)

    thmpp_e_conc = 0.008

    model.addParameter('THMPPabc','Sub1',thmpp_e_conc)
    model.addParameter('THMPPabc','KmSub1',0.0019)
    model.addSubstrate('THMPPabc','Sub2','M_atp_c')

    model.addMetabolite('M_thmpp_c','thmpp',Rxns.partTomM(pmap['M_thmpp_c'],pmap))
    model.addProduct('THMPPabc','Prod1','M_thmpp_c')
    model.addParameter('THMPPabc','KmProd1',0.0019)
    model.addProduct('THMPPabc','Prod2','M_pi_c')
    model.addProduct('THMPPabc','Prod3','M_adp_c')
    
    model.addParameter('THMPPabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    P5PabcRateLaw = Rxns.Enzymatic(1,1)

    RateName = 'P5Pabc_Rate'
    
    enzConc_P5Pabc = Rxns.partTomM(min(pmap['M_PTN_JCVISYN3A_0641_c'],pmap['M_PTN_JCVISYN3A_0642_c'],pmap['M_PTN_JCVISYN3A_0643_c'],pmap['M_PTN_JCVISYN3A_0345_c']),pmap)

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(P5PabcRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('P5Pabc','P5Pabc_Rate','P5P transport')
    model.addParameter(rxnIndx,'Enzyme',enzConc_P5Pabc)

    model.addParameter('P5Pabc','kcatF',1.66)
    model.addParameter('P5Pabc','kcatR',0)

    p5p_e_conc = 0.006

    model.addParameter('P5Pabc','Sub1',p5p_e_conc)
    model.addParameter('P5Pabc','KmSub1',0.0019)
    model.addSubstrate('P5Pabc','Sub2','M_atp_c')

    model.addMetabolite('M_pydx5p_c','pydx5p',Rxns.partTomM(pmap['M_pydx5p_c'],pmap))
    model.addProduct('P5Pabc','Prod1','M_pydx5p_c')
    model.addParameter('P5Pabc','KmProd1',0.0019)
    model.addProduct('P5Pabc','Prod2','M_pi_c')
    model.addProduct('P5Pabc','Prod3','M_adp_c')
    
    model.addParameter('P5Pabc', 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    # Ion Transporters

    Kt6RateLaw = Rxns.Enzymatic(3,4)

    RateName = 'Kt6_Rate'

    model.addMetabolite("M_k_c","potassium",Rxns.partTomM(pmap['M_k_c'],pmap))

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(Kt6RateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('Kt6','Kt6_Rate','K transport')
    
    enzConc_Kt6 = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0686_c'],pmap)
    
    model.addParameter(rxnIndx,'Enzyme',enzConc_Kt6)

    model.addParameter(rxnIndx,'kcatF',3)
    model.addParameter(rxnIndx,'kcatR',0)

    k_e_conc = 10 # 12.66 mM C5Mod-CMRL
    na_c_conc = 10
    na_e_conc = 134.0 # mM - Kim 2/5 update CMRL-C5Mod #6.67 # mM C5mod-CMRL #140

    model.addParameter(rxnIndx,'Sub1',k_e_conc)
    model.addParameter(rxnIndx,'KmSub1',0.46)
    model.addSubstrate(rxnIndx,'Sub2','M_atp_c')
    model.addParameter(rxnIndx,'KmSub2',0.03)
    model.addParameter(rxnIndx,'Sub3',na_c_conc)
    model.addParameter(rxnIndx,'KmSub3',7.47)

    model.addProduct(rxnIndx,'Prod1','M_k_c')
    model.addParameter(rxnIndx,'KmProd1',1.9)
    model.addProduct(rxnIndx,'Prod2','M_adp_c')
    model.addParameter(rxnIndx,'KmProd2',1.5)
    model.addProduct(rxnIndx,'Prod3','M_pi_c')
    model.addParameter(rxnIndx,'KmProd3',10)
    model.addParameter(rxnIndx,'Prod4',na_e_conc)
    model.addParameter(rxnIndx,'KmProd4',12.7)

    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    CA2tRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'CA2_Rate'

    model.addMetabolite("M_ca2_c","calcium",Rxns.partTomM(pmap['M_ca2_c'],pmap))

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(CA2tRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('CA2abc','CA2_Rate','CA2 transport')
    enzConc_cat2 = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0879_c'],pmap)
    
    model.addParameter(rxnIndx,'Enzyme',enzConc_cat2)

    model.addParameter(rxnIndx,'kcatF',9.5)
    model.addParameter(rxnIndx,'kcatR',0)

    ca2_e_conc = 0.68 # mM CMRL-C5mod #0.05

    model.addParameter(rxnIndx,'Sub1',ca2_e_conc)
    model.addParameter(rxnIndx,'KmSub1',0.0075)
    model.addSubstrate(rxnIndx,'Sub2','M_atp_c')
    model.addParameter(rxnIndx,'KmSub2',0.075)

    model.addProduct(rxnIndx,'Prod1','M_ca2_c')
    model.addParameter(rxnIndx,'KmProd1',0.0075)
    model.addProduct(rxnIndx,'Prod2','M_adp_c')
    model.addParameter(rxnIndx,'KmProd2',0.1)
    model.addProduct(rxnIndx,'Prod3','M_pi_c')
    model.addParameter(rxnIndx,'KmProd3',10)

    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    MG2tRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'MG2_Rate'

    model.addMetabolite("M_mg2_c","magnesium",Rxns.partTomM(pmap['M_mg2_c'],pmap))

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(MG2tRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('MG2abc','MG2_Rate','MG2 transport')
    enzConc_mg2 = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0787_c'],pmap)
    
    model.addParameter(rxnIndx,'Enzyme',enzConc_mg2)

    model.addParameter(rxnIndx,'kcatF',22)
    model.addParameter(rxnIndx,'kcatR',0)

    mg2_e_conc = 1 # 1.4 mM - CMRL-C5Mod

    model.addParameter(rxnIndx,'Sub1',mg2_e_conc)
    model.addParameter(rxnIndx,'KmSub1',0.05)
    model.addSubstrate(rxnIndx,'Sub2','M_atp_c')
    model.addParameter(rxnIndx,'KmSub2',2.8)

    model.addProduct(rxnIndx,'Prod1','M_mg2_c')
    model.addParameter(rxnIndx,'KmProd1',0.05)
    model.addProduct(rxnIndx,'Prod2','M_adp_c')
    model.addParameter(rxnIndx,'KmProd2',2.8)
    model.addProduct(rxnIndx,'Prod3','M_pi_c')
    model.addParameter(rxnIndx,'KmProd3',10)

    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
    
    
    enzConc_GluTrans = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0886_c'],pmap)
    enzConc_AaTransSymp = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0876_c'],pmap) + Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0878_c'],pmap)
    enzConc_AaTransAbc = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0168_c'],pmap) #0165
    
    kcat_aa_abc = 0.2 #1/s

    aaTransport = [['R_ARGt2r', 'R_ARG4abc', 'M_arg__L_c', 0.03, 3.4, kcat_aa_abc, 4],
        ['R_ASPt2pr', 'R_ASP4abc', 'M_asp__L_c', 0.03, 9.1, kcat_aa_abc, 0.1],
        ['R_GLYt2r', 'R_GLY4abc', 'M_gly_c', 0.03, 5.1, kcat_aa_abc, 4],
        ['R_ISOt2r', 'R_ISO4abc', 'M_ile__L_c', 0.03, 9.1, kcat_aa_abc, 4],
        ['R_ALAt2r', 'R_ALA4abc', 'M_ala__L_c', 0.03, 5.1, kcat_aa_abc, 4],
        ['R_ASNt2r', 'R_ASN4abc', 'M_asn__L_c', 0.03, 6.1, kcat_aa_abc, 4],
        ['R_LEUt2r', 'R_LEU4abc', 'M_leu__L_c', 0.03, 9.1, kcat_aa_abc, 4.2],
        ['R_HISt2r', 'R_HIS4abc', 'M_his__L_c', 0.03, 3.4, kcat_aa_abc, 4],
        ['R_LYSt2r', 'R_LYS4abc', 'M_lys__L_c', 0.03, 9.1, kcat_aa_abc, 4.2],
        ['R_PROt2r', 'R_PRO4abc', 'M_pro__L_c', 0.03, 3.4, kcat_aa_abc, 4],
        ['R_PHEt2r', 'R_PHE4abc', 'M_phe__L_c', 0.03, 3.4, kcat_aa_abc, 1],
        ['R_THRt2r', 'R_THR4abc', 'M_thr__L_c', 0.03, 5.1, kcat_aa_abc, 1],
        ['R_TRPt2r', 'R_TRP4abc', 'M_trp__L_c', 0.03, 3.4, kcat_aa_abc, 0.5],
        ['R_TYRt2r', 'R_TYR4abc', 'M_tyr__L_c', 0.03, 3.4, kcat_aa_abc, 4],
        ['R_VALt2r', 'R_VAL4abc', 'M_val__L_c', 0.03, 5.1, kcat_aa_abc, 8],
        ['R_SERt2r', 'R_SER4abc', 'M_ser__L_c', 0.03, 6.1, kcat_aa_abc, 1.1],
        ['R_METt2r', 'R_MET4abc', 'M_met__L_c', 0.03, 3.4, kcat_aa_abc, 8],
        ['R_CYSt2r', 'R_CYS4abc', 'M_cys__L_c', 0.03, 3.4, kcat_aa_abc, 2.3],
        ['R_GLUt2pr', 'R_GLU4abc', 'M_glu__L_c', 0.03, 3.4, kcat_aa_abc, 4.3],
        ['R_GLNt2r', 'R_GLN4abc', 'M_gln__L_c', 0.03, 7.1, kcat_aa_abc, 2]] 
    
    for rxn in aaTransport:
        rxnIDsymp = rxn[0] #+ "_zero"
        rxnIDabc = rxn[1]
        metID = rxn[2]
        Km = rxn[3]
        KcatSymp = rxn[4]
        KcatAbc = rxn[5]
        ExtConc = rxn[6]

        if ('glu' not in metID) and ('gln' not in metID) and ('gly' not in metID) and ('ser' not in metID):
#             print(metID)
            model.addMetabolite(metID, metID, Rxns.partTomM(pmap[metID],pmap))
            
        AaTransSympRateLaw = Rxns.Enzymatic(1,1)

        RateName = rxnIDsymp + '_Rate'

        model.addRateForm(RateName, odecell.modelbuilder.RateForm(AaTransSympRateLaw))
        model.updateAvailableForms()
 
        rxnIndx = model.addReaction(rxnIDsymp,RateName,metID + ' Transport')
        
        if 'glu' in metID:
            model.addParameter(rxnIDsymp,'Enzyme',enzConc_GluTrans)
        else:
            model.addParameter(rxnIDsymp,'Enzyme',enzConc_AaTransSymp)

        model.addParameter(rxnIDsymp,'kcatF', KcatSymp)
        model.addParameter(rxnIDsymp,'kcatR',0)
        
        model.addParameter(rxnIDsymp,'Sub1',ExtConc)
        model.addParameter(rxnIDsymp,'KmSub1',Km)
        model.addProduct(rxnIDsymp,'Prod1',metID)
        model.addParameter(rxnIDsymp,'KmProd1',Km)

        model.addParameter(rxnIDsymp, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")
        
        
        AaTransAbcRateLaw = Rxns.Enzymatic(1,1)

        RateName = rxnIDabc + '_Rate'

        model.addRateForm(RateName, odecell.modelbuilder.RateForm(AaTransAbcRateLaw))
        model.updateAvailableForms()
 
        rxnIndx = model.addReaction(rxnIDabc,RateName,metID + ' Transport')
        
        model.addParameter(rxnIDabc,'Enzyme',enzConc_AaTransAbc)

        model.addParameter(rxnIDabc,'kcatF', KcatAbc)
        model.addParameter(rxnIDabc,'kcatR',0)
        
        model.addParameter(rxnIDabc,'Sub1',ExtConc)
        model.addParameter(rxnIDabc,'KmSub1',Km)
        model.addProduct(rxnIDabc,'Prod1',metID)
        model.addParameter(rxnIDabc,'KmProd1',Km)
        
        model.addSubstrate(rxnIDabc,'Sub2','M_atp_c')
        model.addParameter(rxnIDabc,'KmSub2',2.8)

        model.addProduct(rxnIDabc,'Prod2','M_adp_c')
        model.addParameter(rxnIDabc,'KmProd2',2.8)
        model.addProduct(rxnIDabc,'Prod3','M_pi_c')
        model.addParameter(rxnIDabc,'KmProd3',10)

        model.addParameter(rxnIDabc, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")

    ### Missing Metabolic Reactions
    ## GHMT2: h2o_c + methfglu3_c --> 5fthfglu3_c + h_c
    ###No exisiting data other than Keq
    model.addRateForm("GHMT2_Rate", odecell.modelbuilder.RateForm(Rxns.Enzymatic(1,1)))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('GHMT2','GHMT2_Rate','enzymatic reaction GHMT2')
    model.addParameter(rxnIndx,'Enzyme',Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0799_c'],pmap))
    model.addParameter(rxnIndx,'kcatF',640)#Taken from GHMT reaction
    model.addParameter(rxnIndx,'kcatR',23.173)#Taken from GHMT reaction

    model.addSubstrate(rxnIndx,'Sub1','M_methfglu3_c')
    model.addParameter(rxnIndx,'KmSub1',0.0684)#Taken from GHMT reaction, related substrate
    model.addProduct(rxnIndx,'Prod1','M_5fthfglu3_c')
    model.addParameter(rxnIndx,'KmProd1',0.1463)#Taken from GHMT reaction, related substrate
    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")

    ## NADHK: atp_c + nadh_c --> adp_c + h_c + nadph_c
    ###exists in the TSV file but not read, code below should not be needed
    model.addRateForm("NADHK_Rate", odecell.modelbuilder.RateForm(Rxns.Enzymatic(2,2)))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('NADHK','NADHK_Rate','enzymatic reaction NADHK')
    model.addParameter(rxnIndx,'Enzyme',Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0259_c'],pmap))
    model.addParameter(rxnIndx,'kcatF',31.3405)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
    model.addParameter(rxnIndx,'kcatR',104.2685)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv

    model.addSubstrate(rxnIndx,'Sub1','M_nadh_c')
    model.addParameter(rxnIndx,'KmSub1',2.0)#0.1251)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
    model.addSubstrate(rxnIndx,'Sub2','M_atp_c')
    model.addParameter(rxnIndx,'KmSub2',0.4054)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv

    model.addProduct(rxnIndx,'Prod1','M_nadph_c')
    model.addParameter(rxnIndx,'KmProd1',0.3)#7.4674)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
    model.addProduct(rxnIndx,'Prod2','M_adp_c')
    model.addParameter(rxnIndx,'KmProd2',4.624)#From Parameter File, Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv
    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")


    spcConc = Rxns.partTomM(pmap['M_glc__D_c'],pmap)
    model.addMetabolite('M_glc__D_c', 'M_glc__D_c', spcConc)

    GLCKRateLaw = Rxns.Enzymatic(2,3)

    RateName = 'GLCK_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(GLCKRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('GLCK','GLCK_Rate','glucokinase')
    enzConc_glck = Rxns.partTomM(pmap['M_PTN_JCVISYN3A_0250_c'],pmap)
    
    model.addParameter(rxnIndx,'Enzyme',enzConc_glck)

    model.addParameter(rxnIndx,'kcatF',6.3)
    model.addParameter(rxnIndx,'kcatR',0)

    model.addSubstrate(rxnIndx,'Sub1','M_glc__D_c')
    model.addParameter(rxnIndx,'KmSub1',0.78)
    model.addSubstrate(rxnIndx,'Sub2','M_atp_c')
    model.addParameter(rxnIndx,'KmSub2',0.5)

    model.addProduct(rxnIndx,'Prod1','M_g6p_c')
    model.addParameter(rxnIndx,'KmProd1',0.1)
    model.addProduct(rxnIndx,'Prod2','M_adp_c')
    model.addParameter(rxnIndx,'KmProd2',1.0)
    model.addProduct(rxnIndx,'Prod3','M_pi_c')
    model.addParameter(rxnIndx,'KmProd3',10)

    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")


    GLCtRateLaw = '$onoff * $Enzyme * ( ( $kcatF * $Sub1 ) - ( $kcatR * $Prod1 ) )'

    RateName = 'GLCt_Rate'

    model.addRateForm(RateName, odecell.modelbuilder.RateForm(GLCtRateLaw))
    model.updateAvailableForms()

    rxnIndx = model.addReaction('GLCt','GLCt_Rate','glucose uptake')
    enzConc_glct = Rxns.partTomM(pmap['ptsg'],pmap)
    
    model.addParameter(rxnIndx,'Enzyme',enzConc_glct)

    model.addParameter(rxnIndx,'kcatF',4,000)
    model.addParameter(rxnIndx,'kcatR',4,000*40)

    model.addParameter(rxnIndx,'Sub1','M_glc__D_e')

    model.addProduct(rxnIndx,'Prod1','M_glc__D_c')

    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch")

    
    print("Defined reactions")

    return

"""
Author: Zane Thornburg
"""

##### CME Model #####

from GIP_rates import *

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

def constructCME(csim,pmap):
    
    add_CME_species(csim, pmap)
    
    genes_in_model = []
    
    PtnMetDF = pd.read_csv("../model_data/protein_metabolites_frac.csv")
    riboPtnMetDF = pd.read_csv("../model_data/ribo_protein_metabolites.csv")
    memPtnMetDF = pd.read_csv("../model_data/membrane_protein_metabolites.csv")
    trnaMetDF = pd.read_csv("../model_data/trna_metabolites_synthase.csv")
    
    ctRNAcostMap = OrderedDict({"M_alatrna_c":"ALA_cost","M_argtrna_c":"ARG_cost","M_asntrna_c":"ASN_cost",
                            "M_asptrna_c":"ASP_cost","M_cystrna_c":"CYS_cost","M_glutrna_c":"GLU_cost",
                            "M_glntrna_c":"GLN_cost","M_glytrna_c":"GLY_cost","M_histrna_c":"HIS_cost",
                            "M_iletrna_c":"ILE_cost","M_leutrna_c":"LEU_cost","M_lystrna_c":"LYS_cost",
                            "M_mettrna_c":"MET_cost","M_phetrna_c":"PHE_cost","M_protrna_c":"PRO_cost",
                            "M_sertrna_c":"SER_cost","M_thrtrna_c":"THR_cost","M_trptrna_c":"TRP_cost",
                            "M_tyrtrna_c":"TYR_cost","M_valtrna_c":"VAL_cost"})
    
    for index, row in PtnMetDF.iterrows():
        addNamedPtnCME(csim, row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"], pmap, genes_in_model)

    for index, row in memPtnMetDF.iterrows():
        addMembranePtnCME(csim, row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"], pmap, genes_in_model)

    for index, row in riboPtnMetDF.iterrows():
        addRiboPtnCME(csim, row["species"], row["gene"], pmap, genes_in_model)   

    for gene in genomePtnLocDict:
        if gene not in genes_in_model:
            addPtnCME(csim, gene, pmap, genes_in_model)
    
    tRNAadded = []

    for index, row in trnaMetDF.iterrows():
        addtRNACME(csim, row["uncharged"], row["charged"], row["gene"], row["AA"], row["synthase"], tRNAadded, pmap, ctRNAcostMap)


def add_CME_species(csim, pmap):
    
    aaMetIDs = ["M_ala__L_c", "M_arg__L_c","M_glu__L_c","M_gln__L_c",
        "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_gly_c",
        "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
        "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_amet_c"]
    
    glugln = 'M_glutrnagln_c'
    glugln_enz = 'glutrnagln_enz'
    glugln_enz_atp = 'glutrnagln_enz_atp'
    glugln_enz_atp_aa = 'glutrnagln_enz_atp_gln'
    glugln_enz_atp_asn = 'glutrnagln_enz_atp_asn'

    glugln_list = [glugln,glugln_enz,glugln_enz_atp,glugln_enz_atp_aa,glugln_enz_atp_asn]
    
    for key, val in pmap.items():
        
        if 'New_' in key:
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
            
        if (key=='M_atp_c') or (key=='M_adp_c') or (key=='M_amp_c') or (key=='M_ppi_c') or (key=='M_pi_c'):
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
            
        if ('RP_' in key) or ('RP2_' in key):
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
            
        if ('_AA' in key) or ('ATP_' in key) or ('_ATP' in key) or ('_tRNA' in key) or ('P_' in key):
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
            
        if ('M_trna' in key) or ('trna_c' in key):
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
            
        if '_cost' in key:
            
            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
        
        if key in aaMetIDs:

            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
        
        if key in glugln_list:

            csim.defineSpecies([key])
            if val < 0:
                print(key,val)
            csim.addParticles(key,count=int(val))
            
#             print(key,val)
    

def addtRNACME(csim, unchargedMetID, chargedMetID, mmcode, aminoAcid, synthase, tRNAadded, pmap, ctRNAcostMap):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')

    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        # If not, check for a manual AOE connection:
        if mmcode in manGPRPD.MM.values:
            jcvi2ID = "JCVIman_" + mmcode
        else:
            # If not, set an "unknown" code
            jcvi2ID = "JCVIunk_" + mmcode + "_" + str(unkIter)
            unkIter += 1 

#     print(jcvi3AID)

    geneMetID = jcvi3AID + '_gene'

    rnasequence = getRNAsequences(jcvi3AID)
#     print(rnasequence)

    gene = [geneMetID]
#     print(gene)

    new_RNA_ID = 'New_RNA_' + locusNum
    TrscProd = [new_RNA_ID]
#         mRNAdegProd = []

    for i in range(len(rnasequence)):
        TrscProd.append('ATP_trsc')

    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)

    # Add total number of monomers to parameter dict

    N_A = baseCount["A"]

    N_U = baseCount["U"]

    N_C = baseCount["C"]

    N_G = baseCount["G"]

    for i in range(N_A):
        TrscProd.append('ATP_tRNA')

    for i in range(N_U):
        TrscProd.append('UTP_tRNA')


    for i in range(N_G):
        TrscProd.append('GTP_tRNA')


    for i in range(N_C):
        TrscProd.append('CTP_tRNA')


    RNAP_gene = 'RP_' + locusNum

    csim.addReaction(RNAP_gene, tuple(TrscProd), trnaTranscriptRateCME(rnasequence, pmap))

    if unchargedMetID not in tRNAadded:

        tRNA = [unchargedMetID]

        tRNAadded.append(unchargedMetID)

        ctRNA = [chargedMetID]

        tRNAadded.append(chargedMetID)
    
        if 'GLN' not in synthase:
            
            synthase = synthase.split('3A_')[1].lstrip('0')
            synthase_ptn = 'P_' + synthase
    
            synthase_atp = synthase_ptn + '_ATP'

            
            synthase_atp_aa = synthase_atp + '_AA'

            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNA'
            

            csim.addReaction(('M_atp_c',synthase_ptn),(synthase_atp),0.01)
            csim.addReaction((aminoAcid,synthase_atp),(synthase_atp_aa),0.01)
            csim.addReaction((synthase_atp_aa,unchargedMetID),(synthase_atp_aa_trna),0.1)
            csim.addReaction((synthase_atp_aa_trna),(chargedMetID,'M_amp_c','M_ppi_c',synthase_ptn),30)
            
            costID = ctRNAcostMap[chargedMetID]
#             print(costID)
            
            cost_paid = costID + '_paid'
            
            csim.addReaction((chargedMetID,costID),(unchargedMetID,cost_paid),100)
            
#             print('Added charging reactions for '+unchargedMetID)
            
        elif 'GLN' in synthase:
            
            synthase_ptn = 'P_126'
            synthase_atp = synthase_ptn + '_ATP'
            synthase_atp_aa = synthase_atp + '_AA'
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNAgln'
            synth_atp_aa_trna = [synthase_atp_aa_trna]

            
            csim.addReaction((synthase_atp_aa,'M_trnagln_c'),(synthase_atp_aa_trna),0.1)
            csim.addReaction((synthase_atp_aa_trna),('M_glutrnagln_c','M_amp_c','M_ppi_c',synthase_ptn),30)
            
            glugln = 'M_glutrnagln_c'
            glugln_enz = 'glutrnagln_enz'
            glugln_enz_atp = 'glutrnagln_enz_atp'
            glugln_enz_atp_aa = 'glutrnagln_enz_atp_gln'
            glugln_enz_atp_asn = 'glutrnagln_enz_atp_asn'

            
            csim.addReaction((glugln,'P_689'),(glugln_enz),0.1)
            csim.addReaction((glugln_enz,'M_atp_c'),(glugln_enz_atp),0.01)
            csim.addReaction((glugln_enz_atp,'M_gln__L_c'),(glugln_enz_atp_aa),0.01)
            csim.addReaction((glugln_enz_atp_aa),('M_glntrna_c','M_glu__L_c','M_adp_c','M_pi_c','P_689'),30)
            csim.addReaction((glugln_enz_atp,'M_asn__L_c'),(glugln_enz_atp_asn),0.001)
            csim.addReaction((glugln_enz_atp_asn),('M_glntrna_c','M_asp__L_c','M_adp_c','M_pi_c','P_689'),30)

            costID = ctRNAcostMap['M_glntrna_c']
#             print(costID)
            
            cost_paid = costID + '_paid'
            
            csim.addReaction(('M_glntrna_c',costID),('M_trnagln_c',cost_paid),100)
            
            print('GLN added')

    

def addMembranePtnCME(csim, newMetID, mmcode, transcribe, ptnFrac, pmap, genes_in_model):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    genes_in_model.append(jcvi3AID)

    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        # If not, check for a manual AOE connection:
        if mmcode in manGPRPD.MM.values:
            jcvi2ID = "JCVIman_" + mmcode
        else:
            # If not, set an "unknown" code
            jcvi2ID = "JCVIunk_" + mmcode + "_" + str(unkIter)

#     print(mmcode, jcvi2ID, jcvi3AID)

    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)
    
    
    if transcribe:

        geneMetID = jcvi3AID + '_gene'


        rnasequence, aasequence = getSequences(jcvi3AID)
#         print(rnasequence)
#         print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName


        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]

        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')


        baseCount = defaultdict(int)
        for base in set(rnasequence):
            baseCount[base] = rnasequence.count(base)

        # Add total number of monomers to parameter dict

        N_A = baseCount["A"]

        N_U = baseCount["U"]

        N_C = baseCount["C"]

        N_G = baseCount["G"]

        for i in range(N_A):
            TrscProd.append('ATP_mRNA')


        for i in range(N_U):
            TrscProd.append('UTP_mRNA')


        for i in range(N_G):
            TrscProd.append('GTP_mRNA')


        for i in range(N_C):
            TrscProd.append('CTP_mRNA')


        RNAP_gene = 'RP_' + locusNum

        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(newMetID, rnasequence, jcvi2ID, pmap))
    

def addRiboPtnCME(csim, newMetID, mmcode, pmap, genes_in_model):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')

    genes_in_model.append(jcvi3AID)

#     # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        # If not, check for a manual AOE connection:
        if mmcode in manGPRPD.MM.values:
            jcvi2ID = "JCVIman_" + mmcode
        else:
            # If not, set an "unknown" code
            jcvi2ID = "JCVIunk_" + mmcode + "_" + str(unkIter)
            unkIter += 1 

#     print(mmcode, jcvi2ID, jcvi3AID)

    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)
    
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):

        rnaMetID = "M_RNA_" + jcvi3AID + "_c"
        rnaName = "(mRNA) " + ptnName
        
    
    trsc_rate = riboTranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID)

    new_RNA_ID = 'New_mRNA_' + locusNum

    TrscProd = [new_RNA_ID]

    
    for i in range(len(rnasequence)):
        TrscProd.append('ATP_trsc')

        
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)

    # Add total number of monomers to parameter dict

    N_A = baseCount["A"]

    N_U = baseCount["U"]

    N_C = baseCount["C"]

    N_G = baseCount["G"]

    for i in range(N_A):
        TrscProd.append('ATP_mRNA')


    for i in range(N_U):
        TrscProd.append('UTP_mRNA')


    for i in range(N_G):
        TrscProd.append('GTP_mRNA')


    for i in range(N_C):
        TrscProd.append('CTP_mRNA')

    RNAP_gene = 'RP_' + locusNum
    
    csim.addReaction(RNAP_gene, tuple(TrscProd), riboTranscriptRateCME( newMetID, rnasequence, jcvi2ID, pmap))


        

def addNamedPtnCME(csim, newMetID, mmcode, transcribe, ptnFrac, pmap, genes_in_model):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    genes_in_model.append(jcvi3AID)

#     # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        # If not, check for a manual AOE connection:
        if mmcode in manGPRPD.MM.values:
            jcvi2ID = "JCVIman_" + mmcode
        else:
            # If not, set an "unknown" code
            jcvi2ID = "JCVIunk_" + mmcode + "_" + str(unkIter)

#     print(mmcode, jcvi2ID, jcvi3AID)

    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)
    
    
    if transcribe:


        rnasequence, aasequence = getSequences(jcvi3AID)
#         print(rnasequence)
#         print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName


        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]


        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')

        baseCount = defaultdict(int)
        for base in set(rnasequence):
            baseCount[base] = rnasequence.count(base)

        # Add total number of monomers to parameter dict

        N_A = baseCount["A"]

        N_U = baseCount["U"]

        N_C = baseCount["C"]

        N_G = baseCount["G"]

        for i in range(N_A):
            TrscProd.append('ATP_mRNA')


        for i in range(N_U):
            TrscProd.append('UTP_mRNA')


        for i in range(N_G):
            TrscProd.append('GTP_mRNA')


        for i in range(N_C):
            TrscProd.append('CTP_mRNA')


        RNAP_gene = 'RP_' + locusNum

        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(newMetID, rnasequence, jcvi2ID, pmap))

        

def addPtnCME(csim, jcvi3AID, pmap, genes_in_model):
    locusNum = jcvi3AID.split('_')[1]
    mmcode = 'MMSYN1_' + locusNum
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode

#     print(mmcode, jcvi2ID, jcvi3AID)
    
    genes_in_model.append(jcvi3AID)

    ptnMetID = 'P_' + locusNum
    # If the protein is not in the model, add it:

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)

    # Get nucleotide and amino acid sequences, if available
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):

        
        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]
    
        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')

            
        baseCount = defaultdict(int)
        for base in set(rnasequence):
            baseCount[base] = rnasequence.count(base)

        # Add total number of monomers to parameter dict

        N_A = baseCount["A"]

        N_U = baseCount["U"]

        N_C = baseCount["C"]

        N_G = baseCount["G"]
        
        for i in range(N_A):
            TrscProd.append('ATP_mRNA')

            
        for i in range(N_U):
            TrscProd.append('UTP_mRNA')

            
        for i in range(N_G):
            TrscProd.append('GTP_mRNA')

            
        for i in range(N_C):
            TrscProd.append('CTP_mRNA')

        RNAP_gene = 'RP_' + locusNum
    
        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(ptnMetID, rnasequence, jcvi2ID, pmap))

        
        

baseMap = OrderedDict({ "A":"M_atp_c", "U":"M_utp_c", "G":"M_gtp_c", "C":"M_ctp_c" })
# baseMapToMonoP = OrderedDict({ "A":"M_amp_c", "U":"M_ump_c", "G":"M_gmp_c", "C":"M_cmp_c" })

# Global parameters for transcription
rnaPolKcat = 20 # nt/s
rnaPolK0 = 1e-4 #mM
rnaPolKd = 0.1 #mM

rrnaPolKcat = 85 # nt/s

krnadeg = 0.00578/2 # 1/s
# rna_deg_rate = sim.rateConst('RNAdeg', krnadeg, 2)

ptnDegRate = 7.70e-06 # 1/s

ATPconc = 4 #mM
UTPconc = 3 #mM
CTPconc = 1 #mM
GTPconc = 2 #mM

# Cell radius (meters):
# r_cell = 2.5*(10**-7)
r_cell = 2.0*(10**-7) # m

CytoVolume = (4*np.pi/3)*1000*r_cell**3 # L
cellVolume = CytoVolume

subvolume_vol = 1000*(8e-9)**3 # L

# print(cellVolume)

# Avogadro:
avgdr   = 6.022e23 # molec/mol
Avognum = avgdr

NaV = Avognum * subvolume_vol

countToMiliMol = 1000/(avgdr*cellVolume)

RnaPconc = 187*countToMiliMol # mM

# Global parameter for degradation of mRNAs
rnaDegRate = 0.00578/2 # 1/s

# degrad_bind_rate = 11/60/countToMiliMol*1000 #1/M/s
Ecoli_V = 1e-15 #L
degrad_bind_rate = 11*avgdr*Ecoli_V/60/2400 #1/M/s

# Create a map for rna sequence to NTP concentration.
baseMap = OrderedDict({ "A":ATPconc, "U":UTPconc, "G":GTPconc, "C":CTPconc })

# Create Dictionaries to map tRNAs to associated aa abbreviations in protein sequences.
aaMap = OrderedDict({"A":"M_ala__L_c", "R":"M_arg__L_c", 
    "N":"M_asn__L_c", "D":"M_asp__L_c", "C":"M_cys__L_c", "E":"M_glu__L_c", "Q":"M_gln__L_c", "G":"M_gly_c", 
    "H":"M_his__L_c", "I":"M_ile__L_c", "L":"M_leu__L_c", "K":"M_lys__L_c", "M":"M_met__L_c", "F":"M_phe__L_c", 
    "P":"M_pro__L_c", "S":"M_ser__L_c", "T":"M_thr__L_c", "W":"M_trp__L_c", "Y":"M_tyr__L_c", "V":"M_val__L_c",
    "*":"Stop_Codon"})

aaTRNAMap = OrderedDict({"A":"M_alatrna_c", "R":"M_argtrna_c", 
    "N":"M_asntrna_c", "D":"M_asptrna_c", "C":"M_cystrna_c", "E":"M_glutrna_c", "Q":"M_glntrna_c", "G":"M_glytrna_c", 
    "H":"M_histrna_c", "I":"M_iletrna_c", "L":"M_leutrna_c", "K":"M_lystrna_c", "M":"M_mettrna_c", "F":"M_phetrna_c", 
    "P":"M_protrna_c", "S":"M_sertrna_c", "T":"M_thrtrna_c", "W":"M_trptrna_c", "Y":"M_tyrtrna_c", "V":"M_valtrna_c"})

aaTRNAFreeMap = OrderedDict({"A":"M_trnaala_c", "R":"M_trnaarg_c", 
    "N":"M_trnaasn_c", "D":"M_trnaasp_c", "C":"M_trnacys_c", "E":"M_trnaglu_c", "Q":"M_trnagln_c", "G":"M_trnagly_c", 
    "H":"M_trnahis_c", "I":"M_trnaile_c", "L":"M_trnaleu_c", "K":"M_trnalys_c", "M":"M_trnamet_c", "F":"M_trnaphe_c", 
    "P":"M_trnapro_c", "S":"M_trnaser_c", "T":"M_trnathr_c", "W":"M_trnatrp_c", "Y":"M_trnatyr_c", "V":"M_trnaval_c"})


defaultPtnCount = 9

# Global parameters for translation
riboKcat = 12 # 1/s
riboK0 = 4*25e-6 # mM
riboKd = 0.001 # mM

ribo_init = 5*Ecoli_V*avgdr/60/6800

ribosomeConc = 500*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 200*countToMiliMol # mM

# Global parameter for degradation of proteins
# Derived from eLife's model, using average protein half life of 25 hours. 
ptnDegRate = 7.70e-06 # 1/s


### Load all necessary files
# The reconstruction matches reactions with gene-protein-reactions (GPR) that use MMSYN1* IDs.
reconstPD = pd.read_excel("../model_data/reconstruction.xlsx", sheet_name='Reactions')

# The annotation matches MMSYN1* IDs with JCVISYN3* IDs (or "locus tags").
annotatPD = pd.read_excel("../model_data/FBA/Syn3A_annotation_compilation.xlsx",
                         sheet_name="Syn3A_annotation_compilation_condensed")

# The genome data matches "locus tags" with AOE* protein IDs.
# It provides both the gene sequence, needed for transcription reactions in the ODE model,
# and the protein sequence, needed for translation reactions in the model.
# This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP014992.1).

genomeFile2 = '../model_data/syn2.gb'
genome2 = next(SeqIO.parse(genomeFile2, "gb"))

# This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP016816.2).
genomeFile3A = '../model_data/syn3A.gb'
genome3A = next(SeqIO.parse(genomeFile3A, "gb"))

# The proteomics matches AOE IDs with quantitative proteomics data.
proteomPD = pd.read_excel("../model_data/proteomics.xlsx", sheet_name="Proteomics", skiprows=[0] )

genome_syn3A = list(SeqIO.parse(genomeFile3A, "genbank"))
dna3A = genome_syn3A[0]


AOEtoJ2 = dict()
J2toAOE = dict()
genomeLocDict = dict()
genomePtnLocDict = dict()
genomeRnaLocDict = dict()
Locus3A = []

for f in genome2.features:
    if f.type == "CDS":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
        #print(JCVSYN2_tag)
        # Not all entries have an AOE protein_id
        if('protein_id' in f.qualifiers.keys()):
            AOE_locus = f.qualifiers['protein_id'][0]
            AOEtoJ2[AOE_locus] = JCVSYN2_tag
            J2toAOE[JCVSYN2_tag] = AOE_locus
#             genomeLocDict[JCVSYN2_tag] = f.location
        else:
            print("Locus ", JCVSYN2_tag, " has no AOE id!")
    if f.type == "rRNA":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
#         genomeLocDict[JCVSYN2_tag] = f.location
    if f.type == "tRNA":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
#         genomeLocDict[JCVSYN2_tag] = f.location
        
for f in genome3A.features:
    if f.type == "CDS":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        #print(JCVSYN2_tag)
        # Not all entries have an AOE protein_id
        if('protein_id' in f.qualifiers.keys()):
#             AOE_locus = f.qualifiers['protein_id'][0]
#             AOEtoJ2[AOE_locus] = JCVSYN2_tag
#             J2toAOE[JCVSYN2_tag] = AOE_locus
            genomePtnLocDict[JCVSYN3A_tag] = f.location
            genomeLocDict[JCVSYN3A_tag] = f.location
        else:
            print("Locus ", JCVSYN3A_tag, " is pseudo.")
    if f.type == "rRNA":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        genomeRnaLocDict[JCVSYN3A_tag] = f.location
        genomeLocDict[JCVSYN3A_tag] = f.location
    if f.type == "tRNA":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        genomeRnaLocDict[JCVSYN3A_tag] = f.location
        genomeLocDict[JCVSYN3A_tag] = f.location
        

def getSequences(jcvi3AID):
    # returns genomic and protein sequences
    try:
        rnasequence = genomeLocDict[jcvi3AID].extract(genome3A.seq).transcribe()
        
        # Using translation table 4 from NCBI: "Mycoplasma Code"
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4
        aasequence  = genomeLocDict[jcvi3AID].extract(genome3A.seq).transcribe().translate(table=4)
        
    except:
        aasequence  = 0
        rnasequence = 0
    
    return rnasequence, aasequence


def getRNAsequences(jcvi3AID):
    # returns genomic and protein sequences
    try:
        rnasequence = genomeLocDict[jcvi3AID].extract(genome3A.seq).transcribe()
        
    except:
        rnasequence = 0
    
    return rnasequence
        
# Create list of proteins with no proteomics data
# ptnNoQuant = set()

def getPtnCount(newMetID, jcvi2ID):
    
    # Check if protein quantification is available.
    try:
        if jcvi2ID.startswith("JCVIman_"):
            aoeID = manGPRPD.loc[ manGPRPD.MM == jcvi2ID.replace("JCVIman_",""), "AOE" ].values[0]
        else:
            aoeID = J2toAOE[ jcvi2ID ]
        
        ptnCount = max(defaultPtnCount,round(proteomPD.loc[ proteomPD.Protein == aoeID ].iloc[0,21]))
#         
        ptnName  = proteomPD.loc[ proteomPD.Protein == aoeID ].iloc[0,1].replace(
            " [synthetic bacterium JCVI-Syn3.0]","")
        
#         ptnConcentration = ptnCount*countToMiliMol
    except:
#         print("WARNING: No protein count for", newMetID)
#         print("Using default protein concentration.")

        ptnName = newMetID
        ptnCount = defaultPtnCount
#         ptnConcentration = defaultPtnConcentration

#         ptnNoQuant.add(newMetID)
    
    return ptnCount, ptnName

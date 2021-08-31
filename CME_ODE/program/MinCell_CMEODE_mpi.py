#!/usr/bin/env python
# coding: utf-8


# <h2> Set up gene expression reactions for lipid module </h2>


# Import needed modules
from pyLM import *
from pyLM.units import *
import math as math
import numpy as np

import os
from contextlib import redirect_stdout

# SBtab classes source code
from sbtab import SBtab

# SBtab validator
from sbtab import validatorSBtab

# Converter SBtab -> SBML
from sbtab import sbtab2sbml

# Converter SBML -> SBtab
from sbtab import sbml2sbtab

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import in_out as in_out

import time as timer

# Argument parsing for parallel runs
import argparse
ap = argparse.ArgumentParser()

# The process ID
ap.add_argument('-procid','--processid',required=True)
ap.add_argument('-t','--time',required=True)

args = ap.parse_args()
procid=str(args.processid)

# Use the process ID to seed the random number generation
# so that different random numbers are used each replicate
proc_id = int(procid)
currTime = int(timer.time())
np.random.seed(proc_id*429496)

# ## Load in model data such as the genome and proteomics


### Load all necessary files

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
genome3A_DNA = list(SeqIO.parse(genomeFile3A, "genbank"))

# The proteomics matches AOE IDs with quantitative proteomics data.
proteomPD = pd.read_excel("../model_data/proteomics.xlsx", sheet_name="Proteomics", skiprows=[0] )
# _fba_250

# ## Create Definitions to Get Protein Counts and Gene Sequences
# Create list of proteins with no proteomics data
ptnNoQuant = set()

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
        
    except:
        #print("WARNING: No protein count for", newMetID)
        #print("Using default protein concentration.")

        ptnName = newMetID
        ptnCount = defaultPtnCount

        ptnNoQuant.add(newMetID)
    
    return ptnCount, ptnName



def getSequences(jcvi3AID):
    # returns genomic and protein sequences
    try:
        rnasequence = genomePtnLocDict[jcvi3AID].extract(genome3A.seq).transcribe()
        
        # Using translation table 11 from NCBI: "Bacterial, Archaeal and Plant Plastid Code"
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4
        aasequence  = genomePtnLocDict[jcvi3AID].extract(genome3A.seq).transcribe().translate(table=4)
        
    except:
        aasequence  = 0
        rnasequence = 0
    
    return rnasequence, aasequence




def getRNAsequences(jcvi3AID):
    # returns genomic and protein sequences
    try:
        rnasequence = genomeRnaLocDict[jcvi3AID].extract(genome3A.seq).transcribe()
        
    except:
        rnasequence = 0
    
    return rnasequence




AOEtoJ2 = dict()
J2toAOE = dict()
genomePtnLocDict = dict()
genomeRnaLocDict = dict()
Locus3A = []

for f in genome2.features:
    if f.type == "CDS":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
        # Not all entries have an AOE protein_id
        if('protein_id' in f.qualifiers.keys()):
            AOE_locus = f.qualifiers['protein_id'][0]
            AOEtoJ2[AOE_locus] = JCVSYN2_tag
            J2toAOE[JCVSYN2_tag] = AOE_locus
        else:
            print("Locus ", JCVSYN2_tag, " has no AOE id!")
    if f.type == "rRNA":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
    if f.type == "tRNA":
        JCVSYN2_tag = f.qualifiers['locus_tag'][0]
        
for f in genome3A.features:
    if f.type == "CDS":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        # Not all entries have an AOE protein_id
        if('protein_id' in f.qualifiers.keys()):
            genomePtnLocDict[JCVSYN3A_tag] = f.location
        else:
            print("Locus ", JCVSYN3A_tag, " is pseudo.")
    if f.type == "rRNA":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        genomeRnaLocDict[JCVSYN3A_tag] = f.location
    if f.type == "tRNA":
        JCVSYN3A_tag = f.qualifiers['locus_tag'][0]
        Locus3A.append(JCVSYN3A_tag)
        genomeRnaLocDict[JCVSYN3A_tag] = f.location



#genomePtnLocDict


### Initialize the Simulation


# Create our simulation object
sim=CME.CMESimulation(name="Min Cell Genetic Processes")


### Load Lists of Proteins for Metabolism and with Specific Naming in the model


PtnMetDF = pd.read_csv("../model_data/protein_metabolites_frac.csv")
PtnMetDF


memPtnMetDF = pd.read_csv("../model_data/membrane_protein_metabolites.csv")
memPtnMetDF

riboPtnMetDF = pd.read_csv("../model_data/ribo_protein_metabolites.csv")
riboPtnMetDF


RNApolPtnMetDF = pd.read_csv("../model_data/RNApol_proteins.csv")
RNApolPtnMetDF



rrnaMetDF_1 = pd.read_csv("../model_data/rrna_metabolites_1.csv")
rrnaMetDF_1



rrnaMetDF_2 = pd.read_csv("../model_data/rrna_metabolites_2.csv")
rrnaMetDF_2



trnaMetDF = pd.read_csv("../model_data/trna_metabolites_synthase.csv")
trnaMetDF



MetPtnGenes = ['MMSYN1_0621', 'MMSYN1_0419',
              'MMSYN1_0513', 'MMSYN1_0512','MMSYN1_0117','MMSYN1_0139',
              'MMSYN1_0147','MMSYN1_0304','MMSYN1_0420','MMSYN1_0616','MMSYN1_0617',
              'MMSYN1_0218','MMSYN1_0214','MMSYN1_0875',
               'MMSYN1_0836','MMSYN1_0642','MMSYN1_0643','MMSYN1_0641'] # Includes COAabc genes

MetLocusNums = []

for gene in MetPtnGenes:
    locusNum = gene.split('_')[1]
    MetLocusNums.append(locusNum)
    

named_PTN_list = []

for index, row in riboPtnMetDF.iterrows():
    named_PTN_list.append(row["gene"]) 

    
for index, row in PtnMetDF.iterrows():
    named_PTN_list.append(row["gene"])


### Define Transcription and Translation Rate Constants


ModelSpecies = []


defaultPtnCount = 10


baseMap = OrderedDict({ "A":"M_atp_c", "U":"M_utp_c", "G":"M_gtp_c", "C":"M_ctp_c" })

# Global parameters for transcription
rnaPolKcat = 0.155*187/493*20 # nt/s
rnaPolK0 = 1e-4 #mM
rnaPolKd = 0.1 #0.01 #mM

rrnaPolKcat = 85*2 # nt/s
trnaPolKcat = 0.155*187/493*25 # nt/s

krnadeg = 0.00578/2 # 1/s
ptnDegRate = 7.70e-06 # 1/s

ATPconc = 1.04 #mM
UTPconc = 0.68 #mM
CTPconc = 0.34 #mM
GTPconc = 0.68 #mM

# Cell radius (meters):
r_cell = 2.0*(10**-7) # m

CytoVolume = (4*np.pi/3)*1000*r_cell**3 # L
cellVolume = CytoVolume

# Avogadro:
avgdr   = 6.022e23 # molec/mol
Avognum = avgdr

countToMiliMol = 1000/(avgdr*cellVolume)

RnaPconc = 187*countToMiliMol # mM

# Global parameter for degradation of mRNAs
rnaDegRate = 0.00578/2 # 1/s

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

aaCostMap = OrderedDict({"A":"ALA_cost","R":"ARG_cost","N":"ASN_cost","D":"ASP_cost","C":"CYS_cost",
                         "E":"GLU_cost","Q":"GLN_cost","G":"GLY_cost","H":"HIS_cost","I":"ILE_cost",
                         "L":"LEU_cost","K":"LYS_cost","M":"MET_cost","F":"PHE_cost","P":"PRO_cost",
                         "S":"SER_cost","T":"THR_cost","W":"TRP_cost","Y":"TYR_cost","V":"VAL_cost",
                         "FM":"FMET_cost"})

ctRNAcostMap = OrderedDict({"M_alatrna_c":"ALA_cost","M_argtrna_c":"ARG_cost","M_asntrna_c":"ASN_cost",
                            "M_asptrna_c":"ASP_cost","M_cystrna_c":"CYS_cost","M_glutrna_c":"GLU_cost",
                            "M_glntrna_c":"GLN_cost","M_glytrna_c":"GLY_cost","M_histrna_c":"HIS_cost",
                            "M_iletrna_c":"ILE_cost","M_leutrna_c":"LEU_cost","M_lystrna_c":"LYS_cost",
                            "M_mettrna_c":"MET_cost","M_phetrna_c":"PHE_cost","M_protrna_c":"PRO_cost",
                            "M_sertrna_c":"SER_cost","M_thrtrna_c":"THR_cost","M_trptrna_c":"TRP_cost",
                            "M_tyrtrna_c":"TYR_cost","M_valtrna_c":"VAL_cost"}) #,"FM":"FMET_cost"})

# Global parameters for translation
riboKcat = 10 # 1/s
riboK0 = 4*25e-6 # 
riboKd = 0.0001 # 

ribosomeConc = 503*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 150*countToMiliMol # mM

# Global parameter for degradation of proteins
# Derived from eLife's model, using average protein half life of 25 hours. 
ptnDegRate = 7.70e-06 # 1/s


# In[ ]:


# Define how to calculate transcription rate constants as in equation 3 for transcription reactions.
# Uses mature transcript length and proteomics for promoter strength.

def TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    kcat_mod = min(rnaPolKcat*(ptnCount/(180)),2*90)

    
    # Add total number of monomers to parameter dict
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["C"]
    
    NMono_C = baseCount["G"]
    
    NMono_G = baseCount["U"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    

    k_transcription = kcat_mod / ((1+rnaPolK0/RnaPconc)*(rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
    return k_transcription



# Define transcription rate for ribosomal protein-coding genes.

def riboTranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    
    kcat_mod = rnaPolKcat*500.0/180.0#*(339/(180))
    
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["C"]
    
    NMono_C = baseCount["G"]
    
    NMono_G = baseCount["U"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc

    k_transcription = kcat_mod / ((1+rnaPolK0/RnaPconc)*(rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
    return k_transcription




def DegradationRate(rnaMetID, rnasequence):
    # Define how to calculate transcription rate constants as in equation 3 for transcription reactions.
# Uses mature transcript length and proteomics for promoter strength.
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    kcat = (18/452)*88 #1/s # INSTEAD OF 18 or 20

    n_tot = sum(list(baseCount.values()))

    k_deg = kcat / n_tot 
    
    return k_deg



# Define transcription rate for tRNA genes.

def trnaTranscriptRate(rnaMetID, rnasequence):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["C"]
    
    NMono_C = baseCount["G"]
    
    NMono_G = baseCount["U"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    
    k_transcription = trnaPolKcat/ ((1+rnaPolK0/RnaPconc)*(rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
    return k_transcription


# Define transcription rate for tRNA genes.

def rrnaTranscriptRate(rnasequence):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    n_tot = sum(list(baseCount.values()))
    
    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["C"]
    
    NMono_C = baseCount["G"]
    
    NMono_G = baseCount["U"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    
    k_transcription = rrnaPolKcat / ((1+rnaPolK0/RnaPconc)*(rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
    return k_transcription



from translation_rate_start import TranslatRate as TranslatRate

# Define how to calculate translation rate constants as in equation 3 for translation reactions.

def TranslocRate(aasequence):
    """
    Calculate a translocation rate for membrane proteins based on protein amino acid length and 
    a typical secY translocation rate

    Parameters:

        aasequence (string) - the amino acid sequence of the membrane protein

    Returns:
    
        k_translocation (float) - the translocation rate for the protein

    """

    # Check that we know all residues used in the sequence
    if ( set(aasequence) - set(aaMap.keys()) ):
        raise Exception("Unknown residue(s) in Protein sequence {}".format(set(aasequence) - set(aaMap.keys())) )
    
    # Count how many times each residue is used
    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)
    
    ptnLen = sum(list(aaCount.values()))
    
    k_transloc = 50/ptnLen #secyNum*
    
    return k_transloc


# ## Create Definitions to Add RNA, Proteins, and Reactions to the Simulation


def addPtn(jcvi3AID):
    locusNum = jcvi3AID.split('_')[1]
    mmcode = 'MMSYN1_' + locusNum
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode

    #print(mmcode, jcvi2ID, jcvi3AID)
    
    genes_in_model.append(jcvi3AID)

    ptnMetID = 'M_PTN_' + jcvi3AID + '_c'
    # If the protein is not in the model, add it:

    ModelSpecies.append(ptnMetID)

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)


    geneMetID = jcvi3AID + '_gene'
    
    ModelSpecies.append(geneMetID)

    # Get nucleotide and amino acid sequences, if available
    rnasequence, aasequence = getSequences(jcvi3AID)

    if (rnasequence != 0) and (aasequence != 0):

        rnaMetID = "M_RNA_" + jcvi3AID + "_c"
        rnaName = "(mRNA) " + ptnName

        ModelSpecies.append(rnaMetID)
    
        species = []
        species = [geneMetID, rnaMetID, ptnMetID]
        sim.defineSpecies(species)
        
        for index,row in mRNA_counts_DF.iterrows():
        
            if row["LocusTag"] == jcvi3AID:

                avg_mRNA_cnt = row["Count"]

                if avg_mRNA_cnt == 0.0:
                    avg_mRNA_cnt = 0.001

                init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
                continue

        sim.addParticles(species = geneMetID, count = 1)
        sim.addParticles(species = ptnMetID,  count = int(ptnCount))
        sim.addParticles(species = rnaMetID, count = init_mRNA_count)
        
        TrscProd = [geneMetID,rnaMetID]
        mRNAdegProd = []
    
        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
            mRNAdegProd.append('ATP_mRNAdeg')
            
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
            mRNAdegProd.append('AMP_mRNAdeg')
            
        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
            mRNAdegProd.append('UMP_mRNAdeg')
            
        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
            mRNAdegProd.append('GMP_mRNAdeg')
            
        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
            mRNAdegProd.append('CMP_mRNAdeg')
            
        TranslatProd = [rnaMetID,ptnMetID]
        ptnDegProd = []
        
        for i in range(len(aasequence)):
            TranslatProd.append('ATP_translat')
            TranslatProd.append('ATP_translat')
            
        aaCount = defaultdict(int)
        for aa in set(aasequence):
            aaCount[aa] = aasequence.count(aa)

        NMono_A = aaCount["A"]
        NMono_R = aaCount["R"]
        NMono_N = aaCount["N"]
        NMono_D = aaCount["D"]
        NMono_C = aaCount["C"]
        NMono_E = aaCount["E"]
        NMono_Q = aaCount["Q"]
        NMono_H = aaCount["H"]
        NMono_I = aaCount["I"]
        NMono_L = aaCount["L"]
        NMono_K = aaCount["K"]
        NMono_M = max(0,aaCount["M"] - 1)
        NMono_P = aaCount["P"]
        NMono_S = aaCount["S"]
        NMono_T = aaCount["T"]
        NMono_W = aaCount["W"]
        NMono_Y = aaCount["Y"]
        NMono_G = aaCount["G"]
        NMono_F = aaCount["F"]
        NMono_V = aaCount["V"]
        NMono_FM = 1

        NStop = aaCount["*"]

        if NStop > 1:
            print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

        AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
                     ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
                     ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
                     ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
                     ['FM',NMono_FM]]
        
        for aaCost in AaUsed:
            aa_ID = aaCost[0]
            numberUsed = aaCost[1]
            
            aaCostID = aaCostMap[aa_ID]
            
            for i in range(numberUsed):
                TranslatProd.append(aaCostID)
    
        sim.addReaction(geneMetID, tuple(TrscProd), TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID))
        sim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
        sim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence))
        


def addNamedPtn(newMetID, mmcode, transcribe, ptnFrac):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
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


    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)
    
    ptnCount = ptnFrac*ptnCount
    
    ModelSpecies.append(newMetID)

    
    if transcribe:

        geneMetID = jcvi3AID + '_gene'

        ModelSpecies.append(geneMetID)

        rnasequence, aasequence = getSequences(jcvi3AID)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName

        ModelSpecies.append(rnaMetID)

        species = []
        species = [geneMetID, rnaMetID, newMetID]
        sim.defineSpecies(species)
        
        for index,row in mRNA_counts_DF.iterrows():

            if row["LocusTag"] == jcvi3AID:

                avg_mRNA_cnt = row["Count"]

                if avg_mRNA_cnt == 0.0:
                    avg_mRNA_cnt = 0.001

                init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
                continue


        sim.addParticles(species = geneMetID, count = 1)
        sim.addParticles(species = newMetID,  count = int(ptnCount))
        sim.addParticles(species = rnaMetID, count = init_mRNA_count)

        TrscProd = [geneMetID,rnaMetID]
        mRNAdegProd = []

        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
            mRNAdegProd.append('ATP_mRNAdeg')

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
            mRNAdegProd.append('AMP_mRNAdeg')

        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
            mRNAdegProd.append('UMP_mRNAdeg')

        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
            mRNAdegProd.append('GMP_mRNAdeg')

        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
            mRNAdegProd.append('CMP_mRNAdeg')

        TranslatProd = [rnaMetID,newMetID]
        ptnDegProd = []

        for i in range(len(aasequence)):
            TranslatProd.append('ATP_translat')
            TranslatProd.append('ATP_translat')
            
        aaCount = defaultdict(int)
        for aa in set(aasequence):
            aaCount[aa] = aasequence.count(aa)

        NMono_A = aaCount["A"]
        NMono_R = aaCount["R"]
        NMono_N = aaCount["N"]
        NMono_D = aaCount["D"]
        NMono_C = aaCount["C"]
        NMono_E = aaCount["E"]
        NMono_Q = aaCount["Q"]
        NMono_H = aaCount["H"]
        NMono_I = aaCount["I"]
        NMono_L = aaCount["L"]
        NMono_K = aaCount["K"]
        NMono_M = max(0,aaCount["M"] - 1)
        NMono_P = aaCount["P"]
        NMono_S = aaCount["S"]
        NMono_T = aaCount["T"]
        NMono_W = aaCount["W"]
        NMono_Y = aaCount["Y"]
        NMono_G = aaCount["G"]
        NMono_F = aaCount["F"]
        NMono_V = aaCount["V"]
        NMono_FM = 1

        NStop = aaCount["*"]

        if NStop > 1:
            print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

        AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
                     ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
                     ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
                     ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
                     ['FM',NMono_FM]]

        for aaCost in AaUsed:
            aa_ID = aaCost[0]
            numberUsed = aaCost[1]

            aaCostID = aaCostMap[aa_ID]

            for i in range(numberUsed):
                TranslatProd.append(aaCostID)

        sim.addReaction(geneMetID, tuple(TrscProd), TranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID))
        sim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
        sim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))

    
    if not transcribe:
        species = [newMetID]
        sim.defineSpecies(species)
        sim.addParticles(species = newMetID,  count = int(ptnCount))




def addMembranePtn(newMetID, mmcode, transcribe, ptnFrac):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
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


    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)
    
    ptnCount = ptnFrac*ptnCount
    
    ModelSpecies.append(newMetID)

    
    if transcribe:

        geneMetID = jcvi3AID + '_gene'

        ModelSpecies.append(geneMetID)

        rnasequence, aasequence = getSequences(jcvi3AID)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName

        ModelSpecies.append(rnaMetID)
        
        newMetID_cyto = newMetID + '_cyto'

        species = []
        species = [geneMetID, rnaMetID, newMetID, newMetID_cyto]
        sim.defineSpecies(species)
        
        for index,row in mRNA_counts_DF.iterrows():

            if row["LocusTag"] == jcvi3AID:

                avg_mRNA_cnt = row["Count"]

                if avg_mRNA_cnt == 0.0:
                    avg_mRNA_cnt = 0.001

                init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
                continue


        sim.addParticles(species = geneMetID, count = 1)
        sim.addParticles(species = newMetID,  count = int(ptnCount))
        sim.addParticles(species = rnaMetID, count = init_mRNA_count)

        TrscProd = [geneMetID,rnaMetID]
        mRNAdegProd = []

        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
            mRNAdegProd.append('ATP_mRNAdeg')

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
            mRNAdegProd.append('AMP_mRNAdeg')

        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
            mRNAdegProd.append('UMP_mRNAdeg')

        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
            mRNAdegProd.append('GMP_mRNAdeg')

        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
            mRNAdegProd.append('CMP_mRNAdeg')

        TranslatProd = [rnaMetID,newMetID_cyto]
        TranslocProd = [newMetID]
        ptnDegProd = []

        for i in range(len(aasequence)):
            TranslatProd.append('ATP_translat')
            TranslatProd.append('ATP_translat')
    
        for i in range(int(len(aasequence)/10)):
            TranslocProd.append('ATP_transloc')
            
        aaCount = defaultdict(int)
        for aa in set(aasequence):
            aaCount[aa] = aasequence.count(aa)

        NMono_A = aaCount["A"]
        NMono_R = aaCount["R"]
        NMono_N = aaCount["N"]
        NMono_D = aaCount["D"]
        NMono_C = aaCount["C"]
        NMono_E = aaCount["E"]
        NMono_Q = aaCount["Q"]
        NMono_H = aaCount["H"]
        NMono_I = aaCount["I"]
        NMono_L = aaCount["L"]
        NMono_K = aaCount["K"]
        NMono_M = max(0,aaCount["M"] - 1)
        NMono_P = aaCount["P"]
        NMono_S = aaCount["S"]
        NMono_T = aaCount["T"]
        NMono_W = aaCount["W"]
        NMono_Y = aaCount["Y"]
        NMono_G = aaCount["G"]
        NMono_F = aaCount["F"]
        NMono_V = aaCount["V"]
        NMono_FM = 1

        NStop = aaCount["*"]

        if NStop > 1:
            print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

        AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
                     ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
                     ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
                     ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
                     ['FM',NMono_FM]]

        for aaCost in AaUsed:
            aa_ID = aaCost[0]
            numberUsed = aaCost[1]

            aaCostID = aaCostMap[aa_ID]

            for i in range(numberUsed):
                TranslatProd.append(aaCostID)

        sim.addReaction(geneMetID, tuple(TrscProd), TranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID))
        sim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
        sim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))
        sim.addReaction(newMetID_cyto, tuple(TranslocProd), TranslocRate(aasequence))
    
    if not transcribe:
        species = [newMetID]
        sim.defineSpecies(species)
        sim.addParticles(species = newMetID,  count = int(ptnCount))



RiboPtnNames = ['Ribosomal Protein']
RiboPtnLens = ['Gene Length (nt)']
RiboPtnTrscRates = ['Transcription Rate']
RiboPtnTranslatRates = ['Translation Rate']

def addRiboPtn(newMetID, mmcode):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum

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
            unkIter += 1 


    ptnCount, ptnName = getPtnCount(newMetID, jcvi2ID)

    
    ModelSpecies.append(newMetID)

    geneMetID = jcvi3AID + '_gene'

    ModelSpecies.append(geneMetID)
    
    rnasequence, aasequence = getSequences(jcvi3AID)

    if (rnasequence != 0) and (aasequence != 0):

        rnaMetID = "M_RNA_" + jcvi3AID + "_c"
        rnaName = "(mRNA) " + ptnName
        
    ModelSpecies.append(rnaMetID)

    species = []
    species = [geneMetID, rnaMetID, newMetID]
    sim.defineSpecies(species)
    
    for index,row in mRNA_counts_DF.iterrows():
        
        if row["LocusTag"] == jcvi3AID:
            
            avg_mRNA_cnt = row["Count"]
            
            if avg_mRNA_cnt == 0.0:
                avg_mRNA_cnt = 0.001
            
            init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
            
            continue

    sim.addParticles(species = geneMetID, count = 1)
    sim.addParticles(species = newMetID,  count = int(ptnCount))
    sim.addParticles(species = rnaMetID, count = init_mRNA_count)
    
    RiboPtnNames.append(newMetID)
    
    trsc_rate = riboTranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID)
    RiboPtnTrscRates.append(trsc_rate)
    
    translat_rate = TranslatRate(rnaMetID, newMetID, rnasequence, aasequence)
    RiboPtnTranslatRates.append(translat_rate)
    
    rna_length = len(rnasequence)
    RiboPtnLens.append(rna_length)
    
    TrscProd = [geneMetID,rnaMetID]
    mRNAdegProd = []
    
    for i in range(len(rnasequence)):
        TrscProd.append('ATP_trsc')
        mRNAdegProd.append('ATP_mRNAdeg')
        
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
        mRNAdegProd.append('AMP_mRNAdeg')

    for i in range(N_U):
        TrscProd.append('UTP_mRNA')
        mRNAdegProd.append('UMP_mRNAdeg')

    for i in range(N_G):
        TrscProd.append('GTP_mRNA')
        mRNAdegProd.append('GMP_mRNAdeg')

    for i in range(N_C):
        TrscProd.append('CTP_mRNA')
        mRNAdegProd.append('CMP_mRNAdeg')
        
    TranslatProd = [rnaMetID,newMetID]
        
    for i in range(len(aasequence)):
        TranslatProd.append('ATP_translat')
        TranslatProd.append('ATP_translat')
            
    aaCount = defaultdict(int)
    for aa in set(aasequence):
        aaCount[aa] = aasequence.count(aa)

    NMono_A = aaCount["A"]
    NMono_R = aaCount["R"]
    NMono_N = aaCount["N"]
    NMono_D = aaCount["D"]
    NMono_C = aaCount["C"]
    NMono_E = aaCount["E"]
    NMono_Q = aaCount["Q"]
    NMono_H = aaCount["H"]
    NMono_I = aaCount["I"]
    NMono_L = aaCount["L"]
    NMono_K = aaCount["K"]
    NMono_M = max(0,aaCount["M"] - 1)
    NMono_P = aaCount["P"]
    NMono_S = aaCount["S"]
    NMono_T = aaCount["T"]
    NMono_W = aaCount["W"]
    NMono_Y = aaCount["Y"]
    NMono_G = aaCount["G"]
    NMono_F = aaCount["F"]
    NMono_V = aaCount["V"]
    NMono_FM = 1

    NStop = aaCount["*"]

    if NStop > 1:
        print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

    AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
                 ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
                 ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
                 ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
                 ['FM',NMono_FM]]

    for aaCost in AaUsed:
        aa_ID = aaCost[0]
        numberUsed = aaCost[1]

        aaCostID = aaCostMap[aa_ID]

        for i in range(numberUsed):
            TranslatProd.append(aaCostID)
    
    sim.addReaction(geneMetID, tuple(TrscProd), riboTranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID))
    sim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
    sim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))




tRNAadded = []

def addtRNA(unchargedMetID, chargedMetID, mmcode, aminoAcid, synthase):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum

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


    geneMetID = jcvi3AID + '_gene'

    ModelSpecies.append(geneMetID)

    rnasequence = getRNAsequences(jcvi3AID)

    gene = [geneMetID]
    sim.defineSpecies(gene)

    sim.addParticles(species = geneMetID, count = 1)


    if unchargedMetID not in tRNAadded:
        ModelSpecies.append(unchargedMetID)
        tRNA = [unchargedMetID]
        sim.defineSpecies(tRNA)
        tRNAadded.append(unchargedMetID)
        ModelSpecies.append(chargedMetID)
        ctRNA = [chargedMetID]
        sim.defineSpecies(ctRNA)
        tRNAadded.append(chargedMetID)
    
        if 'GLN' not in synthase:
    
            synthase_atp = synthase + '_ATP'
            synth_atp = [synthase_atp]
            sim.defineSpecies(synth_atp)
            sim.addParticles(species = synthase_atp, count = 1)
            
            synthase_atp_aa = synthase_atp + '_AA'
            synth_atp_aa = [synthase_atp_aa]
            sim.defineSpecies(synth_atp_aa)
            sim.addParticles(species = synthase_atp_aa, count = 1)
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNA'
            synth_atp_aa_trna = [synthase_atp_aa_trna]
            sim.defineSpecies(synth_atp_aa_trna)
            sim.addParticles(species = synthase_atp_aa_trna, count = 1)
            
            synthase_ptn = 'M_PTN_' + synthase + '_c'

            sim.addReaction(('M_atp_c',synthase_ptn),(synthase_atp),0.01)
            sim.addReaction((aminoAcid,synthase_atp),(synthase_atp_aa),0.01)
            sim.addReaction((synthase_atp_aa,unchargedMetID),(synthase_atp_aa_trna),0.1)
            sim.addReaction((synthase_atp_aa_trna),(chargedMetID,'M_amp_c','M_ppi_c',synthase_ptn),30)
            
            costID = ctRNAcostMap[chargedMetID]
            
            cost_paid = cost + '_paid'
            
            sim.addReaction((chargedMetID,costID),(unchargedMetID,cost_paid),100)
            
            
        elif 'GLN' in synthase:
            synthase = 'JCVISYN3A_0126'
            synthase_atp = synthase + '_ATP'
            synthase_atp_aa = synthase_atp + '_AA'
            synthase_ptn = 'M_PTN_' + synthase + '_c'

            synthase_atp_aa_trna = synthase_atp_aa + '_tRNAgln'
            synth_atp_aa_trna = [synthase_atp_aa_trna]
            sim.defineSpecies(synth_atp_aa_trna)
            sim.addParticles(species = synthase_atp_aa_trna, count = 1)

            sim.addReaction((synthase_atp_aa,'M_trnagln_c'),(synthase_atp_aa_trna),0.1)
            sim.addReaction((synthase_atp_aa_trna),('M_glutrnagln_c','M_amp_c','M_ppi_c',synthase_ptn),30)
            glugln = 'M_glutrnagln_c'
            glugln_enz = 'glutrnagln_enz'
            glugln_enz_atp = 'glutrnagln_enz_atp'
            glugln_enz_atp_gln = 'glutrnagln_enz_atp_gln'
            glugln_enz_atp_asn = 'glutrnagln_enz_atp_asn'
            gge = [glugln_enz]
            sim.defineSpecies(gge)
            sim.addParticles(species = glugln_enz, count = 1)
            ggea = [glugln_enz_atp]
            sim.defineSpecies(ggea)
            sim.addParticles(species = glugln_enz_atp, count = 1)
            ggeagln = [glugln_enz_atp_gln]
            sim.defineSpecies(ggeagln)
            sim.addParticles(species = glugln_enz_atp_gln, count = 1)
            ggeaasn = [glugln_enz_atp_asn]
            sim.defineSpecies(ggeaasn)
            sim.addParticles(species = glugln_enz_atp_asn, count = 1)

            sim.addReaction((glugln,'M_PTN_JCVISYN3A_0689_c'),(glugln_enz),0.01)
            sim.addReaction((glugln_enz,'M_atp_c'),(glugln_enz_atp),0.01)
            sim.addReaction((glugln_enz_atp,'M_gln__L_c'),(glugln_enz_atp_gln),0.01)
            sim.addReaction((glugln_enz_atp_gln),('M_glntrna_c','M_glu__L_c','M_adp_c','M_pi_c','M_PTN_JCVISYN3A_0689_c'),30)
            sim.addReaction((glugln_enz_atp,'M_asn__L_c'),(glugln_enz_atp_asn),0.001)
            sim.addReaction((glugln_enz_atp_asn),('M_glntrna_c','M_asp__L_c','M_adp_c','M_pi_c','M_PTN_JCVISYN3A_0689_c'),30)
        
            
            costID = ctRNAcostMap['M_glntrna_c']
            
            cost_paid = cost + '_paid'
            
            sim.addReaction(('M_glntrna_c',costID),('M_trnagln_c',cost_paid),100)
            
    
    sim.addParticles(species = unchargedMetID,  count = int(round(5800/29*0.2)))
    sim.addParticles(species = chargedMetID,  count = int(round(5800/29*0.8)))
    TrscProd = [geneMetID,unchargedMetID]
        
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
        
    
    sim.addReaction(geneMetID, tuple(TrscProd), trnaTranscriptRate(unchargedMetID, rnasequence))




def addrRNA():
    
    rrna_gene_locs_1 = []
    rrna_gene_locs_2 = []

    rRNA_species = []
    
    for index, row in rrnaMetDF_1.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
        gene_str = str(genomeRnaLocDict[jcvi3AID])
        start_nt = int(gene_str.split(':')[0].replace('[',''))
        end_nt = int(gene_str.split(':')[1].replace('](-)',''))
        rrna_gene_locs_1.append(start_nt)
        rrna_gene_locs_1.append(end_nt)
            
        rRNA_species.append(newMetID)
        
    rrna_pos_1 = Seq(str(genome3A.seq[min(rrna_gene_locs_1):max(rrna_gene_locs_1)]))
    
    rrna_operon_1 = rrna_pos_1.reverse_complement()
            
    
    
    for index, row in rrnaMetDF_2.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
        gene_str = str(genomeRnaLocDict[jcvi3AID])
        start_nt = int(gene_str.split(':')[0].replace('[',''))
        end_nt = int(gene_str.split(':')[1].replace('](-)',''))
        rrna_gene_locs_2.append(start_nt)
        rrna_gene_locs_2.append(end_nt)
        
    rrna_pos_2 = Seq(str(genome3A.seq[min(rrna_gene_locs_2):max(rrna_gene_locs_2)]))
    
    rrna_operon_2 = rrna_pos_2.reverse_complement()
            
    
    sim.defineSpecies(rRNA_species)
    
    
    for rRNA in rRNA_species:
        
        sim.addParticles(species = rRNA,  count = 1)
            
        ModelSpecies.append(rRNA)
    
    for index, row in rrnaMetDF_1.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
        geneMetID = jcvi3AID + '_gene'
    
        ModelSpecies.append(geneMetID)
        
        gene = [geneMetID]
        sim.defineSpecies(gene)

        sim.addParticles(species = geneMetID, count = 1)


    rnasequence = rrna_operon_1.transcribe()
        
    TrscProd = ['JCVISYN3A_0067_gene','M_rRNA_5S_c','M_rRNA_23S_c','M_rRNA_16S_c']
        
    for i in range(int(len(rnasequence))):
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
        TrscProd.append('ATP_rRNA')

    for i in range(N_U):
        TrscProd.append('UTP_rRNA')

    for i in range(N_G):
        TrscProd.append('GTP_rRNA')

    for i in range(N_C):
        TrscProd.append('CTP_rRNA')

    sim.addReaction('JCVISYN3A_0067_gene', tuple(TrscProd), rrnaTranscriptRate(rnasequence))
        
        
    for index, row in rrnaMetDF_2.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
        geneMetID = jcvi3AID + '_gene'
    
        ModelSpecies.append(geneMetID)
        
        gene = [geneMetID]
        sim.defineSpecies(gene)

        sim.addParticles(species = geneMetID, count = 1)


    rnasequence = rrna_operon_2.transcribe()
        
    TrscProd = ['JCVISYN3A_0532_gene','M_rRNA_5S_c','M_rRNA_23S_c','M_rRNA_16S_c']
        
    for i in range(int(len(rnasequence))):
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
        TrscProd.append('ATP_rRNA')

    for i in range(N_U):
        TrscProd.append('UTP_rRNA')

    for i in range(N_G):
        TrscProd.append('GTP_rRNA')

    for i in range(N_C):
        TrscProd.append('CTP_rRNA')

    sim.addReaction('JCVISYN3A_0532_gene', tuple(TrscProd), rrnaTranscriptRate(rnasequence))
    




NA = 6.022*(10**(23))

from rep_start import addRepInit as addRepInit
from rep_start import addReplication as addReplication

import setICs as setICs
setICs.__main__(sim)



# Populate the model with particles and reactions.

genes_in_model = []

mRNA_counts_DF = pd.read_csv('../model_data/mRNA_counts.csv')

Nuc_counters = ['ATP_trsc','ATP_translat','ATP_mRNAdeg','ATP_ptndeg','ATP_DNArep','ATP_transloc',
               'ATP_mRNA','UTP_mRNA','CTP_mRNA','GTP_mRNA',
               'AMP_mRNAdeg','UMP_mRNAdeg','CMP_mRNAdeg','GMP_mRNAdeg',
               'ATP_tRNA','UTP_tRNA','CTP_tRNA','GTP_tRNA',
               'ATP_rRNA','UTP_rRNA','CTP_rRNA','GTP_rRNA',
               'dATP_DNArep','dTTP_DNArep','dCTP_DNArep','dGTP_DNArep']

sim.defineSpecies(Nuc_counters)

AA_counters = ["ALA_cost","ARG_cost","ASN_cost","ASP_cost","CYS_cost","GLU_cost","GLN_cost","GLY_cost",
               "HIS_cost","ILE_cost","LEU_cost","LYS_cost","MET_cost","PHE_cost","PRO_cost","SER_cost",
               "THR_cost","TRP_cost","TYR_cost","VAL_cost","FMET_cost"]

sim.defineSpecies(AA_counters)

AA_paid = []

for cost in AA_counters:
    
    cost_paid = cost + '_paid'
    AA_paid.append(cost_paid)
    
sim.defineSpecies(AA_paid)

for index, row in PtnMetDF.iterrows():
    addNamedPtn(row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"])
    
for index, row in memPtnMetDF.iterrows():
    addMembranePtn(row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"])

for index, row in riboPtnMetDF.iterrows():
    addRiboPtn(row["species"], row["gene"])
    
for gene in genomePtnLocDict:
    if gene not in genes_in_model:
        addPtn(gene)
        
    
for index, row in trnaMetDF.iterrows():
    addtRNA(row["uncharged"], row["charged"], row["gene"], row["AA"], row["synthase"])
    
addrRNA()

from rep_start.py import addRepInit
from rep_start.py import addReplication

addRepInit(sim)
addReplication(sim,genome3A_DNA,ModelSpecies)



NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell
def mMtoPart(concDict):
    partDict = {}
    for key,val in concDict.items():
        mM = concDict[key]
        particle = (mM/1000)*NA*V
        partDict[key] = round(particle)
    return partDict


fn = 'simulations/batch1/out-' + str(procid) + '.lm'
my_log_file = 'simulations/batch1/log-' + str(procid) + '.log'

try:
    os.remove(fn)

except:
    print('Nothing to delete')



try:
    os.mkdir('fluxes/batch1/'+ 'rep-' + str(procid))#+'fluxes')
    os.mkdir('simulations/batch1/'+'rep-'+str(procid))
    print('Made fluxes folder')
    
except:
    print('Fluxes folder exists')


# In[ ]:


# break
import time
start = time.time()




### Constants
delt = 1.0 # s
odestep = delt/10.0#100.0 # s
write = 60.0 # s
st=int(args.time) #in minutes
simTime = st*60.0 #s

sim.setWriteInterval(write)
sim.setHookInterval(delt)
sim.setSimulationTime(simTime)

sim.save(fn)




"""
File to run the hybrid simulation

Author: David Bianchi
"""
import numpy as np
import hook_start as hook
import sys
import lm as lm
import species_counts as species_counts
from pyLM import *


# Compile with Cython or NO?
cythonBool = False # For smaller system no cython (scipy.ode)


# TODO: Probably not really necessary!!!
IC = np.zeros(25)

mySpecies = species_counts.SpeciesCounts(sim)

totalTime = simTime


with open(my_log_file, 'w') as f, redirect_stdout(f):
    odeHookSolver = hook.MyOwnSolver(IC, delt, odestep, mySpecies, cythonBool, totalTime,str(procid))

    sim.runSolver(filename=fn,solver=odeHookSolver,replicates=1, cudaDevices=None)

    f.close()

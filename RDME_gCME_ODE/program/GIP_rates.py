### Rate Constant Calculations #####

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

##################
# Define how to calculate transcription rate constants for transcription reactions.
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
    
    kcat_mod = min(rnaPolKcat*(ptnCount/(180)),85)

    kcat_mod = max(10,kcat_mod)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    
    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


def TranscriptRateCME(ptnMetID, rnasequence, jcvi2ID, pmap):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    kcat_mod = min(rnaPolKcat*(ptnCount/(180)),85)

    kcat_mod = max(10,kcat_mod)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    atp = partTomM(max(1,pmap['M_atp_c']), pmap)
    ctp = partTomM(max(1,pmap['M_ctp_c']), pmap)
    gtp = partTomM(max(1,pmap['M_gtp_c']), pmap)
    utp = partTomM(max(1,pmap['M_utp_c']), pmap)
    
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


#################
def RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID):
    
    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    if ptnCount>765:
        binding_rate = 4*Ecoli_V*avgdr/11400/60 #/1800/60
    elif ptnCount<45:
        binding_rate = 4*(180/765)*(45/180)*Ecoli_V*avgdr/11400/60 #/1800/60
    else:
        binding_rate = 4*(180/765)*(ptnCount/180)*Ecoli_V*avgdr/11400/60 #/1800/60
    
    return binding_rate
##################

#################
def RNAP_binding_ribo(rnaMetID, ptnMetID, rnasequence, jcvi2ID):
    
    binding_rate = 4*(180/765)*(500/180)*Ecoli_V*avgdr/11400/60 #/1800/60
    
    return binding_rate
##################

##################
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
    
#     kcat_mod = min(rnaPolKcat*(ptnCount/(180)),2*85)
    
    kcat_mod = rnaPolKcat*(500/(180))
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc

    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


##################
# Define transcription rate for ribosomal protein-coding genes.

def riboTranscriptRateCME(ptnMetID, rnasequence, jcvi2ID, pmap):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
#     kcat_mod = min(rnaPolKcat*(ptnCount/(180)),2*85)
    
    kcat_mod = rnaPolKcat*(500/(180))
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    atp = partTomM(max(1,pmap['M_atp_c']), pmap)
    ctp = partTomM(max(1,pmap['M_ctp_c']), pmap)
    gtp = partTomM(max(1,pmap['M_gtp_c']), pmap)
    utp = partTomM(max(1,pmap['M_utp_c']), pmap)
    
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = kcat_mod / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################



##################
# Define transcription rate for tRNA genes.

def trnaTranscriptRate(rnasequence):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    
    k_transcription = rnaPolKcat/ ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


##################
# Define transcription rate for tRNA genes.

def trnaTranscriptRateCME(rnasequence, pmap):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict

    n_tot = sum(list(baseCount.values()))

    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    atp = partTomM(max(1,pmap['M_atp_c']), pmap)
    ctp = partTomM(max(1,pmap['M_ctp_c']), pmap)
    gtp = partTomM(max(1,pmap['M_gtp_c']), pmap)
    utp = partTomM(max(1,pmap['M_utp_c']), pmap)
    
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = rnaPolKcat/ ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


##################
# Define transcription rate for rRNA genes.

def rrnaTranscriptRate(rnasequence):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict
#     paramDict["n"] = sum(list(baseCount.values()))
    n_tot = sum(list(baseCount.values()))
    
    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    
    k_transcription = rrnaPolKcat / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


##################
# Define transcription rate for rRNA genes.

def rrnaTranscriptRateCME(rnasequence, pmap):
    # Add trascription reaction
    
    # Check that we know all bases used in the sequence
    if ( set(rnasequence) - set(baseMap.keys()) ):
        raise Exception("Unknown base(s) in RNA sequence {}".format(set(rnasequence) - set(baseMap.keys())) )
    
    # Count how many times each base is used
    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict
#     paramDict["n"] = sum(list(baseCount.values()))
    n_tot = sum(list(baseCount.values()))
    
    NMono_A = baseCount["A"]
    
    NMono_U = baseCount["U"]
    
    NMono_C = baseCount["C"]
    
    NMono_G = baseCount["G"]
    
    NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    CMono1 = baseMap[ rnasequence[0] ]
    
    CMono2 = baseMap[ rnasequence[1] ]
    
    atp = partTomM(max(1,pmap['M_atp_c']), pmap)
    ctp = partTomM(max(1,pmap['M_ctp_c']), pmap)
    gtp = partTomM(max(1,pmap['M_gtp_c']), pmap)
    utp = partTomM(max(1,pmap['M_utp_c']), pmap)
    
    NMonoSum = NMono_A*rnaPolKd/atp + NMono_C*rnaPolKd/ctp + NMono_U*rnaPolKd/utp + NMono_G*rnaPolKd/gtp
    
    k_transcription = rrnaPolKcat / ((rnaPolKd**2)/(CMono1*CMono2) + NMonoSum + n_tot - 1)
    
#     k_transcription = k_transcription * NaV
    
    return k_transcription
##################


##################
# Define how to calculate translation rate constants as in equation 3 for translation reactions.

def TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence):
    # Add translation reaction
    
    # Considers amino acids up to the first stop codon.
#     aasequence = aasequence[0:aasequence.find("*")]
    
    # Check that we know all residues used in the sequence
    if ( set(aasequence) - set(aaMap.keys()) ):
        raise Exception("Unknown residue(s) in Protein sequence {}".format(set(aasequence) - set(aaMap.keys())) )
    
    # Count how many times each residue is used
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
    NMono_M = aaCount["M"]
    NMono_P = aaCount["P"]
    NMono_S = aaCount["S"]
    NMono_T = aaCount["T"]
    NMono_W = aaCount["W"]
    NMono_Y = aaCount["Y"]
    NMono_G = aaCount["G"]
    NMono_F = aaCount["F"]
    NMono_V = aaCount["V"]
    
    NStop = aaCount["*"]
    
    if NStop > 1:
        print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")
    
    NMonoDict = [NMono_A,NMono_R,NMono_N,NMono_D,NMono_C,NMono_E,NMono_Q,NMono_H,
                 NMono_I,NMono_L,NMono_K,NMono_M,NMono_P,NMono_S,NMono_T,NMono_W,
                 NMono_Y,NMono_G,NMono_F,NMono_V]
    
    NMonoSum = 0
    
    for nmono in range(0,len(NMonoDict)):
        NMonoSum = NMonoSum + NMonoDict[nmono]*riboKd/ctRNAconc
        
    n_tot = sum(list(aaCount.values()))

    baseCount = defaultdict(int)
    for base in set(rnasequence):
        baseCount[base] = rnasequence.count(base)
        
    transcript_length = sum(list(baseCount.values()))
    
#     print(transcript_length)
    
#     ribo_num = max(1,int(transcript_length/300)) #max(1,int(transcript_length/300))
    
#     print(ribo_num)
    
    kcat_mod = riboKcat #*ribo_num #*0.4
    
    k_translation = kcat_mod / ((riboKd**2)/(ctRNAconc**2) + NMonoSum + n_tot - 1)
    
#     k_translation = k_translation * NaV
    
    return k_translation
##################


##################
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
    

#     ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    kcat = 88 #1/s
    
    # The rate form needs specific sequence data for the first two monomers:
#     paramDict["CMonoA"] = baseMap[ rnasequence[0] ]
#     paramDict["CMonoB"] = baseMap[ rnasequence[1] ]
#     paramDict["KDA"] = rnaPolKd # Since we are current;y using the same Kd for all nucleotides
#     paramDict["KDB"] = rnaPolKd 
    
    # Add total number of monomers to parameter dict
    
#     CMono1 = baseMap[ rnasequence[0] ]
    
#     CMono2 = baseMap[ rnasequence[1] ]

    n_tot = sum(list(baseCount.values()))

#     NMono_A = baseCount["A"]
    
#     NMono_U = baseCount["U"]
    
#     NMono_C = baseCount["C"]
    
#     NMono_G = baseCount["G"]
    
#     NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_U]
#     print(NMonoDict)
    
    
#     NMonoSum = NMono_A*rnaPolKd/ATPconc + NMono_C*rnaPolKd/CTPconc + NMono_U*rnaPolKd/UTPconc + NMono_G*rnaPolKd/GTPconc
    

    k_deg = kcat / n_tot 
    
#     k_transcription = k_transcription * NaV
    
    return k_deg
##################


def TranslocRate(aasequence):
    
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
degrad_bind_rate = 11*avgdr*Ecoli_V/60/2400 #1/M/s 7800

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


defaultPtnCount = 1

# Global parameters for translation
riboKcat = 12 # 1/s
riboK0 = 4*25e-6 # mM
riboKd = 0.001 # mM

ribo_init = 30*Ecoli_V*avgdr/60/6800

ribosomeConc = 500*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 200*countToMiliMol # mM

# Global parameter for degradation of proteins
# Derived from eLife's model, using average protein half life of 25 hours. 
ptnDegRate = 7.70e-06 # 1/s


### Load all necessary files
# The reconstruction matches reactions with gene-protein-reactions (GPR) that use MMSYN1* IDs.
reconstPD = pd.read_excel("../../model_data/reconstruction.xlsx", sheet_name='Reactions')

# The annotation matches MMSYN1* IDs with JCVISYN3* IDs (or "locus tags").
annotatPD = pd.read_excel("../../model_data/FBA/Syn3A_annotation_compilation.xlsx",
                         sheet_name="Syn3A_annotation_compilation_condensed")

# The genome data matches "locus tags" with AOE* protein IDs.
# It provides both the gene sequence, needed for transcription reactions in the ODE model,
# and the protein sequence, needed for translation reactions in the model.
# This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP014992.1).

genomeFile2 = '../../model_data/syn2.gb'
genome2 = next(SeqIO.parse(genomeFile2, "gb"))

# This is the NCBI Gene Bank-formated file (https://www.ncbi.nlm.nih.gov/nuccore/CP016816.2).
genomeFile3A = '../../model_data/syn3A.gb'
genome3A = next(SeqIO.parse(genomeFile3A, "gb"))

# The proteomics matches AOE IDs with quantitative proteomics data.
proteomPD = pd.read_excel("../../model_data/proteomics.xlsx", sheet_name="Proteomics", skiprows=[0] )

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


def calcCellVolume(pmap):
    
    SurfaceArea = pmap['CellSA']
    
    cellRadius_calc = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
    cellRadius = min(cellRadius_calc,255e-9)
    
    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
#     print('Volume',cellVolume)
    
    return cellVolume

def partTomM(particles, pmap):
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
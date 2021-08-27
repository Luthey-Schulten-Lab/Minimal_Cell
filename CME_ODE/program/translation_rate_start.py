from collections import defaultdict, OrderedDict
import numpy as np

# Cell radius (meters):
# r_cell = 2.5*(10**-7)
r_cell = 2.0*(10**-7) # m

CytoVolume = (4*np.pi/3)*1000*r_cell**3 # L
cellVolume = CytoVolume

# Avogadro:
avgdr   = 6.022e23 # molec/mol
Avognum = avgdr

countToMiliMol = 1000/(avgdr*cellVolume)

# Global parameters for translation
riboKcat = 12 #10 #10 # 1/s
riboK0 = 4*25e-6 # 
riboKd = 0.0001 # 

ribosomeConc = 503*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 150*countToMiliMol # mM

# Global parameter for degradation of proteins
# Derived from Breuer et al. eLife 2019 model, using average protein half life of 25 hours. 
ptnDegRate = 7.70e-06 # 1/s

# Create Dictionaries to map tRNAs to associated aa abbreviations in protein sequences.
aaMap = OrderedDict({"A":"M_ala__L_c", "R":"M_arg__L_c",
    "N":"M_asn__L_c", "D":"M_asp__L_c", "C":"M_cys__L_c", "E":"M_glu__L_c", "Q":"M_gln__L_c", "G":"M_gly_c",
    "H":"M_his__L_c", "I":"M_ile__L_c", "L":"M_leu__L_c", "K":"M_lys__L_c", "M":"M_met__L_c", "F":"M_phe__L_c",
    "P":"M_pro__L_c", "S":"M_ser__L_c", "T":"M_thr__L_c", "W":"M_trp__L_c", "Y":"M_tyr__L_c", "V":"M_val__L_c",
    "*":"Stop_Codon"})

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
    
    
    ribo_num = max(1,round(transcript_length/125-1)) #max(1,int(transcript_length/300))
    
    ribo_num = min(15,ribo_num)
    
#     poly_size_dict = {'3': 75, '2': 120, '4': 32, '6': 18, '9': 9, '5': 5, '7': 7, '12': 12}
#     poly_size_dict = {'3': 105, '5': 35, '2': 118, '4': 56, '10': 10, '6': 48, '7': 14}
#     poly_size_dict = {'3': 81, '2': 118, '4': 32, '6': 30, '9': 9}
    
#     print(ribo_num)
    if ribo_num > 1:

#         poly_frac = 0

#         for poly_size in range(2,11):

#             try:

#                 poly_frac = poly_frac + ribo_num*poly_size_dict[str(poly_size)]/684  #min(ribo_num,poly_size)

#             except:
        
#                 poly_frac = poly_frac + ribo_num*poly_size/684
        
#         kcat_mod = poly_frac*riboKcat + 0.05*riboKcat #*

        kcat_mod = ribo_num*0.25*riboKcat + 0.25*riboKcat # 503 ribosomes
    
#         kcat_mod = ribo_num*0.40*riboKcat + 0.05*riboKcat # 684 ribosomes
        
    else:
        
        kcat_mod = 0.45*riboKcat
    
    k_translation = kcat_mod / ((1+riboK0/ribosomeConc)*(riboKd**2)/(ctRNAconc**2) + NMonoSum + n_tot - 1)
    
    return k_translation


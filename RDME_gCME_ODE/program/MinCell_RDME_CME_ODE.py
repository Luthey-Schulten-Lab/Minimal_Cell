from jLM.RegionBuilder import RegionBuilder
from jLM.RDME import Sim as RDMESim
from jLM.RDME import File as RDMEFile
import jLM

from jLM.Solvers import makeSolver

from pyLM import CME

from pyLM.units import *

import lm

from lm import MpdRdmeSolver
from lm import IntMpdRdmeSolver

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy  as np

import os

import scipy.ndimage as spnd
import ipywidgets as ipw
import h5py
import itertools
import random
import copy

# import ipyvolume
# from sidecar import Sidecar
# import numpy as np
# from ipywebrtc import WidgetStream, VideoRecorder

import math
import scipy as sp
import scipy.spatial

# import seaborn as sns

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import time

try:
        from tqdm import tqdm
#         print('Imported tqdm')
except:
        def tqdm(x,ascii=False):
                return x


%matplotlib notebook

import warnings
warnings.simplefilter("ignore")

import sys, getopt

delt = 1.0 #s
odestep = 0.1 # s
cythonBool = False
totalTime = 10.0 #s
gpuID = sys.argv[0]
repID = sys.argv[1]

simDir = './sims/cell_' + str(repID)

try:
    os.mkdir('sims')
    print('Created sims directory')
except:
    print('Sims directory exists')
    
try:
    os.mkdir(simDir)
    print('Created sim directory for cell ' + str(repID))
except:
    print('Sim directory already exists')

def initSim():
#     filename='./MinCell_jLM_RDME_CME_ODE_20min.lm'
    filename= simDir + '/MinCell_jLM_RDME_CME_ODE.lm'

    N_edges = 64 # Number of subvolumes making up and edge of the simulation space N x N x N

    N_2 = N_edges/2

    sim = RDMESim("JCVI-syn3A",
                  filename,
                  [N_edges,N_edges,N_edges],
                  8e-9,
                  "extracellular")

    cyto_radius = 2.00e-7/sim.latticeSpacing #m converted to lattice sites (8 nm lattice spacing)
    dna_monomers = 46188

    cyto_vol = (4/3)*np.pi*0.200**3

    cyto_200 = (4/3)*np.pi*0.2**3

    ptn_ratio = (2.3e6*cyto_vol)/(2.3e6*cyto_200)
#     print(ptn_ratio)

#     riboFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/08262020/s1c8_free_origin_rep00001_obst_CG.dat'
    # riboFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/08292020/s1c8_fixed_origin_rep00001_obst_CG.dat'
    riboFile = '../model_data/s1c15/s1c15_coords_nm_adaptive_fitting_s1c15_trans_id_8nm.txt'

#     riboFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/08292020/s1c8_fixed_origin_rep00001_obst_CG.dat' 
#     riboFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/09012020/s1c14_coords_trans_id.txt'

    dnaFile = '../model_data/s1c15/s1c15_base_CG_reps00001_00090/s1c15_base/CG/s1c15_base_rep00057_CG_coords.dat'
    dnaPartFile = '../model_data/s1c15/s1c15_base_CG_reps00001_00090/s1c15_base/CG/s1c15_base_rep00057_FG_nodes.dat'

#     dataFolder = '/home/zane/MinCell/LM_sims/model_data/s1c14/'

#     riboFile = dataFolder + 's1c14_base_rep00001_obst_CG.dat'

#     dnaFile = dataFolder + 's1c14_base_rep00001_CG_coords.dat'

#     dnaPartFile = dataFolder + 's1c14_base_rep00001_FG_nodes.dat'

#     dnaFile = dataFolder + 's1c14_2chrom00001_CG/CG/s1c14_2chrom00001_rep00005_CG_coords.dat'

#     dnaPartFile = dataFolder + 's1c14_2chrom00001_CG/CG/s1c14_2chrom00001_rep00005_FG_nodes.dat'
    
    sim_center = [N_2,N_2,N_2]

    sim.timestep = 30e-6
    sim.simulationTime=totalTime
    sim.latticeWriteInterval=1.0
    sim.speciesWriteInterval=1.0
    replicates = 1
    
#     sim.hookInterval(delt)
    
    pmap = {}
    
    PartIdxMap = {}
    
#     print('Configuration ' + str(rep+1) + '/' + str(len(configs)) + ' initialized')
    print('Simulation Initialized')
    
    return sim, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, riboFile, dnaFile, dnaPartFile, filename, PartIdxMap, pmap

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

gene_list = []
for i in range(len(dna3A.features)):
    if ('product' in dna3A.features[i].qualifiers.keys()):
        #print(i) # This first statement works
        #print(dna.features[i].qualifiers['product'])
        if dna3A.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
            #print(dna.features[i].qualifiers['product'])
            gene_list.append(i)
# gene_list

gene_starts = []

for gene in gene_list:
    
    locusTag = dna3A.features[gene].qualifiers['locus_tag'][0]
    gene_start = dna3A.features[gene].location.start.real
    
    direction = dna3A.features[gene].strand
    
    gene_starts.append([locusTag,gene_start,direction])
    
# gene_starts


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
        
#         ptnConcentration = ptnCount*countToMiliMol
    except:
#         print("WARNING: No protein count for", newMetID)
#         print("Using default protein concentration.")

        ptnName = newMetID
        ptnCount = defaultPtnCount
#         ptnConcentration = defaultPtnConcentration

        ptnNoQuant.add(newMetID)
    
    return ptnCount, ptnName


PtnMetDF = pd.read_csv("../model_data/protein_metabolites_frac.csv")
# PtnMetDF

riboPtnMetDF = pd.read_csv("../model_data/ribo_protein_metabolites.csv")
# riboPtnMetDF

memPtnMetDF = pd.read_csv("../model_data/membrane_protein_metabolites.csv")
# memPtnMetDF

MetPtnGenes = ['MMSYN1_0445', 'MMSYN1_0220',
 'MMSYN1_0131', 'MMSYN1_0727', 'MMSYN1_0607',
 'MMSYN1_0451', 'MMSYN1_0606', 'MMSYN1_0729',
 'MMSYN1_0213', 'MMSYN1_0221', 'MMSYN1_0495',
 'MMSYN1_0494', 'MMSYN1_0493', 'MMSYN1_0726',
 'MMSYN1_0435', 'MMSYN1_0475', 'MMSYN1_0227',
 'MMSYN1_0228', 'MMSYN1_0229', 'MMSYN1_0230',
 'MMSYN1_0316', 'MMSYN1_0262', 'MMSYN1_0800',
 'MMSYN1_0831', 'MMSYN1_0733', 'MMSYN1_0732',
 'MMSYN1_0432', 'MMSYN1_0381', 'MMSYN1_0687',
 'MMSYN1_0688', 'MMSYN1_0689', 'MMSYN1_0012',
 'MMSYN1_0519', 'MMSYN1_0260', 'MMSYN1_0634',
 'MMSYN1_0837', 'MMSYN1_0126', 'MMSYN1_0535',
 'MMSYN1_0613', 'MMSYN1_0308', 'MMSYN1_0061',
 'MMSYN1_0222', 'MMSYN1_0282', 'MMSYN1_0287',
 'MMSYN1_0076', 'MMSYN1_0064', 'MMSYN1_0288',
 'MMSYN1_0528', 'MMSYN1_0529', 'MMSYN1_0163',
 'MMSYN1_0405', 'MMSYN1_0614', 'MMSYN1_0380',
 'MMSYN1_0378', 'MMSYN1_0259', 'MMSYN1_0291',
 'MMSYN1_0823', 'MMSYN1_0684', 'MMSYN1_0390',
 'MMSYN1_0799', 'MMSYN1_0443', 'MMSYN1_0413',
 'MMSYN1_0747', 'MMSYN1_0651', 'MMSYN1_0771',
 'MMSYN1_0772', 'MMSYN1_0773', 'MMSYN1_0819',
 'MMSYN1_0344', 'MMSYN1_0330', 'MMSYN1_0382',
 'MMSYN1_0216', 'MMSYN1_0203', 'MMSYN1_0798',
 'MMSYN1_0537', 'MMSYN1_0347', 'MMSYN1_0140',
 'MMSYN1_0045', 'MMSYN1_0129', 'MMSYN1_0515',
 'MMSYN1_0447', 'MMSYN1_0218', 'MMSYN1_0513',
 'MMSYN1_0139', 'MMSYN1_0420', 'MMSYN1_0616',
 'MMSYN1_0617', 'MMSYN1_0419', 'MMSYN1_0117',
 'MMSYN1_0512', 'MMSYN1_0304', 'MMSYN1_0875',
 'MMSYN1_0214', 'MMSYN1_0147', 'MMSYN1_0115',
 'MMSYN1_0813', 'MMSYN1_0814', 'MMSYN1_0114',
 'MMSYN1_0697', 'MMSYN1_0113']
# MetPtnGenes

MetLocusNums = []

for gene in MetPtnGenes:
    locusNum = gene.split('_')[1]
    MetLocusNums.append(locusNum)
    
# MetLocusNums

rrnaMetDF_1 = pd.read_csv("../model_data/rrna_metabolites_1.csv")
# rrnaMetDF_1

rrnaMetDF_2 = pd.read_csv("../model_data/rrna_metabolites_2.csv")
# rrnaMetDF_2

trnaMetDF = pd.read_csv("../model_data/trna_metabolites_synthase.csv")
# trnaMetDF

named_PTN_list = []

for index, row in riboPtnMetDF.iterrows():
    named_PTN_list.append(row["gene"]) 

    
for index, row in PtnMetDF.iterrows():
    print(row["gene"])
    named_PTN_list.append(row["gene"])

baseMap = OrderedDict({ "A":"M_atp_c", "U":"M_utp_c", "G":"M_gtp_c", "C":"M_ctp_c" })
# baseMapToMonoP = OrderedDict({ "A":"M_amp_c", "U":"M_ump_c", "G":"M_gmp_c", "C":"M_cmp_c" })

# Global parameters for transcription
rnaPolKcat = 20 # nt/s
rnaPolK0 = 1e-4 #mM
rnaPolKd = 0.01 #mM

rrnaPolKcat = 85 # nt/s

krnadeg = 0.00578/2 # 1/s
# rna_deg_rate = sim.rateConst('RNAdeg', krnadeg, 2)

ptnDegRate = 7.70e-06 # 1/s

ATPconc = 1.04 #mM
UTPconc = 0.68 #mM
CTPconc = 0.34 #mM
GTPconc = 0.68 #mM

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


defaultPtnCount = 1

# Global parameters for translation
riboKcat = 12 # 1/s
riboK0 = 4*25e-6 # mM
riboKd = 0.0001 # mM

ribo_init = 5*Ecoli_V*avgdr/60/6800

ribosomeConc = 500*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 150*countToMiliMol # mM

# Global parameter for degradation of proteins
# Derived from eLife's model, using average protein half life of 25 hours. 
ptnDegRate = 7.70e-06 # 1/s

trnaCmeMetDF = pd.read_csv("../model_data/trna_metabolites_synthase.csv")


from diffusion import *
# from MC_CME import *
from MC_RDME import *
from regions_and_complexes import *
from GIP_rates import *
# import hook as hook
# import setICs as setICs

rep = 1
partIdx = 1

sim, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, riboFile, dnaFile, dnaPartFile, filename, PartIdxMap, pmap = initSim()
    
sim, genePoints, ribo_points, ribo_center_points, ext, mem, cyt, ribo, dna, she, cyto_shell, partIdx = buildRegions(sim, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, riboFile, dnaFile, filename, pmap, PartIdxMap, partIdx)

# sim.finalize()

sim, geneEnds, geneStarts, singleStatePtnDict, multiStatePtnDict, degDict, tRNAstateDict, RDME_species_list, partIdx, rtRNA_ID_dict = constructRDME(sim, pmap, genePoints, ribo_points, ribo_center_points, ext, mem, cyt, ribo, dna, she, cyto_shell, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, dnaPartFile, gene_starts, PtnMetDF, riboPtnMetDF, memPtnMetDF, trnaMetDF, genomePtnLocDict, PartIdxMap, partIdx)

import setICs
setICs.__main__(pmap)

import hook as hook
rdmeCmeOdeHookSolver = hook.MyOwnSolver

Solver = makeSolver(IntMpdRdmeSolver, rdmeCmeOdeHookSolver)
solver = Solver(sim, delt, odestep, cythonBool, pmap, totalTime, geneEnds, geneStarts, singleStatePtnDict, multiStatePtnDict, degDict, tRNAstateDict, RDME_species_list, PartIdxMap, rtRNA_ID_dict)

sim.finalize()

sim.run(solver=solver, cudaDevices=[gpuID])

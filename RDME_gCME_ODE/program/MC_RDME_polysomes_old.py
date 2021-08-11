###### Genetic Information Processes Reactions RDME  #######

##################
# Define how to add the particles and reactions for each protein and its corresponding mRNA and gene.

# mrna_diff_list = []

from regions_and_complexes_polysomes import *
from diffusion import *
from GIP_rates import *

from jLM.RegionBuilder import RegionBuilder
import random

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict


def constructRDME(sim, pmap, genePoints, ribo_points, ribo_center_points, ext, mem, cyt, ribo, dna, she, cyto_shell, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, dnaPartFile, gene_starts, PtnMetDF, riboPtnMetDF, memPtnMetDF, trnaMetDF, genomePtnLocDict, PartIdxMap, partIdx):
    
    defaultDiffs(sim, ext, mem, cyt, ribo, dna, she)
    
    mRNA_counts_DF = pd.read_csv('../model_data/mRNA_Counts.csv')
    rrnaMetDF_1 = pd.read_csv("../model_data/rrna_metabolites_1.csv")
    rrnaMetDF_2 = pd.read_csv("../model_data/rrna_metabolites_2.csv")
    
    degradosome = sim.species('Degradosome')
    
    RDME_species_list = ['Degradosome','RNApol','Ribosome1','secY']
    
    singleStatePtnDict = {}
    multiStatePtnDict = {}
    degDict = {}
    tRNAstateDict = {}
    rtRNA_ID_dict = {}
    
    aaCostMap = OrderedDict({"A":"ALA_cost","R":"ARG_cost","N":"ASN_cost","D":"ASP_cost","C":"CYS_cost",
                         "E":"GLU_cost","Q":"GLN_cost","G":"GLY_cost","H":"HIS_cost","I":"ILE_cost",
                         "L":"LEU_cost","K":"LYS_cost","M":"MET_cost","F":"PHE_cost","P":"PRO_cost",
                         "S":"SER_cost","T":"THR_cost","W":"TRP_cost","Y":"TYR_cost","V":"VAL_cost",
                         "FM":"FMET_cost"})
    
    geneOccupancies = readDNAoccupancies(dnaPartFile)
    DNA_map = mapDNA(gene_starts, dna_monomers)
    sim, geneEnds, geneStarts, partIdx = addDNApart(sim, DNA_map, genePoints, geneOccupancies, ext, mem, cyt, ribo, dna, she, RDME_species_list, PartIdxMap, partIdx)
    
    RNApol = sim.species('RNApol')
    PartIdxMap['RNApol'] = partIdx
    partIdx = partIdx + 1
    Ribosome = sim.species('Ribosome1')
    PartIdxMap['Ribosome1'] = partIdx
    partIdx = partIdx + 1
    
    addRNAPpart(sim, pmap, ptn_ratio, RNApol, ext, mem, cyt, ribo, dna, she)
    
    sim, poly_ribo, partIdx = addRiboPart(sim, pmap, partIdx, PartIdxMap, ext, mem, cyt, ribo, dna, she, ribo_center_points, RDME_species_list)
    
    poly_sizes = [len(poly) for poly in poly_ribo]
    Max_Poly_Count = max(poly_sizes)
    print('Maximum polysome size = ', str(Max_Poly_Count))
    
    secY_init = 1.0e6 #5*Ecoli_V*avgdr/60/6800

    secY_on = sim.rateConst('secY_on', secY_init, 2)
    secY_off = sim.rateConst('secY_off', 1e-3, 1)
    
    genes_in_model_RDME = []
    ptnCounts = []
    
#     diffRNA = sim.diffusionConst("mrna_diff", 0.3e-12)
    diffPTN = sim.diffusionConst("ptn_diff", 1e-12)
    diffPTNslow = sim.diffusionConst("ptn_diff", 0.001e-12)
    
#     RNAP_on = sim.rateConst('RNAP_on', 1e6, 2)
    RNAP_off = sim.rateConst('RNAP_off', 1e-3, 1)
    Ribo_on = sim.rateConst('Ribo_on', ribo_init, 2)
    Ribo_off = sim.rateConst('Ribo_off', 1e-3, 1)
    deg_bind_rate = sim.rateConst('RNAdeg', degrad_bind_rate, 2)
    
    pmap['Trsc_Term'] = 0
    
    RDME_species_list.append('Trsc_Term')
    Trsc_Term = sim.species('Trsc_Term')
    PartIdxMap['Trsc_Term'] = partIdx
    partIdx = partIdx + 1
    
    TrscTermination = sim.rateConst('TrscTermination',1e9,2)
    
#     added = 0
#     for ID in genomePtnLocDict:
#         if ID not in genes_in_model:
#             added+=1
#             addPtn(ID, sim, ptn_ratio)
          
            
    for index, row in PtnMetDF.iterrows():
        partIdx = addNamedPtnRDME(sim, pmap, row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"], ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, singleStatePtnDict, multiStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx)

    for index, row in riboPtnMetDF.iterrows():
        partIdx = addRiboPtnRDME(sim, pmap, row["species"], row["gene"], ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, singleStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx)
        
    for index, row in memPtnMetDF.iterrows():
        partIdx = addMembranePtnRDME(sim, pmap, row["species"], row["gene"], row["transcribe"], row["proteomics_fraction"], ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, secY_on, secY_off, degDict, multiStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx)
        
    binding_rate = 4*(180/765)*Ecoli_V*avgdr/11400/60 #/1800/60
    RNAP_on = sim.rateConst('RNAP_on', binding_rate, 2)  
    tRNAadded = []
    for index, row in trnaMetDF.iterrows():
#         addtRNA(row["species"], row["gene"], sim, tRNAadded, RNAP_on)
    
        partIdx = addtRNA(sim, pmap, row["uncharged"], row["charged"], row["gene"], row["AA"], row["synthase"], ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_on, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, multiStatePtnDict, aaCostMap, tRNAstateDict, tRNAadded, RDME_species_list, PartIdxMap, partIdx, rtRNA_ID_dict)
        
    print('Added tRNA')

    for jcvi3AID in genomePtnLocDict:
        if jcvi3AID not in genes_in_model_RDME:
    #         added+=1
            partIdx = addPtnRDME(jcvi3AID, sim, ptn_ratio, pmap, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, ext, mem, cyt, ribo, dna, she, degDict, singleStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx)
        
    print('Added all proteins')
    print(partIdx)
        
    binding_rate = 4*Ecoli_V*avgdr/11400/60 #/1800/60
    RNAP_on = sim.rateConst('RNAP_on_rRNA', binding_rate, 2)

    addrRNA(sim,pmap,ext, mem, cyt, ribo, dna, she,rrnaMetDF_1,rrnaMetDF_2,RNApol,RNAP_on,RNAP_off,RDME_species_list)
    print('Added rRNA')
        
    
    pmap['P_227'] = 182
    pmap['P_227'] = pmap['M_lpl_PdhC_c']+pmap['M_dhlpl_PdhC_c']+pmap['M_acdhlpl_PdhC_c']
    
    multiStatePtnDict['ptsg']['Type'] = 'ptsg'
    multiStatePtnDict['ptsg']['States'] = ['ptsg','ptsg_P','S_779','ptsg_C']
    
    multiStatePtnDict['ptsi']['States'] = ['ptsi','ptsi_P']
    multiStatePtnDict['ptsh']['States'] = ['ptsh','ptsh_P']
    multiStatePtnDict['crr']['States'] = ['crr','crr_P']
    
    multiStatePtnDict['M_trdox_c']['States'] = ['M_trdox_c','M_trdrd_c']
    multiStatePtnDict['M_lpl_PdhC_c']['States'] = ['M_lpl_PdhC_c','M_dhlpl_PdhC_c','M_acdhlpl_PdhC_c']
    
    multiStatePtnDict['M_apoACP_c']['States'] = ['M_apoACP_c','M_ACP_c','M_ACP_R_c']
    
    synthase_ptn = 'P_126'
    synthase_atp = synthase_ptn + '_ATP'
    synthase_atp_aa = synthase_atp + '_AA'
    
    print(multiStatePtnDict)
    
    multiStatePtnDict['P_126']['States'] = ['P_126','P_126_ATP','P_126_ATP_AA','P_126_ATP_AA_tRNAgln','P_126_ATP_AA_tRNA']
    
    tRNAstateDict['M_trnamet_c'].append('M_fmettrna_c')

#     binding_rate = 4*Ecoli_V*avgdr/205/60
#     RNAP_on = sim.rateConst('RNAP_on', binding_rate, 2)
#     addrRNA(sim)

    Nuc_counters = ['ATP_trsc','ATP_translat','ATP_mRNAdeg','ATP_ptndeg','ATP_DNArep','ATP_transloc',
               'ATP_mRNA','UTP_mRNA','CTP_mRNA','GTP_mRNA',
               'AMP_mRNAdeg','UMP_mRNAdeg','CMP_mRNAdeg','GMP_mRNAdeg',
               'ATP_tRNA','UTP_tRNA','CTP_tRNA','GTP_tRNA',
               'ATP_rRNA','UTP_rRNA','CTP_rRNA','GTP_rRNA',
               'dATP_DNArep','dTTP_DNArep','dCTP_DNArep','dGTP_DNArep']

#     sim.defineSpecies(Nuc_counters)
    for cost in Nuc_counters:
        pmap[cost] = 0

    AA_counters = ["ALA_cost","ARG_cost","ASN_cost","ASP_cost","CYS_cost","GLU_cost","GLN_cost","GLY_cost",
                   "HIS_cost","ILE_cost","LEU_cost","LYS_cost","MET_cost","PHE_cost","PRO_cost","SER_cost",
                   "THR_cost","TRP_cost","TYR_cost","VAL_cost","FMET_cost"]

#     sim.defineSpecies(AA_counters)
    
    for cost in AA_counters:
        
        pmap[cost] = 0

    AA_paid = []

    for cost in AA_counters:

        cost_paid = cost + '_paid'
        pmap[cost_paid] = 0
#         AA_paid.append(cost_paid)

    
    return sim, geneEnds, geneStarts, singleStatePtnDict, multiStatePtnDict, degDict, tRNAstateDict, RDME_species_list, partIdx, rtRNA_ID_dict
    
    

def addPtnRDME(jcvi3AID, sim, ptn_ratio, pmap, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, ext, mem, cyt, ribo, dna, she, degDict, singleStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx):
    locusNum = jcvi3AID.split('_')[1]
    mmcode = 'MMSYN1_' + locusNum
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode
        
    locusNum = jcvi3AID.split('_')[1].lstrip('0')

#     print(mmcode, jcvi2ID, jcvi3AID)
    
    genes_in_model_RDME.append(jcvi3AID)
    
    # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
    ptnMetID = 'P_' + locusNum
    RDME_species_list.append(ptnMetID)
    
    # If the protein is not in the model, add it:

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    ptnCount = max(2,ptnCount)
    
#     ptnCounts.append(ptnCount)

#     print(ptnMetID, ptnCount)

    geneMetID = 'g_' + locusNum

    # Get nucleotide and amino acid sequences, if available
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):
        
        aaCount = defaultdict(int)
        for aa in set(aasequence):
            aaCount[aa] = aasequence.count(aa)
        ptn_len = sum(list(aaCount.values()))
        
        if ptn_len <= 134:
            max_poly_size = 1
        elif ptn_len > 134:
            max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

        rnaMetID = "R_" + locusNum
        rnaName = "(mRNA) " + ptnName
        RDME_species_list.append(rnaMetID)
        
        pmap['New_mRNA_' + locusNum] = 0
    
        species = []
        species = [geneMetID, rnaMetID, ptnMetID]
        gene = sim.species(geneMetID)
        
        RNAP_gene = sim.species('RP_' + locusNum)
        PartIdxMap['RP_' + locusNum] = partIdx
        partIdx = partIdx + 1
        pmap['RP_' + locusNum] = 0
        RDME_species_list.append('RP_' + locusNum)
        
        RNA = sim.species(rnaMetID)
        PartIdxMap[rnaMetID] = partIdx
        partIdx = partIdx + 1
        
        Deg_RNA = sim.species('D_' + locusNum)
        PartIdxMap['D_' + locusNum] = partIdx
        partIdx = partIdx + 1
        pmap['D_' + locusNum] = 0
        RDME_species_list.append('D_' + locusNum)
        
        PTN = sim.species(ptnMetID)
        PartIdxMap[ptnMetID] = partIdx
        partIdx = partIdx + 1
        
        for mp in range(Max_Poly_Count):
            
            for pn in range(mp+1): #range(max_poly_size):
                
                if (max_poly_size == 1) or (mp == 0):
                    
                    Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                    Ribo_RNA = sim.species(Ribo_ID)
                    PartIdxMap[Ribo_ID] = partIdx
                    partIdx = partIdx + 1
                    pmap[Ribo_ID] = 0
                    RDME_species_list.append(Ribo_ID)
            
                    break
                
                elif max_poly_size > 1:
                    
                    if pn == 0:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)
            
                    elif (pn+1 == max_poly_size) or (pn == mp):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)
                        
                        break
                        
                    else:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)
        
        for index,row in mRNA_counts_DF.iterrows():
        
            if row["LocusTag"] == jcvi3AID:

                avg_mRNA_cnt = row["Count"]

                if avg_mRNA_cnt == 0.0:
                    avg_mRNA_cnt = 0.001

                init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
                continue
        
        sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
        pmap[rnaMetID] = init_mRNA_count
        sim.distributeNumber(PTN, sim.region('cytoplasm'), int(4*ptnCount*ptn_ratio/5))
        sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))
        pmap[ptnMetID] = int(ptnCount*ptn_ratio)
        
        diffRNA = sim.diffusionConst('mrna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))
#         mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

        RNA.diffusionRate(cyt, diffRNA)
        PTN.diffusionRate(cyt, diffPTN)
        
        RNA.diffusionRate(mem,sim.diffusionZero)
        RNA.diffusionRate(ext,sim.diffusionZero)
        RNA.diffusionRate(ribo,diffRNA)
        RNA.diffusionRate(dna,diffRNA)
        RNA.diffusionRate(she,diffRNA)

        PTN.diffusionRate(mem,sim.diffusionZero)
        PTN.diffusionRate(ext,sim.diffusionZero)
        PTN.diffusionRate(ribo,diffPTN)
        PTN.diffusionRate(dna,diffPTN)
        PTN.diffusionRate(she,diffPTN)
        
        sim.transitionRate(RNA, dna, cyt, diffRNA)
        sim.transitionRate(RNA, cyt, dna, diffRNA)
        sim.transitionRate(RNA, cyt, ribo, diffRNA)
        sim.transitionRate(RNA, ribo, cyt, diffRNA)
        sim.transitionRate(RNA, ribo, she, diffRNA)
        sim.transitionRate(RNA, she, ribo, diffRNA)
        sim.transitionRate(RNA, she, cyt, diffRNA)
        sim.transitionRate(RNA, cyt, she, diffRNA)
        sim.transitionRate(RNA, she, dna, diffRNA)
        sim.transitionRate(RNA, dna, she, diffRNA)
        sim.transitionRate(RNA, dna, ribo, diffRNA)
        sim.transitionRate(RNA, ribo, dna, diffRNA)
        
        sim.transitionRate(PTN, dna, cyt, diffPTN)
        sim.transitionRate(PTN, cyt, dna, diffPTN)
#         sim.transitionRate(PTN, cyt, ribo, diffPTN)
        sim.transitionRate(PTN, ribo, cyt, diffPTN)
        sim.transitionRate(PTN, ribo, she, diffPTN)
#         sim.transitionRate(PTN, she, ribo, diffPTN)
        sim.transitionRate(PTN, she, cyt, diffPTN)
        sim.transitionRate(PTN, cyt, she, diffPTN)
        sim.transitionRate(PTN, she, dna, diffPTN)
        sim.transitionRate(PTN, dna, she, diffPTN)
        sim.transitionRate(PTN, dna, ribo, diffPTN)
        sim.transitionRate(PTN, ribo, dna, diffPTN)
        
#         trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)
        
        translat_rate = sim.rateConst(locusNum + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
        
        if max_poly_size > 1:
            
            aasequence_poly = aasequence[-134:]
        
            translat_rate_poly = sim.rateConst(locusNum + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)
        
        rna_deg_rate = sim.rateConst(locusNum + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)
        
        RNAP_on = sim.rateConst(locusNum + 'RNAP_on', RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)
        
        dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)
        
        dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)
        
#         dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)
        
#         dna.addReaction([gene], [gene, RNA], trsc_rate)
        
        she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)
        
        she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)
        
#         ribo.addReaction([RNA], [], rna_deg_rate)

        for mp in range(Max_Poly_Count):
        
            Ribosome = sim.species('Ribosome'+str(mp+1))
            
            for pn in range(mp+1): #range(max_poly_size):
                
                if (max_poly_size == 1) or (mp == 0):
                    
                    Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                    Ribo_RNA = sim.species(Ribo_ID)
                    
                    ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                    ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)
                    
                    break
                
                elif max_poly_size > 1:
                    
                    if pn == 0:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)
                        
                        ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)
        
                        Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)
                
                        Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                        ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate)
            
                    elif (pn+1 == max_poly_size) or (pn+1 == mp+1):
                    
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)

                        ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate_poly)
                        
                        break
                        
                    else:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)
                        
                        Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)
                
                        Ribo_RNA_Next = sim.species(Ribo_ID_Next)
                        
                        ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate_poly)
        

#         Trsc_Term = sim.species('Trsc_Term')
        
#         dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)
        

        baseCount = defaultdict(int)
        for base in set(rnasequence):
            baseCount[base] = rnasequence.count(base)

        # Add total number of monomers to parameter dict

        N_A = baseCount["A"]

        N_U = baseCount["U"]

        N_C = baseCount["C"]

        N_G = baseCount["G"]
        
        gene_deg_dict = {}
        
        gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G
        
        gene_deg_dict['AMP_mRNAdeg'] = N_A
        gene_deg_dict['UMP_mRNAdeg'] = N_U
        gene_deg_dict['CMP_mRNAdeg'] = N_C
        gene_deg_dict['GMP_mRNAdeg'] = N_G
        
        degDict['D_' + locusNum] = gene_deg_dict
        
        
        gene_translat_dict = {}
        
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
        
        translatATP = 0
        
        for aaCost in AaUsed:
            
            aa_ID = aaCost[0]
            numberUsed = aaCost[1]

            aaCostID = aaCostMap[aa_ID]
            
            gene_translat_dict[aaCostID] = numberUsed
            
            translatATP = translatATP + 2*numberUsed
            
            
        gene_translat_dict['ATP_translat'] = translatATP
        
        singleStatePtnDict[ptnMetID] = gene_translat_dict    


        return partIdx
##################



def addMembranePtnRDME(sim, pmap, newMetID, mmcode, transcribe, ptnFrac, ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, secY_on, secY_off, degDict, multiStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx):
    
    locusNum = mmcode.split('_')[1]
    
#     locusNum = jcvi3AID.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum

#     mmcode = 'MMSYN1_' + locusNum
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode

    print(mmcode, jcvi2ID, jcvi3AID)
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
    if 'ptsg' in newMetID:
        ptnMetID = newMetID
    else:
        ptnMetID = 'P_' + locusNum
#     
    
    # If the protein is not in the model, add it:

#     ModelSpecies.append(ptnMetID)

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    ptnCount = max(2,ptnCount)
    
#     ptnCount = ptnCount*ptnFrac
    
#     ptnCounts.append(ptnCount)

    print(ptnMetID, ptnCount)
    
    if transcribe:
        
        RDME_species_list.append(ptnMetID)
        
        genes_in_model_RDME.append(jcvi3AID)

        geneMetID = 'g_' + locusNum

#         ModelSpecies.append(geneMetID)

        # Get nucleotide and amino acid sequences, if available
        rnasequence, aasequence = getSequences(jcvi3AID)
    #     print(rnasequence)
    #     print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):
            
            aaCount = defaultdict(int)
            for aa in set(aasequence):
                aaCount[aa] = aasequence.count(aa)
            ptn_len = sum(list(aaCount.values()))

            if ptn_len <= 134:
                max_poly_size = 1
            elif ptn_len > 134:
                max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

            rnaMetID = "R_" + locusNum
            rnaName = "(mRNA) " + ptnName
            RDME_species_list.append(rnaMetID)
            
            pmap['New_mRNA_' + locusNum] = 0

#             ModelSpecies.append(rnaMetID)

#             species = []
#             species = [geneMetID, rnaMetID, ptnMetID]
#             gene = sim.species(geneMetID)
#             RNAP_gene = sim.species('RNAP_' + locusNum)
#             RNA = sim.species(rnaMetID)
#             Ribo_RNA = sim.species('Ribo_' + locusNum)
#             Deg_RNA = sim.species('Deg_' + locusNum)
#             PTN = sim.species(ptnMetID)
            
            species = []
            species = [geneMetID, rnaMetID, ptnMetID]
            gene = sim.species(geneMetID)

            RNAP_gene = sim.species('RP_' + locusNum)
            PartIdxMap['RP_' + locusNum] = partIdx
            partIdx = partIdx + 1
            pmap['RP_' + locusNum] = 0
            RDME_species_list.append('RP_' + locusNum)

            RNA = sim.species(rnaMetID)
            PartIdxMap[rnaMetID] = partIdx
            partIdx = partIdx + 1

            Deg_RNA = sim.species('D_' + locusNum)
            PartIdxMap['D_' + locusNum] = partIdx
            partIdx = partIdx + 1
            pmap['D_' + locusNum] = 0
            RDME_species_list.append('D_' + locusNum)

            PTN = sim.species(ptnMetID)
            PartIdxMap[ptnMetID] = partIdx
            partIdx = partIdx + 1

            for mp in range(Max_Poly_Count):

                for pn in range(mp+1): #range(max_poly_size):

                    if (max_poly_size == 1) or (mp == 0):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)

                        break

                    elif max_poly_size > 1:

                        if pn == 0:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                        elif (pn+1 == max_poly_size) or (pn == mp):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                            break

                        else:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)
            
            cytoPtnMetID = ptnMetID + '_C'
            cPTN = sim.species(cytoPtnMetID)
            PartIdxMap[cytoPtnMetID] = partIdx
            partIdx = partIdx + 1
            pmap[cytoPtnMetID] = 0
            RDME_species_list.append(cytoPtnMetID)

#             ModelSpecies.append('Ribo_' + locusNum)

            for index,row in mRNA_counts_DF.iterrows():
        
                if row["LocusTag"] == jcvi3AID:

                    avg_mRNA_cnt = row["Count"]

                    if avg_mRNA_cnt == 0.0:
                        avg_mRNA_cnt = 0.001

                    init_mRNA_count = np.random.poisson(avg_mRNA_cnt)

                    continue
        
            sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
            pmap[rnaMetID] = init_mRNA_count
            sim.distributeNumber(PTN, sim.region('membrane'), int(ptnCount*ptn_ratio))
            ptnCount = ptnCount*ptnFrac*ptn_ratio
            pmap[ptnMetID] = int(ptnCount)
#             sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))

            diffRNA = sim.diffusionConst('mrna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))
#             mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

            RNA.diffusionRate(cyt, diffRNA)
            cPTN.diffusionRate(cyt, diffPTN)

            RNA.diffusionRate(mem,sim.diffusionZero)
            RNA.diffusionRate(ext,sim.diffusionZero)
            RNA.diffusionRate(ribo,diffRNA)
            RNA.diffusionRate(dna,diffRNA)
            RNA.diffusionRate(she,diffRNA)

            PTN.diffusionRate(mem,sim.diffusionZero)
            PTN.diffusionRate(ext,sim.diffusionZero)
            PTN.diffusionRate(ribo,sim.diffusionZero)
            PTN.diffusionRate(dna,sim.diffusionZero)
            PTN.diffusionRate(she,diffPTN)
            PTN.diffusionRate(mem,diffPTNslow)
            
            cPTN.diffusionRate(mem,sim.diffusionZero)
            cPTN.diffusionRate(ext,sim.diffusionZero)
            cPTN.diffusionRate(ribo,diffPTN)
            cPTN.diffusionRate(dna,diffPTN)
            cPTN.diffusionRate(she,diffPTN)

            sim.transitionRate(RNA, dna, cyt, diffRNA)
            sim.transitionRate(RNA, cyt, dna, diffRNA)
            sim.transitionRate(RNA, cyt, ribo, diffRNA)
            sim.transitionRate(RNA, ribo, cyt, diffRNA)
            sim.transitionRate(RNA, ribo, she, diffRNA)
            sim.transitionRate(RNA, she, ribo, diffRNA)
            sim.transitionRate(RNA, she, cyt, diffRNA)
            sim.transitionRate(RNA, cyt, she, diffRNA)
            sim.transitionRate(RNA, she, dna, diffRNA)
            sim.transitionRate(RNA, dna, she, diffRNA)
            sim.transitionRate(RNA, dna, ribo, diffRNA)
            sim.transitionRate(RNA, ribo, dna, diffRNA)

#             sim.transitionRate(PTN, dna, cyt, diffPTN)
#             sim.transitionRate(PTN, cyt, dna, diffPTN)
#     #         sim.transitionRate(PTN, cyt, ribo, diffPTN)
#             sim.transitionRate(PTN, ribo, cyt, diffPTN)
#             sim.transitionRate(PTN, ribo, she, diffPTN)
#     #         sim.transitionRate(PTN, she, ribo, diffPTN)
#             sim.transitionRate(PTN, she, cyt, diffPTN)
#             sim.transitionRate(PTN, cyt, she, diffPTN)
#             sim.transitionRate(PTN, she, dna, diffPTN)
#             sim.transitionRate(PTN, dna, she, diffPTN)
#             sim.transitionRate(PTN, dna, ribo, diffPTN)
#             sim.transitionRate(PTN, ribo, dna, diffPTN)
            sim.transitionRate(PTN, she, mem, diffPTN)
#             sim.transitionRate(PTN, mem, she, diffPTN)

            sim.transitionRate(cPTN, dna, cyt, diffPTN)
            sim.transitionRate(cPTN, cyt, dna, diffPTN)
    #         sim.transitionRate(cPTN, cyt, ribo, diffPTN)
            sim.transitionRate(cPTN, ribo, cyt, diffPTN)
            sim.transitionRate(cPTN, ribo, she, diffPTN)
    #         sim.transitionRate(cPTN, she, ribo, diffPTN)
            sim.transitionRate(cPTN, she, cyt, diffPTN)
            sim.transitionRate(cPTN, cyt, she, diffPTN)
            sim.transitionRate(cPTN, she, dna, diffPTN)
            sim.transitionRate(cPTN, dna, she, diffPTN)
            sim.transitionRate(cPTN, dna, ribo, diffPTN)
            sim.transitionRate(cPTN, ribo, dna, diffPTN)

#             trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)

            translat_rate = sim.rateConst(locusNum + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
    
            if max_poly_size > 1:

                aasequence_poly = aasequence[-134:]

                translat_rate_poly = sim.rateConst(locusNum + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)

            rna_deg_rate = sim.rateConst(locusNum + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)

            RNAP_on = sim.rateConst(locusNum + 'RNAP_on', RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)

            dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)

            dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)

#             dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)

    #         dna.addReaction([gene], [gene, RNA], trsc_rate)

            she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)

            she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)

    #         ribo.addReaction([RNA], [], rna_deg_rate)

#             ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)
        
#             Trsc_Term = sim.species('Trsc_Term')

#             dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)

    #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

#             ribo.addReaction([Ribo_RNA], [RNA, Ribosome, cPTN], translat_rate)
        
            for mp in range(Max_Poly_Count):

                Ribosome = sim.species('Ribosome'+str(mp+1))

                for pn in range(mp+1): #range(max_poly_size):

                    if (max_poly_size == 1) or (mp == 0):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)

                        ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

            #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                        ribo.addReaction([Ribo_RNA], [RNA, Ribosome, cPTN], translat_rate)

                        break

                    elif max_poly_size > 1:

                        if pn == 0:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

            #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                            Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)

                            Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                            ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, cPTN], translat_rate)

                        elif (pn+1 == max_poly_size) or (pn+1 == mp+1):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([Ribo_RNA], [RNA, Ribosome, cPTN], translat_rate_poly)

                            break

                        else:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)

                            Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                            ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, cPTN], translat_rate_poly)

            print(species)

#             ptnMetID = 'PTN_' + locusNum
#             PTN = sim.species(ptnMetID)

            secyID = 'S_' + locusNum
            secy_bound = sim.species(secyID)
            PartIdxMap[secyID] = partIdx
            partIdx = partIdx + 1
            pmap[secyID] = 0
            RDME_species_list.append(secyID)

            secy = sim.species('secY')

            # Get nucleotide and amino acid sequences, if available
#             rnasequence, aasequence = getSequences(jcvi3AID)

            transloc_rate = sim.rateConst(locusNum + '_insertion', TranslocRate(aasequence), 1)

            she.addReaction([secy,cPTN],[secy_bound],secY_on)
            she.addReaction([secy_bound],[secy,PTN],transloc_rate)

            secy_bound.diffusionRate(mem,sim.diffusionZero)
            secy_bound.diffusionRate(ext,sim.diffusionZero)
            secy_bound.diffusionRate(ribo,sim.diffusionZero)
            secy_bound.diffusionRate(dna,sim.diffusionZero)
            secy_bound.diffusionRate(cyt,sim.diffusionZero)
            secy_bound.diffusionRate(she,sim.diffusionZero)
            
            baseCount = defaultdict(int)
            for base in set(rnasequence):
                baseCount[base] = rnasequence.count(base)

            # Add total number of monomers to parameter dict

            N_A = baseCount["A"]

            N_U = baseCount["U"]

            N_C = baseCount["C"]

            N_G = baseCount["G"]

            gene_deg_dict = {}

            gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G

            gene_deg_dict['AMP_mRNAdeg'] = N_A
            gene_deg_dict['UMP_mRNAdeg'] = N_U
            gene_deg_dict['CMP_mRNAdeg'] = N_C
            gene_deg_dict['GMP_mRNAdeg'] = N_G

            degDict['D_' + locusNum] = gene_deg_dict


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

            translatATP = 0
            
            gene_ptn_dict = {}
            gene_ptn_dict['Type'] = 'membrane'
            gene_ptn_dict['States'] = [ptnMetID,ptnMetID + '_C',secyID]
            
            gene_translat_dict = {}

            for aaCost in AaUsed:

                aa_ID = aaCost[0]
                numberUsed = aaCost[1]

                aaCostID = aaCostMap[aa_ID]

                gene_translat_dict[aaCostID] = numberUsed

                translatATP = translatATP + 2*numberUsed


            gene_translat_dict['ATP_translat'] = translatATP
            
            gene_ptn_dict['Translat_costs'] = gene_translat_dict
            
            gene_ptn_dict['Transloc_cost'] = int(translatATP/20)

            multiStatePtnDict[ptnMetID] = gene_ptn_dict  

            
    else:
        
#         PTN = sim.species(ptnMetID)
#         sim.distributeNumber(PTN, sim.region('membrane'), int(ptnCount*ptn_ratio))
        ptnCount = ptnCount*ptnFrac*ptn_ratio
        pmap[ptnMetID] = ptnCount
        
#         PTN.diffusionRate(mem,sim.diffusionZero)
#         PTN.diffusionRate(ext,sim.diffusionZero)
#         PTN.diffusionRate(ribo,sim.diffusionZero)
#         PTN.diffusionRate(dna,sim.diffusionZero)
#         PTN.diffusionRate(she,diffPTN)
#         PTN.diffusionRate(mem,diffPTN)
        
    return partIdx


def addNamedPtnRDME(sim, pmap, newMetID, mmcode, transcribe, ptnFrac, ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, singleStatePtnDict, multiStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx):
    
    locusNum = mmcode.split('_')[1]
    
#     locusNum = jcvi3AID.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum

#     mmcode = 'MMSYN1_' + locusNum
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode

    print(mmcode, jcvi2ID, jcvi3AID)
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
#     ptnMetID = 'PTN_' + locusNum
    ptnMetID = newMetID
    
    # If the protein is not in the model, add it:

#     ModelSpecies.append(ptnMetID)

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    ptnCount = max(2,ptnCount)
    
#     ptnCount = ptnCount*ptnFrac
    
#     ptnCounts.append(ptnCount)

    print(ptnMetID, ptnCount)
    
    if transcribe:
        
        RDME_species_list.append(ptnMetID)
        
        genes_in_model_RDME.append(jcvi3AID)

        geneMetID = 'g_' + locusNum

#         ModelSpecies.append(geneMetID)

        # Get nucleotide and amino acid sequences, if available
        rnasequence, aasequence = getSequences(jcvi3AID)
    #     print(rnasequence)
    #     print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):
            
            aaCount = defaultdict(int)
            for aa in set(aasequence):
                aaCount[aa] = aasequence.count(aa)
            ptn_len = sum(list(aaCount.values()))

            if ptn_len <= 134:
                max_poly_size = 1
            elif ptn_len > 134:
                max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

            rnaMetID = "R_" + locusNum
            rnaName = "(mRNA) " + ptnName
            RDME_species_list.append(rnaMetID)
            
            pmap['New_mRNA_' + locusNum] = 0

#             ModelSpecies.append(rnaMetID)

#             species = []
#             species = [geneMetID, rnaMetID, ptnMetID]
#             gene = sim.species(geneMetID)
#             RNAP_gene = sim.species('RNAP_' + locusNum)
#             RNA = sim.species(rnaMetID)
#             Ribo_RNA = sim.species('Ribo_' + locusNum)
#             Deg_RNA = sim.species('Deg_' + locusNum)
#             PTN = sim.species(ptnMetID)
            
            species = []
            species = [geneMetID, rnaMetID, ptnMetID]
            gene = sim.species(geneMetID)

            RNAP_gene = sim.species('RP_' + locusNum)
            PartIdxMap['RP_' + locusNum] = partIdx
            partIdx = partIdx + 1
            pmap['RP_' + locusNum] = 0
            RDME_species_list.append('RP_' + locusNum)

            RNA = sim.species(rnaMetID)
            PartIdxMap[rnaMetID] = partIdx
            partIdx = partIdx + 1

            Deg_RNA = sim.species('D_' + locusNum)
            PartIdxMap['D_' + locusNum] = partIdx
            partIdx = partIdx + 1
            pmap['D_' + locusNum] = 0
            RDME_species_list.append('D_' + locusNum)

            PTN = sim.species(ptnMetID)
            PartIdxMap[ptnMetID] = partIdx
            partIdx = partIdx + 1

            for mp in range(Max_Poly_Count):

                for pn in range(mp+1): #range(max_poly_size):

                    if (max_poly_size == 1) or (mp == 0):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)

                        break

                    elif max_poly_size > 1:

                        if pn == 0:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                        elif (pn+1 == max_poly_size) or (pn == mp):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                            break

                        else:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)
            
#             cytoPtnMetID = ptnMetID + '_cyto'
#             cPTN = sim.species(cytoPtnMetID)

#             ModelSpecies.append('Ribo_' + locusNum)

            for index,row in mRNA_counts_DF.iterrows():
        
                if row["LocusTag"] == jcvi3AID:

                    avg_mRNA_cnt = row["Count"]

                    if avg_mRNA_cnt == 0.0:
                        avg_mRNA_cnt = 0.001

                    init_mRNA_count = np.random.poisson(avg_mRNA_cnt)

                    continue

            sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
            pmap[rnaMetID] = init_mRNA_count
            sim.distributeNumber(PTN, sim.region('cytoplasm'), int(4*ptnCount*ptn_ratio/5))
            sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))
            ptnCount = ptnCount*ptnFrac*ptn_ratio
            pmap[ptnMetID] = int(ptnCount)
#             sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))

            diffRNA = sim.diffusionConst('mrna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))
#             mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

            RNA.diffusionRate(cyt, diffRNA)
            PTN.diffusionRate(cyt, diffPTN)

            RNA.diffusionRate(mem,sim.diffusionZero)
            RNA.diffusionRate(ext,sim.diffusionZero)
            RNA.diffusionRate(ribo,diffRNA)
            RNA.diffusionRate(dna,diffRNA)
            RNA.diffusionRate(she,diffRNA)

            PTN.diffusionRate(mem,sim.diffusionZero)
            PTN.diffusionRate(ext,sim.diffusionZero)
            PTN.diffusionRate(ribo,diffPTN)
            PTN.diffusionRate(dna,diffPTN)
            PTN.diffusionRate(she,diffPTN)

            sim.transitionRate(RNA, dna, cyt, diffRNA)
            sim.transitionRate(RNA, cyt, dna, diffRNA)
            sim.transitionRate(RNA, cyt, ribo, diffRNA)
            sim.transitionRate(RNA, ribo, cyt, diffRNA)
            sim.transitionRate(RNA, ribo, she, diffRNA)
            sim.transitionRate(RNA, she, ribo, diffRNA)
            sim.transitionRate(RNA, she, cyt, diffRNA)
            sim.transitionRate(RNA, cyt, she, diffRNA)
            sim.transitionRate(RNA, she, dna, diffRNA)
            sim.transitionRate(RNA, dna, she, diffRNA)
            sim.transitionRate(RNA, dna, ribo, diffRNA)
            sim.transitionRate(RNA, ribo, dna, diffRNA)
            
            sim.transitionRate(PTN, dna, cyt, diffPTN)
            sim.transitionRate(PTN, cyt, dna, diffPTN)
    #         sim.transitionRate(PTN, cyt, ribo, diffPTN)
            sim.transitionRate(PTN, ribo, cyt, diffPTN)
            sim.transitionRate(PTN, ribo, she, diffPTN)
    #         sim.transitionRate(PTN, she, ribo, diffPTN)
            sim.transitionRate(PTN, she, cyt, diffPTN)
            sim.transitionRate(PTN, cyt, she, diffPTN)
            sim.transitionRate(PTN, she, dna, diffPTN)
            sim.transitionRate(PTN, dna, she, diffPTN)
            sim.transitionRate(PTN, dna, ribo, diffPTN)
            sim.transitionRate(PTN, ribo, dna, diffPTN)

#           
#             trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)

            translat_rate = sim.rateConst(locusNum + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
    
            if max_poly_size > 1:

                aasequence_poly = aasequence[-134:]

                translat_rate_poly = sim.rateConst(locusNum + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)

            rna_deg_rate = sim.rateConst(locusNum + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)

            RNAP_on = sim.rateConst(locusNum + 'RNAP_on', RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)

            dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)

            dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)

#             dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)

    #         dna.addReaction([gene], [gene, RNA], trsc_rate)

            she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)

            she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)

    #         ribo.addReaction([RNA], [], rna_deg_rate)

#             ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)
        
#             Trsc_Term = sim.species('Trsc_Term')

#             dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)

    #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

#             ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

            for mp in range(Max_Poly_Count):
        
                Ribosome = sim.species('Ribosome'+str(mp+1))

                for pn in range(mp+1): #range(max_poly_size):

                    if (max_poly_size == 1) or (mp == 0):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)

                        ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

            #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                        ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

                        break

                    elif max_poly_size > 1:

                        if pn == 0:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

            #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                            Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)

                            Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                            ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate)

                        elif (pn+1 == max_poly_size) or (pn+1 == mp+1):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate_poly)

                            break

                        else:

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)

                            Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                            ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate_poly)

            print(species)
            
            
            baseCount = defaultdict(int)
            for base in set(rnasequence):
                baseCount[base] = rnasequence.count(base)

            # Add total number of monomers to parameter dict

            N_A = baseCount["A"]

            N_U = baseCount["U"]

            N_C = baseCount["C"]

            N_G = baseCount["G"]

            gene_deg_dict = {}

            gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G

            gene_deg_dict['AMP_mRNAdeg'] = N_A
            gene_deg_dict['UMP_mRNAdeg'] = N_U
            gene_deg_dict['CMP_mRNAdeg'] = N_C
            gene_deg_dict['GMP_mRNAdeg'] = N_G

            degDict['D_' + locusNum] = gene_deg_dict


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

            translatATP = 0
            
            if ptnFrac == 1.0:
                
                print(ptnMetID,ptnFrac)
                
                gene_translat_dict = {}
                
                for aaCost in AaUsed:
            
                    aa_ID = aaCost[0]
                    numberUsed = aaCost[1]

                    aaCostID = aaCostMap[aa_ID]

                    gene_translat_dict[aaCostID] = numberUsed

                    translatATP = translatATP + 2*numberUsed


                gene_translat_dict['ATP_translat'] = translatATP

                singleStatePtnDict[ptnMetID] = gene_translat_dict
                
            else:
            
                gene_ptn_dict = {}
                gene_ptn_dict['Type'] = 'cytoplasm'
                gene_ptn_dict['States'] = [ptnMetID]

                gene_translat_dict = {}

                for aaCost in AaUsed:

                    aa_ID = aaCost[0]
                    numberUsed = aaCost[1]

                    aaCostID = aaCostMap[aa_ID]

                    gene_translat_dict[aaCostID] = numberUsed

                    translatATP = translatATP + 2*numberUsed


                gene_translat_dict['ATP_translat'] = translatATP

                gene_ptn_dict['Translat_costs'] = gene_translat_dict

                multiStatePtnDict[ptnMetID] = gene_ptn_dict
                

    else:
        
#         PTN = sim.species(ptnMetID)
#         sim.distributeNumber(PTN, sim.region('membrane'), int(ptnCount*ptn_ratio))
        ptnCount = ptnCount*ptnFrac*ptn_ratio
        pmap[ptnMetID] = int(ptnCount)
        
#         PTN.diffusionRate(mem,sim.diffusionZero)
#         PTN.diffusionRate(ext,sim.diffusionZero)
#         PTN.diffusionRate(ribo,sim.diffusionZero)
#         PTN.diffusionRate(dna,sim.diffusionZero)
#         PTN.diffusionRate(she,diffPTN)
#         PTN.diffusionRate(mem,diffPTN)
        
    return partIdx


def addRiboPtnRDME(sim, pmap, newMetID, mmcode, ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, singleStatePtnDict, aaCostMap, RDME_species_list, PartIdxMap, partIdx):
    
    locusNum = mmcode.split('_')[1]
    
#     locusNum = jcvi3AID.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum

#     mmcode = 'MMSYN1_' + locusNum
    # Checks if a translation to JCVISYN2* code is available
    try:
        jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
    except:
        jcvi2ID = "JCVIunk_" + mmcode

    print(mmcode, jcvi2ID, jcvi3AID)
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
#     ptnMetID = 'PTN_' + locusNum
    ptnMetID = newMetID
    RDME_species_list.append(ptnMetID)
    
    # If the protein is not in the model, add it:

#     ModelSpecies.append(ptnMetID)

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
    
    ptnCount = max(2,ptnCount)
    
#     ptnCount = ptnCount*ptnFrac
    
#     ptnCounts.append(ptnCount)

    print(ptnMetID, ptnCount)
    
    genes_in_model_RDME.append(jcvi3AID)

    geneMetID = 'g_' + locusNum

#     ModelSpecies.append(geneMetID)

    # Get nucleotide and amino acid sequences, if available
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):
        
        aaCount = defaultdict(int)
        for aa in set(aasequence):
            aaCount[aa] = aasequence.count(aa)
        ptn_len = sum(list(aaCount.values()))
        
        if ptn_len <= 134:
            max_poly_size = 1
        elif ptn_len > 134:
            max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

        rnaMetID = "R_" + locusNum
        rnaName = "(mRNA) " + ptnName
        RDME_species_list.append(rnaMetID)
        
        pmap['New_mRNA_' + locusNum] = 0

#         ModelSpecies.append(rnaMetID)

#         species = []
#         species = [geneMetID, rnaMetID, ptnMetID]
#         gene = sim.species(geneMetID)
#         RNAP_gene = sim.species('RNAP_' + locusNum)
#         RNA = sim.species(rnaMetID)
#         Ribo_RNA = sim.species('Ribo_' + locusNum)
#         Deg_RNA = sim.species('Deg_' + locusNum)
#         PTN = sim.species(ptnMetID)
        
        species = []
        species = [geneMetID, rnaMetID, ptnMetID]
        gene = sim.species(geneMetID)

        RNAP_gene = sim.species('RP_' + locusNum)
        PartIdxMap['RP_' + locusNum] = partIdx
        partIdx = partIdx + 1
        pmap['RP_' + locusNum] = 0
        RDME_species_list.append('RP_' + locusNum)

        RNA = sim.species(rnaMetID)
        PartIdxMap[rnaMetID] = partIdx
        partIdx = partIdx + 1

        Deg_RNA = sim.species('D_' + locusNum)
        PartIdxMap['D_' + locusNum] = partIdx
        partIdx = partIdx + 1
        pmap['D_' + locusNum] = 0
        RDME_species_list.append('D_' + locusNum)

        PTN = sim.species(ptnMetID)
        PartIdxMap[ptnMetID] = partIdx
        partIdx = partIdx + 1

        for mp in range(Max_Poly_Count):

            for pn in range(mp+1): #range(max_poly_size):

                if (max_poly_size == 1) or (mp == 0):

                    Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                    Ribo_RNA = sim.species(Ribo_ID)
                    PartIdxMap[Ribo_ID] = partIdx
                    partIdx = partIdx + 1
                    pmap[Ribo_ID] = 0
                    RDME_species_list.append(Ribo_ID)

                    break

                elif max_poly_size > 1:

                    if pn == 0:

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)

                    elif (pn+1 == max_poly_size) or (pn == mp):

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)

                        break

                    else:

                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)

                        Ribo_RNA = sim.species(Ribo_ID)
                        PartIdxMap[Ribo_ID] = partIdx
                        partIdx = partIdx + 1
                        pmap[Ribo_ID] = 0
                        RDME_species_list.append(Ribo_ID)

#         cytoPtnMetID = ptnMetID + '_cyto'
#         cPTN = sim.species(cytoPtnMetID)

#         ModelSpecies.append('Ribo_' + locusNum)

        for index,row in mRNA_counts_DF.iterrows():
        
            if row["LocusTag"] == jcvi3AID:

                avg_mRNA_cnt = row["Count"]

                if avg_mRNA_cnt == 0.0:
                    avg_mRNA_cnt = 0.001

                init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
                continue
        
        sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
        pmap[rnaMetID] = init_mRNA_count
        sim.distributeNumber(PTN, sim.region('cytoplasm'), int(4*ptnCount*ptn_ratio/5))
        sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))
        ptnCount = ptnCount*ptn_ratio
        pmap[ptnMetID] = int(ptnCount)
#             sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))

        diffRNA = sim.diffusionConst('mrna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))
#         mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

        RNA.diffusionRate(cyt, diffRNA)
        PTN.diffusionRate(cyt, diffPTN)

        RNA.diffusionRate(mem,sim.diffusionZero)
        RNA.diffusionRate(ext,sim.diffusionZero)
        RNA.diffusionRate(ribo,diffRNA)
        RNA.diffusionRate(dna,diffRNA)
        RNA.diffusionRate(she,diffRNA)

        PTN.diffusionRate(mem,sim.diffusionZero)
        PTN.diffusionRate(ext,sim.diffusionZero)
        PTN.diffusionRate(ribo,diffPTN)
        PTN.diffusionRate(dna,diffPTN)
        PTN.diffusionRate(she,diffPTN)

        sim.transitionRate(RNA, dna, cyt, diffRNA)
        sim.transitionRate(RNA, cyt, dna, diffRNA)
        sim.transitionRate(RNA, cyt, ribo, diffRNA)
        sim.transitionRate(RNA, ribo, cyt, diffRNA)
        sim.transitionRate(RNA, ribo, she, diffRNA)
        sim.transitionRate(RNA, she, ribo, diffRNA)
        sim.transitionRate(RNA, she, cyt, diffRNA)
        sim.transitionRate(RNA, cyt, she, diffRNA)
        sim.transitionRate(RNA, she, dna, diffRNA)
        sim.transitionRate(RNA, dna, she, diffRNA)
        sim.transitionRate(RNA, dna, ribo, diffRNA)
        sim.transitionRate(RNA, ribo, dna, diffRNA)
        
        sim.transitionRate(PTN, dna, cyt, diffPTN)
        sim.transitionRate(PTN, cyt, dna, diffPTN)
#         sim.transitionRate(PTN, cyt, ribo, diffPTN)
        sim.transitionRate(PTN, ribo, cyt, diffPTN)
        sim.transitionRate(PTN, ribo, she, diffPTN)
#         sim.transitionRate(PTN, she, ribo, diffPTN)
        sim.transitionRate(PTN, she, cyt, diffPTN)
        sim.transitionRate(PTN, cyt, she, diffPTN)
        sim.transitionRate(PTN, she, dna, diffPTN)
        sim.transitionRate(PTN, dna, she, diffPTN)
        sim.transitionRate(PTN, dna, ribo, diffPTN)
        sim.transitionRate(PTN, ribo, dna, diffPTN)

#             trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)

        translat_rate = sim.rateConst(locusNum + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
    
        if max_poly_size > 1:
            
            aasequence_poly = aasequence[-134:]
        
            translat_rate_poly = sim.rateConst(locusNum + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)

        rna_deg_rate = sim.rateConst(locusNum + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)

        RNAP_on = sim.rateConst(locusNum + 'RNAP_on', RNAP_binding_ribo(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)

        dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)

        dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)

#             dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)

#         dna.addReaction([gene], [gene, RNA], trsc_rate)

        she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)

        she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)

#         ribo.addReaction([RNA], [], rna_deg_rate)

#         ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

#         Trsc_Term = sim.species('Trsc_Term')

#         dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)

#         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

#         ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

        for mp in range(Max_Poly_Count):
        
            Ribosome = sim.species('Ribosome'+str(mp+1))
            
            for pn in range(mp+1): #range(max_poly_size):
                
                if (max_poly_size == 1) or (mp == 0):
                    
                    Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                    Ribo_RNA = sim.species(Ribo_ID)
                    
                    ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                    ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)
                    
                    break
                
                elif max_poly_size > 1:
                    
                    if pn == 0:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)
                        
                        ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)
        
                        Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)
                
                        Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                        ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate)
            
                    elif (pn+1 == max_poly_size) or (pn+1 == mp+1):
                    
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)

                        ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate_poly)
                        
                        break
                        
                    else:
                        
                        Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+1)
                
                        Ribo_RNA = sim.species(Ribo_ID)
                        
                        Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum + '_' + str(pn+2)
                
                        Ribo_RNA_Next = sim.species(Ribo_ID_Next)
                        
                        ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate_poly)

        print(species)
        
        baseCount = defaultdict(int)
        for base in set(rnasequence):
            baseCount[base] = rnasequence.count(base)

        # Add total number of monomers to parameter dict

        N_A = baseCount["A"]

        N_U = baseCount["U"]

        N_C = baseCount["C"]

        N_G = baseCount["G"]
        
        gene_deg_dict = {}
        
        gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G
        
        gene_deg_dict['AMP_mRNAdeg'] = N_A
        gene_deg_dict['UMP_mRNAdeg'] = N_U
        gene_deg_dict['CMP_mRNAdeg'] = N_C
        gene_deg_dict['GMP_mRNAdeg'] = N_G
        
        degDict['D_' + locusNum] = gene_deg_dict
        
        
        gene_translat_dict = {}
        
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
        
        translatATP = 0
        
        for aaCost in AaUsed:
            
            aa_ID = aaCost[0]
            numberUsed = aaCost[1]

            aaCostID = aaCostMap[aa_ID]
            
            gene_translat_dict[aaCostID] = numberUsed
            
            translatATP = translatATP + 2*numberUsed
            
            
        gene_translat_dict['ATP_translat'] = translatATP
        
        singleStatePtnDict[ptnMetID] = gene_translat_dict

            
#     else:
        
#         PTN = sim.species(ptnMetID)
# #         sim.distributeNumber(PTN, sim.region('membrane'), int(ptnCount*ptn_ratio))
#         ptnCount = ptnCount*ptnFrac*ptn_ratio
#         pmap[ptnMetID] = ptnCount
        
#         PTN.diffusionRate(mem,sim.diffusionZero)
#         PTN.diffusionRate(ext,sim.diffusionZero)
#         PTN.diffusionRate(ribo,sim.diffusionZero)
#         PTN.diffusionRate(dna,sim.diffusionZero)
#         PTN.diffusionRate(she,diffPTN)
#         PTN.diffusionRate(mem,diffPTN)
        
    return partIdx

##################
def addRNAPpart(sim, pmap, ptn_ratio, RNApol, ext, mem, cyt, ribo, dna, she):
    diffRNApol = sim.diffusionConst("RNAP_diff", 0.22e-12)

    sim.distributeNumber(RNApol, sim.region('cytoplasm'), 187*ptn_ratio)
    pmap['RNApol'] = int(187*ptn_ratio)

    RNApol.diffusionRate(mem,sim.diffusionZero)
    RNApol.diffusionRate(ext,sim.diffusionZero)
    RNApol.diffusionRate(ribo,sim.diffusionZero)
    RNApol.diffusionRate(dna,diffRNApol)
    RNApol.diffusionRate(cyt,diffRNApol)
    RNApol.diffusionRate(she,diffRNApol)

    sim.transitionRate(RNApol, cyt, dna, diffRNApol)
    sim.transitionRate(RNApol, dna, cyt, diffRNApol)
    sim.transitionRate(RNApol, cyt, she, diffRNApol)
    sim.transitionRate(RNApol, she, cyt, diffRNApol)
    sim.transitionRate(RNApol, she, dna, diffRNApol)
    sim.transitionRate(RNApol, dna, she, diffRNApol)
    
    print('Added RNAP particles')
    
    return sim
##################


##################
def addRiboPart(sim, pmap, partIdx, PartIdxMap, ext, mem, cyt, ribo, dna, she, ribo_center_points, RDME_species_list):
    
    polysomes = 0

    poly_ribo = []
    poly_dict = {}

    for i in range(ribo_center_points.shape[0]-1):

    #     poly_cluster = [i]
        poly_dict[str(i)] = []

        for j in range(i+1,ribo_center_points.shape[0]):

    #         dist = np.sum(np.abs(ribo_center_points[i,:]-ribo_center_points[j,:]))

            coord1 = ribo_center_points[i]
            coord2 = ribo_center_points[j]

            x1 = int(coord1[0])
            y1 = int(coord1[1])
            z1 = int(coord1[2])

            x2 = int(coord2[0])
            y2 = int(coord2[1])
            z2 = int(coord2[2])

            radius = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**(1/2)
            if radius < 2.3:
    #             print(i,j)
                polysomes = polysomes + 1

                poly_ribo.append([i,j])
                poly_dict[str(i)].append(j)

    print(polysomes)
    print(poly_ribo)
    
    for cluster in poly_ribo:
    
        complete_cluster = False

        while complete_cluster == False:

            cluster_size = len(cluster)

            for pair in poly_ribo:

                if (pair != cluster):

    #                 if pair[0] == cluster[0]:

    #                     print(cluster)
    #                     cluster.append(pair[1])
    #                     poly_ribo.remove(pair)
    #                     print(cluster)

    #                 else:

    #                 print(cluster)

                    for ribo1 in cluster:

                        if ribo1 in pair:

                            for ribo2 in pair:

                                if ribo2 not in cluster:

                                    cluster.append(ribo2)
    #                                 print(cluster)

    #                 poly_ribo.remove(pair)

            if len(cluster) == cluster_size:

                complete_cluster = True

    for cluster1 in poly_ribo:
    
        for cluster2 in poly_ribo:

            if cluster1 != cluster2:

                same_cluster_check = all(ribo in cluster1 for ribo in cluster2)

                if same_cluster_check:
    #                 print(same_cluster_check,cluster1,cluster2)
                    poly_ribo.remove(cluster2)
        
    total_poly_clust = 0
    single_poly_ribos = []

    for cluster in poly_ribo:

        total_poly_clust = total_poly_clust + len(cluster)

        for riboP in cluster:

            if riboP not in single_poly_ribos:

                single_poly_ribos.append(ribo)

    print(total_poly_clust)
    print(len(single_poly_ribos))
    print(poly_ribo)
    
    print('Polysomes assigned')
    
    ribos = 0
    
    for i in range(len(ribo_center_points)):
    
        Ribosome = sim.species('Ribosome1')

        for cluster in poly_ribo:

            if i in cluster:

                poly_size = len(cluster)

                Ribosome = sim.species('Ribosome'+str(poly_size))
                
                if 'Ribosome'+str(poly_size) not in PartIdxMap:
                    
                    PartIdxMap['Ribosome'+str(poly_size)] = partIdx
                    partIdx = partIdx + 1
                    RDME_species_list.append('Ribosome'+str(poly_size))
                    
                    print('Particle Index for Ribosome' + str(poly_size) + ' = ' + str(partIdx))

        coord = ribo_center_points[i]

        x = int(coord[0])
        y = int(coord[1])
        z = int(coord[2])

        Ribosome.placeParticle(x,y,z,1)

        ribos += 1
        
#     pmap['Ribosome'] = ribos

    print('Added ribosome particles')
    return sim, poly_ribo, partIdx
##################


##################

# Define how to add the particles and reactions for each tRNA and its corresponding gene.

def addtRNA(sim, pmap, unchargedMetID, chargedMetID, mmcode, aminoAcid, synthase, ptn_ratio, ext, mem, cyt, ribo, dna, she, genes_in_model_RDME, mRNA_counts_DF, diffPTN, diffPTNslow, RNAP_on, RNAP_off, Ribo_on, Ribo_off, deg_bind_rate, RNApol, poly_ribo, Max_Poly_Count, degradosome, TrscTermination, degDict, multiStatePtnDict, aaCostMap, tRNAstateDict, tRNAadded, RDME_species_list, PartIdxMap, partIdx, rtRNA_ID_dict):
    
    locusNum = mmcode.split('_')[1]

    jcvi3AID = 'JCVISYN3A_' + locusNum
    
    locusNum = jcvi3AID.split('_')[1].lstrip('0')
    
    new_RNA_ID = 'New_RNA_' + locusNum
    pmap[new_RNA_ID] = 0
    
    rtRNA_ID_dict[new_RNA_ID] = unchargedMetID

    # Checks if a translation to JCVISYN2* code is available
#     try:
#         jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode ].iloc[0, 13].strip()
#     except:
#         # If not, check for a manual AOE connection:
#         if mmcode in manGPRPD.MM.values:
#             jcvi2ID = "JCVIman_" + mmcode
#         else:
#             # If not, set an "unknown" code
#             jcvi2ID = "JCVIunk_" + mmcode + "_" + str(unkIter)
#             unkIter += 1 

#     print(jcvi3AID)

    geneMetID = 'g_' + locusNum

    rnasequence = getRNAsequences(jcvi3AID)
#     print(rnasequence)

#     gene = [geneMetID]
#     print(gene)
#     sim.defineSpecies(gene)

#     sim.addParticles(species = geneMetID, count = 1)

#     charged_tRNA = sim.species(chargedMetID)
    tRNA = sim.species(unchargedMetID)
    
    if unchargedMetID not in tRNAadded:
        
        print('Adding '+unchargedMetID)
        
        RDME_species_list.append(unchargedMetID)
        PartIdxMap[unchargedMetID] = partIdx
        partIdx = partIdx + 1
        
        diffRNA = sim.diffusionConst('trna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))
        
        if (unchargedMetID == 'M_trnaleu_c') or (unchargedMetID == 'M_trnamet_c'):
            
            sim.distributeNumber(tRNA, sim.region('cytoplasm'), int(3*round(5800/29)))
        
        elif (unchargedMetID == 'M_trnathr_c') or (unchargedMetID == 'M_trnatrp_c') or (unchargedMetID == 'M_trnalys_c') or (unchargedMetID == 'M_trnaarg_c') or (unchargedMetID == 'M_trnaser_c'):
            
            sim.distributeNumber(tRNA, sim.region('cytoplasm'), int(2*round(5800/29)))
            
        else:
            
            sim.distributeNumber(tRNA, sim.region('cytoplasm'), int(round(5800/29)))

        tRNA.diffusionRate(mem,sim.diffusionZero)
        tRNA.diffusionRate(ext,sim.diffusionZero)
        tRNA.diffusionRate(cyt,diffRNA)
        tRNA.diffusionRate(ribo,diffRNA)
        tRNA.diffusionRate(dna,diffRNA)
        tRNA.diffusionRate(she,diffRNA)
        
        sim.transitionRate(tRNA, dna, cyt, diffRNA)
        sim.transitionRate(tRNA, cyt, dna, diffRNA)
        sim.transitionRate(tRNA, cyt, ribo, diffRNA)
        sim.transitionRate(tRNA, ribo, cyt, diffRNA)
        sim.transitionRate(tRNA, ribo, she, diffRNA)
        sim.transitionRate(tRNA, she, ribo, diffRNA)
        sim.transitionRate(tRNA, she, cyt, diffRNA)
        sim.transitionRate(tRNA, cyt, she, diffRNA)
        sim.transitionRate(tRNA, she, dna, diffRNA)
        sim.transitionRate(tRNA, dna, she, diffRNA)
        sim.transitionRate(tRNA, dna, ribo, diffRNA)
        sim.transitionRate(tRNA, ribo, dna, diffRNA)

        tRNAadded.append(unchargedMetID)
        tRNAadded.append(chargedMetID)
        
        print(synthase)
    
        if 'GLN' not in synthase:
            
            locusNum_synth = synthase.split('3A_')[1] #jcvi3AID.split('_')[1]
            
            mmcode_synth = 'MMSYN1_' + locusNum_synth
            jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            
            synthase = synthase.split('3A_')[1].lstrip('0')
            synthase_ptn = 'P_' + synthase
            
            locusNum_synth = synthase #jcvi3AID.split('_')[1]
            
            print(synthase,synthase_ptn)
    
            synthase_atp = synthase_ptn + '_ATP'
            if synthase_atp not in pmap:
                pmap[synthase_atp] = 0
            
            synthase_atp_aa = synthase_atp + '_AA'
            if synthase_atp_aa not in pmap:
                pmap[synthase_atp_aa] = 0
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNA'
            if synthase_atp_aa_trna not in pmap:
                pmap[synthase_atp_aa_trna] = 0
                
            tRNAstateDict[unchargedMetID] = [unchargedMetID,chargedMetID,synthase_atp_aa_trna]
                

            # Checks if a translation to JCVISYN2* code is available
            try:
                jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode_synth ].iloc[0, 13].strip()
            except:
                jcvi2ID = "JCVIunk_" + mmcode_synth

        #     print(mmcode, jcvi2ID, jcvi3AID)

            genes_in_model_RDME.append(jcvi3AID_synth)

            # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
            ptnMetID = 'P_' + locusNum_synth
            RDME_species_list.append(ptnMetID)

            # If the protein is not in the model, add it:

            ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)
            
            ptnCount = max(2,ptnCount)

        #     ptnCounts.append(ptnCount)

        #     print(ptnMetID, ptnCount)

            geneMetID = 'g_' + locusNum_synth

            # Get nucleotide and amino acid sequences, if available
            print(jcvi3AID_synth)
            rnasequence, aasequence = getSequences(jcvi3AID_synth)
            print(rnasequence)
            print(aasequence)

            if (rnasequence != 0) and (aasequence != 0):
                
                aaCount = defaultdict(int)
                for aa in set(aasequence):
                    aaCount[aa] = aasequence.count(aa)
                ptn_len = sum(list(aaCount.values()))

                if ptn_len <= 134:
                    max_poly_size = 1
                elif ptn_len > 134:
                    max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

                rnaMetID = "R_" + locusNum_synth
                rnaName = "(mRNA) " + ptnName
                RDME_species_list.append(rnaMetID)
                
                pmap['New_mRNA_' + locusNum_synth] = 0

                species = []
                species = [geneMetID, rnaMetID, ptnMetID]
                gene = sim.species(geneMetID)

                RNAP_gene = sim.species('RP_' + locusNum_synth)
                PartIdxMap['RP_' + locusNum_synth] = partIdx
                partIdx = partIdx + 1
                pmap['RP_' + locusNum_synth] = 0
                RDME_species_list.append('RP_' + locusNum_synth)

                RNA = sim.species(rnaMetID)
                PartIdxMap[rnaMetID] = partIdx
                partIdx = partIdx + 1

                Deg_RNA = sim.species('D_' + locusNum_synth)
                PartIdxMap['D_' + locusNum_synth] = partIdx
                partIdx = partIdx + 1
                pmap['D_' + locusNum_synth] = 0
                RDME_species_list.append('D_' + locusNum_synth)

                PTN = sim.species(ptnMetID)
                PartIdxMap[ptnMetID] = partIdx
                partIdx = partIdx + 1
                
                print('Added stuff')

                for mp in range(Max_Poly_Count):

                    for pn in range(mp+1): #range(max_poly_size):

                        if (max_poly_size == 1) or (mp == 0):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                            break

                        elif max_poly_size > 1:

                            if pn == 0:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                            elif (pn+1 == max_poly_size) or (pn == mp):

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                                break

                            else:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                for index,row in mRNA_counts_DF.iterrows():

                    if row["LocusTag"] == jcvi3AID_synth:

                        avg_mRNA_cnt = row["Count"]

                        if avg_mRNA_cnt == 0.0:
                            avg_mRNA_cnt = 0.001

                        init_mRNA_count = np.random.poisson(avg_mRNA_cnt)

                        continue

                sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
                pmap[rnaMetID] = init_mRNA_count
                sim.distributeNumber(PTN, sim.region('cytoplasm'), int(4*ptnCount*ptn_ratio/5))
                sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))
                pmap[ptnMetID] = int(ptnCount*ptn_ratio)

                diffRNA = sim.diffusionConst('mrna_' + locusNum_synth + '_diff', RNA_diff_coeff(rnasequence))
        #         mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

                RNA.diffusionRate(cyt, diffRNA)
                PTN.diffusionRate(cyt, diffPTN)

                RNA.diffusionRate(mem,sim.diffusionZero)
                RNA.diffusionRate(ext,sim.diffusionZero)
                RNA.diffusionRate(ribo,diffRNA)
                RNA.diffusionRate(dna,diffRNA)
                RNA.diffusionRate(she,diffRNA)

                PTN.diffusionRate(mem,sim.diffusionZero)
                PTN.diffusionRate(ext,sim.diffusionZero)
                PTN.diffusionRate(ribo,diffPTN)
                PTN.diffusionRate(dna,diffPTN)
                PTN.diffusionRate(she,diffPTN)

                sim.transitionRate(RNA, dna, cyt, diffRNA)
                sim.transitionRate(RNA, cyt, dna, diffRNA)
                sim.transitionRate(RNA, cyt, ribo, diffRNA)
                sim.transitionRate(RNA, ribo, cyt, diffRNA)
                sim.transitionRate(RNA, ribo, she, diffRNA)
                sim.transitionRate(RNA, she, ribo, diffRNA)
                sim.transitionRate(RNA, she, cyt, diffRNA)
                sim.transitionRate(RNA, cyt, she, diffRNA)
                sim.transitionRate(RNA, she, dna, diffRNA)
                sim.transitionRate(RNA, dna, she, diffRNA)
                sim.transitionRate(RNA, dna, ribo, diffRNA)
                sim.transitionRate(RNA, ribo, dna, diffRNA)

                sim.transitionRate(PTN, dna, cyt, diffPTN)
                sim.transitionRate(PTN, cyt, dna, diffPTN)
        #         sim.transitionRate(PTN, cyt, ribo, diffPTN)
                sim.transitionRate(PTN, ribo, cyt, diffPTN)
                sim.transitionRate(PTN, ribo, she, diffPTN)
        #         sim.transitionRate(PTN, she, ribo, diffPTN)
                sim.transitionRate(PTN, she, cyt, diffPTN)
                sim.transitionRate(PTN, cyt, she, diffPTN)
                sim.transitionRate(PTN, she, dna, diffPTN)
                sim.transitionRate(PTN, dna, she, diffPTN)
                sim.transitionRate(PTN, dna, ribo, diffPTN)
                sim.transitionRate(PTN, ribo, dna, diffPTN)
                
                print('Diffusion Coefficients')

        #         trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)

                translat_rate = sim.rateConst(locusNum_synth + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
            
                if max_poly_size > 1:
            
                    aasequence_poly = aasequence[-134:]

                    translat_rate_poly = sim.rateConst(locusNum_synth + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)

                rna_deg_rate = sim.rateConst(locusNum_synth + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)

                RNAP_on = sim.rateConst(locusNum_synth + 'RNAP_on', RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)

                dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)

                dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)

        #         dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)

        #         dna.addReaction([gene], [gene, RNA], trsc_rate)

                she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)

                she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)

        #         ribo.addReaction([RNA], [], rna_deg_rate)

#                 ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

#                 ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

                for mp in range(Max_Poly_Count):
        
                    Ribosome = sim.species('Ribosome'+str(mp+1))

                    for pn in range(mp+1): #range(max_poly_size):

                        if (max_poly_size == 1) or (mp == 0):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

                #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                            ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

                            break

                        elif max_poly_size > 1:

                            if pn == 0:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

                #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                                Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+2)

                                Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                                ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate)

                            elif (pn+1 == max_poly_size) or (pn+1 == mp+1):

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate_poly)

                                break

                            else:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+2)

                                Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                                ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate_poly)
                                
#                 Trsc_Term = sim.species('Trsc_Term')

#                 dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)

                print('Added rxns')
    
                baseCount = defaultdict(int)
                for base in set(rnasequence):
                    baseCount[base] = rnasequence.count(base)

                # Add total number of monomers to parameter dict

                N_A = baseCount["A"]

                N_U = baseCount["U"]

                N_C = baseCount["C"]

                N_G = baseCount["G"]

                gene_deg_dict = {}

                gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G

                gene_deg_dict['AMP_mRNAdeg'] = N_A
                gene_deg_dict['UMP_mRNAdeg'] = N_U
                gene_deg_dict['CMP_mRNAdeg'] = N_C
                gene_deg_dict['GMP_mRNAdeg'] = N_G

                degDict['D_' + locusNum_synth] = gene_deg_dict
                
                
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

                translatATP = 0

                gene_ptn_dict = {}
                gene_ptn_dict['Type'] = 'cytoplasm'
                gene_ptn_dict['States'] = [ptnMetID, synthase_atp, synthase_atp_aa, synthase_atp_aa_trna]

                gene_translat_dict = {}

                for aaCost in AaUsed:

                    aa_ID = aaCost[0]
                    numberUsed = aaCost[1]

                    aaCostID = aaCostMap[aa_ID]

                    gene_translat_dict[aaCostID] = numberUsed

                    translatATP = translatATP + 2*numberUsed


                gene_translat_dict['ATP_translat'] = translatATP

                gene_ptn_dict['Translat_costs'] = gene_translat_dict
                
                print(ptnMetID,gene_ptn_dict)

                multiStatePtnDict[ptnMetID] = gene_ptn_dict 

        
        elif 'GLN' in synthase:
            
            locusNum_synth = '0126' #jcvi3AID.split('_')[1]
            
            mmcode_synth = 'MMSYN1_' + locusNum_synth
            jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            
            locusNum_synth = '126' #jcvi3AID.split('_')[1]
            
            synthase_ptn = 'P_126'
    
            synthase_atp = synthase_ptn + '_ATP'
            if synthase_atp not in pmap:
                pmap[synthase_atp] = 0
            
            synthase_atp_aa = synthase_atp + '_AA'
            if synthase_atp_aa not in pmap:
                pmap[synthase_atp_aa] = 0
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNAgln'
            if synthase_atp_aa_trna not in pmap:
                pmap[synthase_atp_aa_trna] = 0
                
                
#             locusNum_synth = '126' #jcvi3AID.split('_')[1]
#             mmcode_synth = 'MMSYN1_' + locusNum_synth
#             jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            
            glugln = 'M_glutrnagln_c'
            glugln_enz = 'glutrnagln_enz'
            glugln_enz_atp = 'glutrnagln_enz_atp'
            glugln_enz_atp_aa = 'glutrnagln_enz_atp_gln'
            
            pmap[glugln] = 0
            pmap[glugln_enz] = 0
            pmap[glugln_enz_atp] = 0
            pmap[glugln_enz_atp_aa] = 0
            
            tRNAstateDict[unchargedMetID] = [unchargedMetID,chargedMetID,synthase_atp_aa_trna,glugln,glugln_enz,glugln_enz_atp,glugln_enz_atp_aa]
            
#             locusNum_synth = synthase #jcvi3AID.split('_')[1]
#             mmcode_synth = 'MMSYN1_' + locusNum_synth
#             jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            # Checks if a translation to JCVISYN2* code is available
    
#             locusNum_synth = '689' #jcvi3AID.split('_')[1]
#             mmcode_synth = 'MMSYN1_' + locusNum_synth
#             jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            
            locusNum_synth = '0689' #jcvi3AID.split('_')[1]
            
            mmcode_synth = 'MMSYN1_' + locusNum_synth
            jcvi3AID_synth = 'JCVISYN3A_' + locusNum_synth
            
            locusNum_synth = '689'
        
            try:
                jcvi2ID = annotatPD.loc[ annotatPD.iloc[:,5] == mmcode_synth ].iloc[0, 13].strip()
            except:
                jcvi2ID = "JCVIunk_" + mmcode_synth

        #     print(mmcode, jcvi2ID, jcvi3AID)

            genes_in_model_RDME.append(jcvi3AID_synth)

            # We name proteins after their locus tag from the gen bank entry to be M_PTN_JCVISYN3A_XXXX_c.
            ptnMetID = 'P_' + locusNum_synth
            RDME_species_list.append(ptnMetID)

            # If the protein is not in the model, add it:

            ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)

        #     ptnCounts.append(ptnCount)

        #     print(ptnMetID, ptnCount)

            geneMetID = 'g_' + locusNum_synth

            # Get nucleotide and amino acid sequences, if available
            rnasequence, aasequence = getSequences(jcvi3AID_synth)
        #     print(rnasequence)
        #     print(aasequence)

            if (rnasequence != 0) and (aasequence != 0):
                
                aaCount = defaultdict(int)
                for aa in set(aasequence):
                    aaCount[aa] = aasequence.count(aa)
                ptn_len = sum(list(aaCount.values()))

                if ptn_len <= 134:
                    max_poly_size = 1
                elif ptn_len > 134:
                    max_poly_size = min(Max_Poly_Count,int(ptn_len/134))

                rnaMetID = "R_" + locusNum_synth
                rnaName = "(mRNA) " + ptnName
                RDME_species_list.append(rnaMetID)
                
                pmap['New_mRNA_' + locusNum_synth] = 0

                species = []
                species = [geneMetID, rnaMetID, ptnMetID]
                gene = sim.species(geneMetID)

                RNAP_gene = sim.species('RP_' + locusNum_synth)
                PartIdxMap['RP_' + locusNum_synth] = partIdx
                partIdx = partIdx + 1
                pmap['RP_' + locusNum_synth] = 0
                RDME_species_list.append('RP_' + locusNum_synth)

                RNA = sim.species(rnaMetID)
                PartIdxMap[rnaMetID] = partIdx
                partIdx = partIdx + 1

                Deg_RNA = sim.species('D_' + locusNum_synth)
                PartIdxMap['D_' + locusNum_synth] = partIdx
                partIdx = partIdx + 1
                pmap['D_' + locusNum_synth] = 0
                RDME_species_list.append('D_' + locusNum_synth)

                PTN = sim.species(ptnMetID)
                PartIdxMap[ptnMetID] = partIdx
                partIdx = partIdx + 1

                for mp in range(Max_Poly_Count):

                    for pn in range(mp+1): #range(max_poly_size):

                        if (max_poly_size == 1) or (mp == 0):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)
                            PartIdxMap[Ribo_ID] = partIdx
                            partIdx = partIdx + 1
                            pmap[Ribo_ID] = 0
                            RDME_species_list.append(Ribo_ID)

                            break

                        elif max_poly_size > 1:

                            if pn == 0:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                            elif (pn+1 == max_poly_size) or (pn == mp):

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                                break

                            else:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)
                                PartIdxMap[Ribo_ID] = partIdx
                                partIdx = partIdx + 1
                                pmap[Ribo_ID] = 0
                                RDME_species_list.append(Ribo_ID)

                for index,row in mRNA_counts_DF.iterrows():

                    if row["LocusTag"] == jcvi3AID_synth:

                        avg_mRNA_cnt = row["Count"]

                        if avg_mRNA_cnt == 0.0:
                            avg_mRNA_cnt = 0.001

                        init_mRNA_count = np.random.poisson(avg_mRNA_cnt)

                        continue

                sim.distributeNumber(RNA, sim.region('cytoplasm'), init_mRNA_count)
                pmap[rnaMetID] = init_mRNA_count
                sim.distributeNumber(PTN, sim.region('cytoplasm'), int(4*ptnCount*ptn_ratio/5))
                sim.distributeNumber(PTN, sim.region('outer_cytoplasm'), int(ptnCount*ptn_ratio/5))
                pmap[ptnMetID] = int(ptnCount*ptn_ratio)

                diffRNA = sim.diffusionConst('mrna_' + locusNum_synth + '_diff', RNA_diff_coeff(rnasequence))
        #         mrna_diff_list.append(RNA_diff_coeff(rnasequence))  

                RNA.diffusionRate(cyt, diffRNA)
                PTN.diffusionRate(cyt, diffPTN)

                RNA.diffusionRate(mem,sim.diffusionZero)
                RNA.diffusionRate(ext,sim.diffusionZero)
                RNA.diffusionRate(ribo,diffRNA)
                RNA.diffusionRate(dna,diffRNA)
                RNA.diffusionRate(she,diffRNA)

                PTN.diffusionRate(mem,sim.diffusionZero)
                PTN.diffusionRate(ext,sim.diffusionZero)
                PTN.diffusionRate(ribo,diffPTN)
                PTN.diffusionRate(dna,diffPTN)
                PTN.diffusionRate(she,diffPTN)

                sim.transitionRate(RNA, dna, cyt, diffRNA)
                sim.transitionRate(RNA, cyt, dna, diffRNA)
                sim.transitionRate(RNA, cyt, ribo, diffRNA)
                sim.transitionRate(RNA, ribo, cyt, diffRNA)
                sim.transitionRate(RNA, ribo, she, diffRNA)
                sim.transitionRate(RNA, she, ribo, diffRNA)
                sim.transitionRate(RNA, she, cyt, diffRNA)
                sim.transitionRate(RNA, cyt, she, diffRNA)
                sim.transitionRate(RNA, she, dna, diffRNA)
                sim.transitionRate(RNA, dna, she, diffRNA)
                sim.transitionRate(RNA, dna, ribo, diffRNA)
                sim.transitionRate(RNA, ribo, dna, diffRNA)

                sim.transitionRate(PTN, dna, cyt, diffPTN)
                sim.transitionRate(PTN, cyt, dna, diffPTN)
        #         sim.transitionRate(PTN, cyt, ribo, diffPTN)
                sim.transitionRate(PTN, ribo, cyt, diffPTN)
                sim.transitionRate(PTN, ribo, she, diffPTN)
        #         sim.transitionRate(PTN, she, ribo, diffPTN)
                sim.transitionRate(PTN, she, cyt, diffPTN)
                sim.transitionRate(PTN, cyt, she, diffPTN)
                sim.transitionRate(PTN, she, dna, diffPTN)
                sim.transitionRate(PTN, dna, she, diffPTN)
                sim.transitionRate(PTN, dna, ribo, diffPTN)
                sim.transitionRate(PTN, ribo, dna, diffPTN)

        #         trsc_rate = sim.rateConst(locusNum + '_trsc', TranscriptRate(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 1)

                translat_rate = sim.rateConst(locusNum_synth + '_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence), 1)
            
                if max_poly_size > 1:
            
                    aasequence_poly = aasequence[-134:]

                    translat_rate_poly = sim.rateConst(locusNum + '_poly_translat', TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence_poly), 1)

                rna_deg_rate = sim.rateConst(locusNum_synth + '_RNAdeg', DegradationRate(rnaMetID, rnasequence), 1)

                RNAP_on = sim.rateConst(locusNum_synth + 'RNAP_on', RNAP_binding(rnaMetID, ptnMetID, rnasequence, jcvi2ID), 2)

                dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)

                dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)

        #         dna.addReaction([RNAP_gene], [gene, RNApol, RNA], trsc_rate)

        #         dna.addReaction([gene], [gene, RNA], trsc_rate)

                she.addReaction([RNA, degradosome], [Deg_RNA], deg_bind_rate)

                she.addReaction([Deg_RNA], [degradosome], rna_deg_rate)

        #         ribo.addReaction([RNA], [], rna_deg_rate)

#                 ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

        #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

#                 ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)
            
                for mp in range(Max_Poly_Count):
        
                    Ribosome = sim.species('Ribosome'+str(mp+1))

                    for pn in range(mp+1): #range(max_poly_size):

                        if (max_poly_size == 1) or (mp == 0):

                            Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                            Ribo_RNA = sim.species(Ribo_ID)

                            ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

                #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                            ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate)

                            break

                        elif max_poly_size > 1:

                            if pn == 0:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                ribo.addReaction([RNA, Ribosome], [Ribo_RNA], Ribo_on)

                #         ribo.addReaction([Ribo_RNA], [RNA, Ribosome], Ribo_off)

                                Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+2)

                                Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                                ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate)

                            elif (pn+1 == max_poly_size) or (pn+1 == mp+1):

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                ribo.addReaction([Ribo_RNA], [RNA, Ribosome, PTN], translat_rate_poly)

                                break

                            else:

                                Ribo_ID = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+1)

                                Ribo_RNA = sim.species(Ribo_ID)

                                Ribo_ID_Next = 'RB' + str(mp+1) + '_' + locusNum_synth + '_' + str(pn+2)

                                Ribo_RNA_Next = sim.species(Ribo_ID_Next)

                                ribo.addReaction([Ribo_RNA], [Ribo_RNA_Next, PTN], translat_rate_poly)

#                 Trsc_Term = sim.species('Trsc_Term')

#                 dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)


                baseCount = defaultdict(int)
                for base in set(rnasequence):
                    baseCount[base] = rnasequence.count(base)

                # Add total number of monomers to parameter dict

                N_A = baseCount["A"]

                N_U = baseCount["U"]

                N_C = baseCount["C"]

                N_G = baseCount["G"]

                gene_deg_dict = {}

                gene_deg_dict['ATP_mRNAdeg'] = N_A + N_U + N_C + N_G

                gene_deg_dict['AMP_mRNAdeg'] = N_A
                gene_deg_dict['UMP_mRNAdeg'] = N_U
                gene_deg_dict['CMP_mRNAdeg'] = N_C
                gene_deg_dict['GMP_mRNAdeg'] = N_G

                degDict['D_' + locusNum_synth] = gene_deg_dict
                
                
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

                translatATP = 0

                gene_ptn_dict = {}
                gene_ptn_dict['Type'] = 'cytoplasm'
                gene_ptn_dict['States'] = [ptnMetID,glugln_enz,glugln_enz_atp,glugln_enz_atp_aa]

                gene_translat_dict = {}

                for aaCost in AaUsed:

                    aa_ID = aaCost[0]
                    numberUsed = aaCost[1]

                    aaCostID = aaCostMap[aa_ID]

                    gene_translat_dict[aaCostID] = numberUsed

                    translatATP = translatATP + 2*numberUsed


                gene_translat_dict['ATP_translat'] = translatATP

                gene_ptn_dict['Translat_costs'] = gene_translat_dict

                multiStatePtnDict[ptnMetID] = gene_ptn_dict
            
    
    if unchargedMetID in pmap:
        pmap[unchargedMetID] = pmap[unchargedMetID] + int(round(5800/29*0.2))
        print(unchargedMetID,pmap[unchargedMetID])
    if unchargedMetID not in pmap:
        pmap[unchargedMetID] = int(round(5800/29*0.2))
        print(unchargedMetID,pmap[unchargedMetID])
        
    if chargedMetID in pmap:
        pmap[chargedMetID] = pmap[chargedMetID] + int(round(5800/29*0.8))
        print(chargedMetID,pmap[chargedMetID])
    if chargedMetID not in pmap:
        pmap[chargedMetID] = int(round(5800/29*0.8))
        print(chargedMetID,pmap[chargedMetID])
        
    geneMetID = 'g_' + locusNum
    gene = sim.species(geneMetID)
    RNAP_gene = sim.species('RP_' + locusNum)
    pmap['RP_' + locusNum] = 0
    RDME_species_list.append('RP_' + locusNum)
    PartIdxMap['RP_' + locusNum] = partIdx
    partIdx = partIdx + 1
    
    binding_rate = 4*(180*3/612)*Ecoli_V*avgdr/11400/60 #/1800/60
    RNAP_on = sim.rateConst('RNAP_on_tRNA', binding_rate, 2)
    
    dna.addReaction([gene, RNApol], [RNAP_gene], RNAP_on)
    
    dna.addReaction([RNAP_gene], [gene, RNApol], RNAP_off)
    
#     Trsc_Term = sim.species('Trsc_Term')

#     dna.addReaction([RNAP_gene, Trsc_Term], [gene], TrscTermination)
    
#     trna_trsc = sim.rateConst('trsc_'+locusNum,  trnaTranscriptRate(rnasequence), 1)
    
#     dna.addReaction([RNAP_gene], [gene, RNApol, tRNA], trna_trsc)
    
    return partIdx
##################


##################
# Define how to add the particles and reactions for each rRNA and its corresponding gene.

# Define how to add the particles and reactions for each rRNA and its corresponding gene.

def addrRNA(sim,pmap,ext, mem, cyt, ribo, dna, she,rrnaMetDF_1,rrnaMetDF_2,RNApol,RNAP_on,RNAP_off,RDME_species_list):
    
    rrna_gene_locs_1 = []
    rrna_gene_locs_2 = []

    rRNA_species = []
    
    for index, row in rrnaMetDF_1.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
#         rnasequence = getRNAsequences(jcvi3AID)
        gene_str = str(genomeRnaLocDict[jcvi3AID])
        start_nt = int(gene_str.split(':')[0].replace('[',''))
        print(start_nt)
        end_nt = int(gene_str.split(':')[1].replace('](-)',''))
        print(end_nt)
        rrna_gene_locs_1.append(start_nt)
        rrna_gene_locs_1.append(end_nt)
            
        rRNA_species.append(newMetID)
        
    rrna_pos_1 = Seq(str(genome3A.seq[min(rrna_gene_locs_1):max(rrna_gene_locs_1)]))
    
    rrna_operon_1 = rrna_pos_1.reverse_complement()
            
    print(rrna_operon_1)
    
    
    for index, row in rrnaMetDF_2.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
#         rnasequence = getRNAsequences(jcvi3AID)
        gene_str = str(genomeRnaLocDict[jcvi3AID])
        start_nt = int(gene_str.split(':')[0].replace('[',''))
        print(start_nt)
        end_nt = int(gene_str.split(':')[1].replace('](-)',''))
        print(end_nt)
        rrna_gene_locs_2.append(start_nt)
        rrna_gene_locs_2.append(end_nt)
        
    rrna_pos_2 = Seq(str(genome3A.seq[min(rrna_gene_locs_2):max(rrna_gene_locs_2)]))
    
    rrna_operon_2 = rrna_pos_2.reverse_complement()
            
    print(rrna_operon_2)
        
        
    for index, row in rrnaMetDF_1.iterrows():
        newMetID = row["species"]
        jcvi3AID = row["gene"]
        
        rRNA = sim.species(newMetID)
        pmap[newMetID] = 0
        RDME_species_list.append(newMetID)

        print(jcvi3AID)
        locusNum = jcvi3AID.split('3A_')[1]
        
        rnasequence = getRNAsequences(jcvi3AID)
        
        diffRNA = sim.diffusionConst('rrna_' + locusNum + '_diff', RNA_diff_coeff(rnasequence))

        rRNA.diffusionRate(mem,sim.diffusionZero)
        rRNA.diffusionRate(ext,sim.diffusionZero)
        rRNA.diffusionRate(cyt,diffRNA)
        rRNA.diffusionRate(ribo,diffRNA)
        rRNA.diffusionRate(dna,diffRNA)
        rRNA.diffusionRate(she,diffRNA)
        
        sim.transitionRate(rRNA, dna, cyt, diffRNA)
        sim.transitionRate(rRNA, cyt, dna, diffRNA)
#         sim.transitionRate(rRNA, cyt, ribo, diffRNA)
        sim.transitionRate(rRNA, ribo, cyt, diffRNA)
        sim.transitionRate(rRNA, ribo, she, diffRNA)
#         sim.transitionRate(rRNA, she, ribo, diffRNA)
        sim.transitionRate(rRNA, she, cyt, diffRNA)
        sim.transitionRate(rRNA, cyt, she, diffRNA)
        sim.transitionRate(rRNA, she, dna, diffRNA)
        sim.transitionRate(rRNA, dna, she, diffRNA)
        sim.transitionRate(rRNA, dna, ribo, diffRNA)
        sim.transitionRate(rRNA, ribo, dna, diffRNA)
    
    rRNA_5S = sim.species('M_rRNA_5S_c')
    rRNA_16S = sim.species('M_rRNA_16S_c')
    rRNA_23S = sim.species('M_rRNA_23S_c')
    
    rnasequence = rrna_operon_1.transcribe()
    
    gene_16S = sim.species('g_69')
    gene_23S = sim.species('g_68')
    RNAP_16S = sim.species('RP_69')
    RNAP_23S = sim.species('RP_68')
    pmap['RP_69'] = 0
    pmap['RP_68'] = 0
        
    trsc_rate = sim.rateConst('rRNA_1',rrnaTranscriptRate(rnasequence),1)
    
    dna.addReaction([gene_16S, RNApol], [RNAP_16S], RNAP_on)
    dna.addReaction([RNAP_16S], [gene_16S, RNApol], RNAP_off)
    dna.addReaction([RNAP_16S], [rRNA_5S,rRNA_16S,rRNA_23S,RNApol,gene_16S], trsc_rate)
    
    dna.addReaction([gene_23S, RNApol], [RNAP_23S], RNAP_on)
    dna.addReaction([RNAP_23S], [gene_23S, RNApol], RNAP_off)
    dna.addReaction([RNAP_23S], [rRNA_5S,rRNA_16S,rRNA_23S,RNApol,gene_23S], trsc_rate)
    
    gene_16S = sim.species('g_534')
    gene_23S = sim.species('g_533')
    RNAP_16S = sim.species('RP_534')
    RNAP_23S = sim.species('RP_533')
    
    pmap['RP_534'] = 0
    pmap['RP_533'] = 0
    
    rnasequence = rrna_operon_2.transcribe()
        
    trsc_rate = sim.rateConst('rRNA_2',rrnaTranscriptRate(rnasequence),1)
    
    dna.addReaction([gene_16S, RNApol], [RNAP_16S], RNAP_on)
    dna.addReaction([RNAP_16S], [gene_16S, RNApol], RNAP_off)
    dna.addReaction([RNAP_16S], [rRNA_5S,rRNA_16S,rRNA_23S,RNApol,gene_16S], trsc_rate)
    
    dna.addReaction([gene_23S, RNApol], [RNAP_23S], RNAP_on)
    dna.addReaction([RNAP_23S], [gene_23S, RNApol], RNAP_off)
    dna.addReaction([RNAP_23S], [rRNA_5S,rRNA_16S,rRNA_23S,RNApol,gene_23S], trsc_rate)
    
##################


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
degrad_bind_rate = 11*avgdr*Ecoli_V/60/2400 #7800 #1/M/s

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

ribo_init = 40*Ecoli_V*avgdr/60/6800

ribosomeConc = 500*countToMiliMol # mM

# Concentration of charged tRNA
ctRNAconc = 150*countToMiliMol # mM

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
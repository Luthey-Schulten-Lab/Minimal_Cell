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
    
    PtnMetDF = pd.read_csv("../../model_data/protein_metabolites_frac.csv")
    riboPtnMetDF = pd.read_csv("../../model_data/ribo_protein_metabolites.csv")
    memPtnMetDF = pd.read_csv("../../model_data/membrane_protein_metabolites.csv")
    trnaMetDF = pd.read_csv("../../model_data/trna_metabolites_synthase.csv")
    
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

#     addPtn(csim, 'JCVISYN3A_0227') ### FIX ####
#     addPtnCME(csim, 'JCVISYN3A_0227', pmap, genes_in_model)
    
    tRNAadded = []

    for index, row in trnaMetDF.iterrows():
    #     addtRNA(row["uncharged"], row["charged"], row["gene"])
        addtRNACME(csim, row["uncharged"], row["charged"], row["gene"], row["AA"], row["synthase"], tRNAadded, pmap, ctRNAcostMap)
    #     addtRNA(row["species"], row["gene"])

#     addrRNACME(csim)

#     addRepInit()
#     addReplication()

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
            
#             print(key,val)
            
#         if ('UTP_' in key) or ('GTP_' in key) or ('CTP_' in key):
            
#             csim.defineSpecies([key])
#             if val < 0:
#                 print(key,val)
#             csim.addParticles(key,count=int(val))
            
#             print(key,val)
            
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

#     ModelSpecies.append(geneMetID)

    rnasequence = getRNAsequences(jcvi3AID)
#     print(rnasequence)

    gene = [geneMetID]
#     print(gene)
#     sim.defineSpecies(gene)

#     sim.addParticles(species = geneMetID, count = 1)

    new_RNA_ID = 'New_RNA_' + locusNum
    TrscProd = [new_RNA_ID]
#         mRNAdegProd = []

    for i in range(len(rnasequence)):
        TrscProd.append('ATP_trsc')
#             mRNAdegProd.append('ATP_mRNAdeg')

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
#             mRNAdegProd.append('AMP_mRNAdeg')

    for i in range(N_U):
        TrscProd.append('UTP_tRNA')
#             mRNAdegProd.append('UMP_mRNAdeg')

    for i in range(N_G):
        TrscProd.append('GTP_tRNA')
#             mRNAdegProd.append('GMP_mRNAdeg')

    for i in range(N_C):
        TrscProd.append('CTP_tRNA')
#             mRNAdegProd.append('CMP_mRNAdeg')

    RNAP_gene = 'RP_' + locusNum

    csim.addReaction(RNAP_gene, tuple(TrscProd), trnaTranscriptRateCME(rnasequence, pmap))

    if unchargedMetID not in tRNAadded:
#         ModelSpecies.append(unchargedMetID)
        tRNA = [unchargedMetID]
#         sim.defineSpecies(tRNA)
#         sim.addParticles(species = newMetID,  count = int(round(7000/20*0.2)))
        tRNAadded.append(unchargedMetID)
#         ModelSpecies.append(chargedMetID)
        ctRNA = [chargedMetID]
#         sim.defineSpecies(ctRNA)
#         sim.addParticles(species = newMetID,  count = int(round(7000/20*0.2)))
        tRNAadded.append(chargedMetID)
    
        if 'GLN' not in synthase:
            
            synthase = synthase.split('3A_')[1].lstrip('0')
            synthase_ptn = 'P_' + synthase
    
            synthase_atp = synthase_ptn + '_ATP'
#             if synthase_atp not in pmap:
#                 pmap[synthase_atp] = 1
#                 synth_atp = [synthase_atp]
#                 sim.defineSpecies(synth_atp)
#                 sim.addParticles(species = synthase_atp, count = 1)
            
            synthase_atp_aa = synthase_atp + '_AA'
#             if synthase_atp_aa not in pmap:
#                 pmap[synthase_atp_aa] = 1
#                 synth_atp_aa = [synthase_atp_aa]
#                 sim.defineSpecies(synth_atp_aa)
#                 sim.addParticles(species = synthase_atp_aa, count = 1)
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNA'
#             if synthase_atp_aa_trna not in pmap:
#                 pmap[synthase_atp_aa_trna] = 1
#                 synth_atp_aa_trna = [synthase_atp_aa_trna]
#                 sim.defineSpecies(synth_atp_aa_trna)
#                 sim.addParticles(species = synthase_atp_aa_trna, count = 1)
            

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
#             synthase_ptn = 'M_PTN_' + synthase + '_c'
            
            synthase_atp_aa_trna = synthase_atp_aa + '_tRNAgln'
            synth_atp_aa_trna = [synthase_atp_aa_trna]
#             sim.defineSpecies(synth_atp_aa_trna)
#             sim.addParticles(species = synthase_atp_aa_trna, count = 1)
            
            csim.addReaction((synthase_atp_aa,'M_trnagln_c'),(synthase_atp_aa_trna),0.1)
            csim.addReaction((synthase_atp_aa_trna),('M_glutrnagln_c','M_amp_c','M_ppi_c',synthase_ptn),30)
            
            glugln = 'M_glutrnagln_c'
            glugln_enz = 'glutrnagln_enz'
            glugln_enz_atp = 'glutrnagln_enz_atp'
            glugln_enz_atp_aa = 'glutrnagln_enz_atp_gln'
            glugln_enz_atp_asn = 'glutrnagln_enz_atp_asn'
#             gge = [glugln_enz]
#             sim.defineSpecies(gge)
#             sim.addParticles(species = glugln_enz, count = 1)
#             ggea = [glugln_enz_atp]
#             sim.defineSpecies(ggea)
#             sim.addParticles(species = glugln_enz_atp, count = 1)
#             ggeaaa = [glugln_enz_atp_aa]
#             sim.defineSpecies(ggeaaa)
#             sim.addParticles(species = glugln_enz_atp_aa, count = 1)
            
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
            
            print('GLN is dumb')
#             print('Added charging reactions for '+unchargedMetID)
    
#     sim.addParticles(species = unchargedMetID,  count = int(round(5800/29*0.2)))
#     sim.addParticles(species = chargedMetID,  count = int(round(5800/29*0.8)))
#     TrscProd = [geneMetID,unchargedMetID]
        
#     for i in range(len(rnasequence)):
#         TrscProd.append('ATP_trsc')
        
#     baseCount = defaultdict(int)
#     for base in set(rnasequence):
#         baseCount[base] = rnasequence.count(base)

#     # Add total number of monomers to parameter dict

#     N_A = baseCount["A"]

#     N_U = baseCount["U"]

#     N_C = baseCount["C"]

#     N_G = baseCount["G"]

#     for i in range(N_A):
#         TrscProd.append('ATP_tRNA')

#     for i in range(N_U):
#         TrscProd.append('UTP_tRNA')

#     for i in range(N_G):
#         TrscProd.append('GTP_tRNA')

#     for i in range(N_C):
#         TrscProd.append('CTP_tRNA')
        
    
#     sim.addReaction(geneMetID, tuple(TrscProd), trnaTranscriptRate(unchargedMetID, rnasequence))


# def addrRNACME(csim):
    
#     rrna_gene_locs_1 = []
#     rrna_gene_locs_2 = []

#     rRNA_species = []
    
#     for index, row in rrnaMetDF_1.iterrows():
#         newMetID = row["species"]
#         jcvi3AID = row["gene"]
        
# #         rnasequence = getRNAsequences(jcvi3AID)
#         gene_str = str(genomeRnaLocDict[jcvi3AID])
#         start_nt = int(gene_str.split(':')[0].replace('[',''))
#         print(start_nt)
#         end_nt = int(gene_str.split(':')[1].replace('](-)',''))
#         print(end_nt)
#         rrna_gene_locs_1.append(start_nt)
#         rrna_gene_locs_1.append(end_nt)
            
#         rRNA_species.append(newMetID)
        
#     rrna_pos_1 = Seq(str(genome3A.seq[min(rrna_gene_locs_1):max(rrna_gene_locs_1)]))
    
#     rrna_operon_1 = rrna_pos_1.reverse_complement()
            
#     print(rrna_operon_1)
    
    
#     for index, row in rrnaMetDF_2.iterrows():
#         newMetID = row["species"]
#         jcvi3AID = row["gene"]
        
# #         rnasequence = getRNAsequences(jcvi3AID)
#         gene_str = str(genomeRnaLocDict[jcvi3AID])
#         start_nt = int(gene_str.split(':')[0].replace('[',''))
#         print(start_nt)
#         end_nt = int(gene_str.split(':')[1].replace('](-)',''))
#         print(end_nt)
#         rrna_gene_locs_2.append(start_nt)
#         rrna_gene_locs_2.append(end_nt)
        
#     rrna_pos_2 = Seq(str(genome3A.seq[min(rrna_gene_locs_2):max(rrna_gene_locs_2)]))
    
#     rrna_operon_2 = rrna_pos_2.reverse_complement()
            
#     print(rrna_operon_2)
    
#     print(rRNA_species)
    
#     csim.defineSpecies(rRNA_species)
    
    
#     for rRNA in rRNA_species:
        
#         csim.addParticles(species = rRNA,  count = 1)
            
#         ModelSpecies.append(rRNA)
    
#     for index, row in rrnaMetDF_1.iterrows():
#         newMetID = row["species"]
#         jcvi3AID = row["gene"]
        
#         geneMetID = jcvi3AID + '_gene'
    
#         ModelSpecies.append(geneMetID)
        
#         gene = [geneMetID]
#         print(gene)
#         csim.defineSpecies(gene)

#         csim.addParticles(species = geneMetID, count = 1)

#         print(jcvi3AID)

#     rnasequence = rrna_operon_1.transcribe()
        
#     TrscProd = ['JCVISYN3A_0067_gene','M_rRNA_5S_c','M_rRNA_23S_c','M_rRNA_16S_c']
        
#     for i in range(int(len(rnasequence))):
#         TrscProd.append('ATP_trsc')
        
#     baseCount = defaultdict(int)
#     for base in set(rnasequence):
#         baseCount[base] = rnasequence.count(base)

#     # Add total number of monomers to parameter dict

#     N_A = baseCount["A"]

#     N_U = baseCount["U"]

#     N_C = baseCount["C"]

#     N_G = baseCount["G"]

#     for i in range(N_A):
#         TrscProd.append('ATP_rRNA')

#     for i in range(N_U):
#         TrscProd.append('UTP_rRNA')

#     for i in range(N_G):
#         TrscProd.append('GTP_rRNA')

#     for i in range(N_C):
#         TrscProd.append('CTP_rRNA')

#     csim.addReaction('JCVISYN3A_0067_gene', tuple(TrscProd), rrnaTranscriptRate(rnasequence))
        
        
#     for index, row in rrnaMetDF_2.iterrows():
#         newMetID = row["species"]
#         jcvi3AID = row["gene"]
        
#         geneMetID = jcvi3AID + '_gene'
    
#         ModelSpecies.append(geneMetID)
        
#         gene = [geneMetID]
#         print(gene)
#         csim.defineSpecies(gene)

#         csim.addParticles(species = geneMetID, count = 1)

#         print(jcvi3AID)

#     rnasequence = rrna_operon_2.transcribe()
        
#     TrscProd = ['JCVISYN3A_0532_gene','M_rRNA_5S_c','M_rRNA_23S_c','M_rRNA_16S_c']
        
#     for i in range(int(len(rnasequence))):
#         TrscProd.append('ATP_trsc')
        
#     baseCount = defaultdict(int)
#     for base in set(rnasequence):
#         baseCount[base] = rnasequence.count(base)

#     # Add total number of monomers to parameter dict

#     N_A = baseCount["A"]

#     N_U = baseCount["U"]

#     N_C = baseCount["C"]

#     N_G = baseCount["G"]

#     for i in range(N_A):
#         TrscProd.append('ATP_rRNA')

#     for i in range(N_U):
#         TrscProd.append('UTP_rRNA')

#     for i in range(N_G):
#         TrscProd.append('GTP_rRNA')

#     for i in range(N_C):
#         TrscProd.append('CTP_rRNA')

#     csim.addReaction('JCVISYN3A_0532_gene', tuple(TrscProd), rrnaTranscriptRate(rnasequence))
    
    

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
    
#     ptnCount = ptnFrac*ptnCount
    
#     ModelSpecies.append(newMetID)

#     print(newMetID, ptnCount)
    
    if transcribe:

        geneMetID = jcvi3AID + '_gene'

#         ModelSpecies.append(geneMetID)

        rnasequence, aasequence = getSequences(jcvi3AID)
#         print(rnasequence)
#         print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName

#         ModelSpecies.append(rnaMetID)
        
#         cytoPtnMetID = jcvi3AID + '_cyto'
#         ModelSpecies.append(cytoPtnMetID)

#         species = []
#         species = [geneMetID, rnaMetID, newMetID, cytoPtnMetID]
#         print(species)
#         csim.defineSpecies(species)
        
#         for index,row in mRNA_counts_DF.iterrows():

#             if row["LocusTag"] == jcvi3AID:

#                 avg_mRNA_cnt = row["Count"]

#                 if avg_mRNA_cnt == 0.0:
#                     avg_mRNA_cnt = 0.001

#                 init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
#                 continue

#         csim.addParticles(species = geneMetID, count = 1)
#         csim.addParticles(species = newMetID,  count = int(ptnCount))
#         csim.addParticles(species = rnaMetID, count = init_mRNA_count)

        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]
#         mRNAdegProd = []

        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
#             mRNAdegProd.append('ATP_mRNAdeg')

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
#             mRNAdegProd.append('AMP_mRNAdeg')

        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
#             mRNAdegProd.append('UMP_mRNAdeg')

        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
#             mRNAdegProd.append('GMP_mRNAdeg')

        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
#             mRNAdegProd.append('CMP_mRNAdeg')

#         TranslatProd = [rnaMetID,cytoPtnMetID]
#         ptnDegProd = []

#         for i in range(len(aasequence)):
#             TranslatProd.append('ATP_translat')
#             TranslatProd.append('ATP_translat')
#     #         ptnDegProd.append('ATP_ptndeg')
            
#         aaCount = defaultdict(int)
#         for aa in set(aasequence):
#             aaCount[aa] = aasequence.count(aa)

#         NMono_A = aaCount["A"]
#         NMono_R = aaCount["R"]
#         NMono_N = aaCount["N"]
#         NMono_D = aaCount["D"]
#         NMono_C = aaCount["C"]
#         NMono_E = aaCount["E"]
#         NMono_Q = aaCount["Q"]
#         NMono_H = aaCount["H"]
#         NMono_I = aaCount["I"]
#         NMono_L = aaCount["L"]
#         NMono_K = aaCount["K"]
#         NMono_M = max(0,aaCount["M"] - 1)
#         NMono_P = aaCount["P"]
#         NMono_S = aaCount["S"]
#         NMono_T = aaCount["T"]
#         NMono_W = aaCount["W"]
#         NMono_Y = aaCount["Y"]
#         NMono_G = aaCount["G"]
#         NMono_F = aaCount["F"]
#         NMono_V = aaCount["V"]
#         NMono_FM = 1

#         NStop = aaCount["*"]

#         if NStop > 1:
#             print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

#         AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
#                      ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
#                      ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
#                      ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
#                      ['FM',NMono_FM]]

#         for aaCost in AaUsed:
#             aa_ID = aaCost[0]
#             numberUsed = aaCost[1]

#             aaCostID = aaCostMap[aa_ID]
# #             print(aa_ID,aaCostID)

#             for i in range(numberUsed):
#                 TranslatProd.append(aaCostID)

        RNAP_gene = 'RP_' + locusNum

        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(newMetID, rnasequence, jcvi2ID, pmap))
#         csim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
#         csim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))


#         ptnLen = sum(list(aaCount.values()))

#         translocProd = [newMetID]

#         translocATP = int(ptnLen/10)
        
#         secY = 'M_PTN_JCVISYN3A_0652_c'

#         for i in range(translocATP):
#             translocProd.append('ATP_transloc')

#         sim.addReaction((cytoPtnMetID),tuple(translocProd),TranslocRate(aasequence)) #secY,
    
#     if not transcribe:
#         continue
#         species = [newMetID]
#         sim.defineSpecies(species)
#         sim.addParticles(species = newMetID,  count = int(ptnCount))
    

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

#     print(newMetID, ptnCount)
    
#     ModelSpecies.append(newMetID)

#     geneMetID = jcvi3AID + '_gene'

#     ModelSpecies.append(geneMetID)
    
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):

        rnaMetID = "M_RNA_" + jcvi3AID + "_c"
        rnaName = "(mRNA) " + ptnName
        
#     ModelSpecies.append(rnaMetID)

#     species = []
#     species = [geneMetID, rnaMetID, newMetID]
#     print(species)
#     sim.defineSpecies(species)
    
#     for index,row in mRNA_counts_DF.iterrows():
        
#         if row["LocusTag"] == jcvi3AID:
            
#             avg_mRNA_cnt = row["Count"]
            
#             if avg_mRNA_cnt == 0.0:
#                 avg_mRNA_cnt = 0.001
            
#             init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
            
#             continue

#     sim.addParticles(species = geneMetID, count = 1)
#     sim.addParticles(species = newMetID,  count = int(ptnCount))
#     sim.addParticles(species = rnaMetID, count = init_mRNA_count)
    
#     RiboPtnNames.append(newMetID)
    
    trsc_rate = riboTranscriptRate(rnaMetID, newMetID, rnasequence, jcvi2ID)
#     RiboPtnTrscRates.append(trsc_rate)
    
#     translat_rate = TranslatRate(rnaMetID, newMetID, rnasequence, aasequence)
#     RiboPtnTranslatRates.append(translat_rate)
    
#     rna_length = len(rnasequence)
#     RiboPtnLens.append(rna_length)

    new_RNA_ID = 'New_mRNA_' + locusNum

    TrscProd = [new_RNA_ID]
#     mRNAdegProd = []
    
    for i in range(len(rnasequence)):
        TrscProd.append('ATP_trsc')
#         mRNAdegProd.append('ATP_mRNAdeg')
        
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
#         mRNAdegProd.append('AMP_mRNAdeg')

    for i in range(N_U):
        TrscProd.append('UTP_mRNA')
#         mRNAdegProd.append('UMP_mRNAdeg')

    for i in range(N_G):
        TrscProd.append('GTP_mRNA')
#         mRNAdegProd.append('GMP_mRNAdeg')

    for i in range(N_C):
        TrscProd.append('CTP_mRNA')
#         mRNAdegProd.append('CMP_mRNAdeg')
        
#     TranslatProd = [rnaMetID,newMetID]
        
#     for i in range(len(aasequence)):
#         TranslatProd.append('ATP_translat')
#         TranslatProd.append('ATP_translat')
            
#     aaCount = defaultdict(int)
#     for aa in set(aasequence):
#         aaCount[aa] = aasequence.count(aa)

#     NMono_A = aaCount["A"]
#     NMono_R = aaCount["R"]
#     NMono_N = aaCount["N"]
#     NMono_D = aaCount["D"]
#     NMono_C = aaCount["C"]
#     NMono_E = aaCount["E"]
#     NMono_Q = aaCount["Q"]
#     NMono_H = aaCount["H"]
#     NMono_I = aaCount["I"]
#     NMono_L = aaCount["L"]
#     NMono_K = aaCount["K"]
#     NMono_M = max(0,aaCount["M"] - 1)
#     NMono_P = aaCount["P"]
#     NMono_S = aaCount["S"]
#     NMono_T = aaCount["T"]
#     NMono_W = aaCount["W"]
#     NMono_Y = aaCount["Y"]
#     NMono_G = aaCount["G"]
#     NMono_F = aaCount["F"]
#     NMono_V = aaCount["V"]
#     NMono_FM = 1

#     NStop = aaCount["*"]

#     if NStop > 1:
#         print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

#     AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
#                  ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
#                  ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
#                  ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
#                  ['FM',NMono_FM]]

#     for aaCost in AaUsed:
#         aa_ID = aaCost[0]
#         numberUsed = aaCost[1]

#         aaCostID = aaCostMap[aa_ID]

#         for i in range(numberUsed):
#             TranslatProd.append(aaCostID)

    RNAP_gene = 'RP_' + locusNum
    
    csim.addReaction(RNAP_gene, tuple(TrscProd), riboTranscriptRateCME( newMetID, rnasequence, jcvi2ID, pmap))
#     csim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
#     csim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))

        

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
    
#     ptnCount = ptnFrac*ptnCount
    
#     ModelSpecies.append(newMetID)

#     print(newMetID, ptnCount)
    
    if transcribe:

#         geneMetID = jcvi3AID + '_gene'

#         ModelSpecies.append(geneMetID)

        rnasequence, aasequence = getSequences(jcvi3AID)
#         print(rnasequence)
#         print(aasequence)

        if (rnasequence != 0) and (aasequence != 0):

            rnaMetID = "M_RNA_" + jcvi3AID + "_c"
            rnaName = "(mRNA) " + ptnName

#         ModelSpecies.append(rnaMetID)

#         species = []
#         species = [geneMetID, rnaMetID, newMetID]
#         print(species)
#         sim.defineSpecies(species)
        
#         for index,row in mRNA_counts_DF.iterrows():

#             if row["LocusTag"] == jcvi3AID:

#                 avg_mRNA_cnt = row["Count"]

#                 if avg_mRNA_cnt == 0.0:
#                     avg_mRNA_cnt = 0.001

#                 init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
#                 continue


#         csim.addParticles(species = geneMetID, count = 1)
#         csim.addParticles(species = newMetID,  count = int(ptnCount))
#         csim.addParticles(species = rnaMetID, count = init_mRNA_count)

        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]
#         mRNAdegProd = []

        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
#             mRNAdegProd.append('ATP_mRNAdeg')

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
#             mRNAdegProd.append('AMP_mRNAdeg')

        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
#             mRNAdegProd.append('UMP_mRNAdeg')

        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
#             mRNAdegProd.append('GMP_mRNAdeg')

        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
#             mRNAdegProd.append('CMP_mRNAdeg')

#         TranslatProd = [rnaMetID,newMetID]
#         ptnDegProd = []

#         for i in range(len(aasequence)):
#             TranslatProd.append('ATP_translat')
#             TranslatProd.append('ATP_translat')
#     #         ptnDegProd.append('ATP_ptndeg')
            
#         aaCount = defaultdict(int)
#         for aa in set(aasequence):
#             aaCount[aa] = aasequence.count(aa)

#         NMono_A = aaCount["A"]
#         NMono_R = aaCount["R"]
#         NMono_N = aaCount["N"]
#         NMono_D = aaCount["D"]
#         NMono_C = aaCount["C"]
#         NMono_E = aaCount["E"]
#         NMono_Q = aaCount["Q"]
#         NMono_H = aaCount["H"]
#         NMono_I = aaCount["I"]
#         NMono_L = aaCount["L"]
#         NMono_K = aaCount["K"]
#         NMono_M = max(0,aaCount["M"] - 1)
#         NMono_P = aaCount["P"]
#         NMono_S = aaCount["S"]
#         NMono_T = aaCount["T"]
#         NMono_W = aaCount["W"]
#         NMono_Y = aaCount["Y"]
#         NMono_G = aaCount["G"]
#         NMono_F = aaCount["F"]
#         NMono_V = aaCount["V"]
#         NMono_FM = 1

#         NStop = aaCount["*"]

#         if NStop > 1:
#             print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

#         AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
#                      ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
#                      ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
#                      ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
#                      ['FM',NMono_FM]]

#         for aaCost in AaUsed:
#             aa_ID = aaCost[0]
#             numberUsed = aaCost[1]

#             aaCostID = aaCostMap[aa_ID]
# #             print(aa_ID,aaCostID)

#             for i in range(numberUsed):
#                 TranslatProd.append(aaCostID)

        RNAP_gene = 'RP_' + locusNum

        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(newMetID, rnasequence, jcvi2ID, pmap))
#         csim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
#         csim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, newMetID, rnasequence, aasequence))

#         sim.addReaction(newMetID,tuple(ptnDegProd),ptnDegRate)
    
#     if not transcribe:
#         continue
#         species = [newMetID]
#         sim.defineSpecies(species)
#         sim.addParticles(species = newMetID,  count = int(ptnCount))

        

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

#     ptnMetID = 'M_PTN_' + jcvi3AID + '_c'
    ptnMetID = 'P_' + locusNum
    # If the protein is not in the model, add it:

#     ModelSpecies.append(ptnMetID)

    ptnCount, ptnName = getPtnCount(ptnMetID, jcvi2ID)

#     print(ptnMetID, ptnCount)

#     geneMetID = jcvi3AID + '_gene'
    
#     ModelSpecies.append(geneMetID)

    # Get nucleotide and amino acid sequences, if available
    rnasequence, aasequence = getSequences(jcvi3AID)
#     print(rnasequence)
#     print(aasequence)

    if (rnasequence != 0) and (aasequence != 0):

#         rnaMetID = "M_RNA_" + jcvi3AID + "_c"
#         rnaMetID = "RNA_" + locusNum
#         rnaName = "(mRNA) " + ptnName

#         ModelSpecies.append(rnaMetID)
    
#         species = []
#         species = [geneMetID, rnaMetID, ptnMetID]
#         print(species)
#         csim.defineSpecies(species)
        
#         for index,row in mRNA_counts_DF.iterrows():
        
#             if row["LocusTag"] == jcvi3AID:

#                 avg_mRNA_cnt = row["Count"]

#                 if avg_mRNA_cnt == 0.0:
#                     avg_mRNA_cnt = 0.001

#                 init_mRNA_count = np.random.poisson(avg_mRNA_cnt)
                
#                 continue

#         csim.addParticles(species = geneMetID, count = 1)
#         csim.addParticles(species = ptnMetID,  count = int(ptnCount))
#         csim.addParticles(species = rnaMetID, count = init_mRNA_count)

#         csim.defineSpecies(['New_mRNA_' + locusNum])
        
        new_RNA_ID = 'New_mRNA_' + locusNum

        TrscProd = [new_RNA_ID]
#         mRNAdegProd = []
    
        for i in range(len(rnasequence)):
            TrscProd.append('ATP_trsc')
#             mRNAdegProd.append('ATP_mRNAdeg')
            
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
#             mRNAdegProd.append('AMP_mRNAdeg')
            
        for i in range(N_U):
            TrscProd.append('UTP_mRNA')
#             mRNAdegProd.append('UMP_mRNAdeg')
            
        for i in range(N_G):
            TrscProd.append('GTP_mRNA')
#             mRNAdegProd.append('GMP_mRNAdeg')
            
        for i in range(N_C):
            TrscProd.append('CTP_mRNA')
#             mRNAdegProd.append('CMP_mRNAdeg')
            
#         TranslatProd = [rnaMetID,ptnMetID]
#         ptnDegProd = []
        
#         for i in range(len(aasequence)):
#             TranslatProd.append('ATP_translat')
#             TranslatProd.append('ATP_translat')
# #             ptnDegProd.append('ATP_ptndeg')
            
#         aaCount = defaultdict(int)
#         for aa in set(aasequence):
#             aaCount[aa] = aasequence.count(aa)

#         NMono_A = aaCount["A"]
#         NMono_R = aaCount["R"]
#         NMono_N = aaCount["N"]
#         NMono_D = aaCount["D"]
#         NMono_C = aaCount["C"]
#         NMono_E = aaCount["E"]
#         NMono_Q = aaCount["Q"]
#         NMono_H = aaCount["H"]
#         NMono_I = aaCount["I"]
#         NMono_L = aaCount["L"]
#         NMono_K = aaCount["K"]
#         NMono_M = max(0,aaCount["M"] - 1)
#         NMono_P = aaCount["P"]
#         NMono_S = aaCount["S"]
#         NMono_T = aaCount["T"]
#         NMono_W = aaCount["W"]
#         NMono_Y = aaCount["Y"]
#         NMono_G = aaCount["G"]
#         NMono_F = aaCount["F"]
#         NMono_V = aaCount["V"]
#         NMono_FM = 1

#         NStop = aaCount["*"]

#         if NStop > 1:
#             print("EXTRA STOP CODON: MISTAKE IN TRANSLATION")

#         AaUsed = [['A',NMono_A],['R',NMono_R],['N',NMono_N],['D',NMono_D],['C',NMono_C],
#                      ['E',NMono_E],['Q',NMono_Q],['H',NMono_H],['I',NMono_I],['L',NMono_L],
#                      ['K',NMono_K],['M',NMono_M],['P',NMono_P],['S',NMono_S],['T',NMono_T],
#                      ['W',NMono_W],['Y',NMono_Y],['G',NMono_G],['F',NMono_F],['V',NMono_V],
#                      ['FM',NMono_FM]]
        
#         for aaCost in AaUsed:
#             aa_ID = aaCost[0]
#             numberUsed = aaCost[1]
            
#             aaCostID = aaCostMap[aa_ID]
            
#             for i in range(numberUsed):
#                 TranslatProd.append(aaCostID)

        RNAP_gene = 'RP_' + locusNum
    
        csim.addReaction(RNAP_gene, tuple(TrscProd), TranscriptRateCME(ptnMetID, rnasequence, jcvi2ID, pmap))
#         csim.addReaction(rnaMetID, tuple(mRNAdegProd), DegradationRate(rnaMetID, rnasequence))
#         csim.addReaction(rnaMetID, tuple(TranslatProd), TranslatRate(rnaMetID, ptnMetID, rnasequence, aasequence))
        
#         sim.addReaction(ptnMetID,tuple(ptnDegProd),ptnDegRate)


# # REPLICATION INITIATION WRITTEN BY COLE CROTTY AND MODIFIED BY ZANE THORNBURG
# NA = 6.022*(10**(23))

# def addRepInit(csim,pmap):
#     k_high = 7800*1000/NA/cellVolume
#     k_low = 35*1000/NA/cellVolume
#     k_on = 100*1000/NA/cellVolume #molecule^-1 sec^-1
#     k_off = 0.55 #sec^-1
    
#     helicase_removal_rate = 600 #s^-1
    
#     nonOricSpec = ['High_Affinity_Site','High_Affinity_Bound','High_Affinity_Site_oriC','High_Affinity_Bound_oriC','Low_Affinity_Site_1','Low_Affinity_Site_2','Low_Affinity_Bound_1','Low_Affinity_Bound_2']
#     csim.defineSpecies(nonOricSpec)
#     csim.addParticles(species = 'High_Affinity_Site', count = 16)
#     csim.addParticles(species = 'High_Affinity_Bound', count = 0)
#     csim.addReaction(('High_Affinity_Site', 'M_DnaA_c'),'High_Affinity_Bound',k_high)
    
#     csim.addParticles(species = 'High_Affinity_Site_oriC', count = 1)
#     csim.addParticles(species = 'Low_Affinity_Site_1', count = 0)
#     csim.addParticles(species = 'Low_Affinity_Site_2', count = 0)
    
#     csim.addReaction(('High_Affinity_Site_oriC', 'M_DnaA_c'),('High_Affinity_Bound_oriC','Low_Affinity_Site_1'),k_high)
#     csim.addReaction(('Low_Affinity_Site_1', 'M_DnaA_c'),('Low_Affinity_Bound_1','Low_Affinity_Site_2'),k_low)
#     csim.addReaction(('Low_Affinity_Site_2', 'M_DnaA_c'),('Low_Affinity_Bound_2','ssDNAunboundSite_1'),k_low)

#     species = []

#     for i in range(30):         #loop adds 30 terms for unbound sites
#         term = 'ssDNAunboundSite_'
#         term = term + str(i+1)
#         species.append(term)
#     for j in range(30):         #adds 30 more terms of Bound sites
#         bnd = 'ssDNABoundSite_'
#         bnd = bnd + str(j+1)
#         species.append(bnd)
#     for k in range(30):
#         unbnd = 'ssDNA_Unbinding_'
#         unbnd = unbnd + str(k+1)
#         species.append(unbnd)
        
#     species.append('Initiator_C')
#     species.append('Initiator_CC')

#     csim.defineSpecies(species)

#     for specy in species:
#         csim.addParticles(species = specy, count = 0) #gives all particles count of 0 

#     csim.addParticles(species = 'ssDNAunboundSite_1', count = 0) #starts at 1(only site existing @ start
    
#     # Bidning and unbinding reactions are constructed so that DnaA in the middle of the filament cannot unbind and
#     # DnaA only bind next to the last DnaA in the filament.

#     # Add the binding reactions for each of the 30 binding events in filament formation.
#     for i in range (1,30):
#         csim.addReaction(('ssDNAunboundSite_' + str(i),'M_DnaA_c'), ('ssDNAunboundSite_' + str(i+1),'ssDNABoundSite_' + str(i)),k_on)
#         csim.addReaction(('ssDNAunboundSite_' + str(i+1), 'ssDNABoundSite_' + str(i)),('ssDNAunboundSite_' + str(i),'M_DnaA_c'),k_off)

#     csim.addReaction(('ssDNAunboundSite_30', 'M_DnaA_c'),('ssDNABoundSite_30','ssDNA_Unbinding_30','Initiator_C','Initiator_CC'),k_on)
    
#     # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
#     for i in range (2,31):
#         csim.addReaction(('ssDNA_Unbinding_' + str(i), 'ssDNABoundSite_' + str(i)),('ssDNA_Unbinding_' + str(i-1),'M_DnaA_c'),helicase_removal_rate)
        
#     csim.addReaction(('ssDNA_Unbinding_1', 'ssDNABoundSite_1'),'M_DnaA_c',helicase_removal_rate)



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

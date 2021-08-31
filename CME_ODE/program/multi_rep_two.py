import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, OrderedDict

NA = 6.022*(10**(23))
r_cell = 2.0*(10**-7) # m

CytoVolume = (4*np.pi/3)*1000*r_cell**3 # L
cellVolume = CytoVolume

countToMiliMol = 1000/(NA*cellVolume)

def getConc(particles,specDict):
    """
    Convert particle counts to mM concentrations for the ODE Solver
        
    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's

    cellVolume = specDict['CellV']
    cellVolume = cellVolume * 1e-19

    conc = (particles*1000.0)/(NA*cellVolume)

    return conc

def addRepInitTwo(sim,specDict):
    k_high = 7800*1000/NA/cellVolume
    k_low = 35*1000/NA/cellVolume
    k_on = 100*1000/NA/cellVolume #molecule^-1 sec^-1
    k_on_2 = k_on/2
    k_off = 0.55 #sec^-1
    
    helicase_removal_rate = 600 #s^-1
    
    #Chromosome 1
    nonOricSpec = ['High_Affinity_Site','High_Affinity_Bound','High_Affinity_Site_oriC',
                   'High_Affinity_Bound_oriC','Low_Affinity_Site_1','Low_Affinity_Site_2',
                   'Low_Affinity_Bound_1','Low_Affinity_Bound_2']
    
    sim.addReaction(('High_Affinity_Site', 'M_DnaA_c'),'High_Affinity_Bound',k_high)
    
    sim.addReaction(('High_Affinity_Site_oriC', 'M_DnaA_c'),('High_Affinity_Bound_oriC','Low_Affinity_Site_1'),k_high)
    sim.addReaction(('Low_Affinity_Site_1', 'M_DnaA_c'),('Low_Affinity_Bound_1','Low_Affinity_Site_2'),k_low)
    sim.addReaction(('Low_Affinity_Site_2', 'M_DnaA_c'),('Low_Affinity_Bound_2','ssDNAunboundSite_1'),k_low)
    
    species = []

    for i in range(30):         #loop adds 30 terms for unbound sites
        term = 'ssDNAunboundSite_'
        term = term + str(i+1)
        species.append(term)
    for j in range(30):         #adds 30 more terms of Bound sites
        bnd = 'ssDNABoundSite_'
        bnd = bnd + str(j+1)
        species.append(bnd)
    for k in range(30):
        unbnd = 'ssDNA_Unbinding_'
        unbnd = unbnd + str(k+1)
        species.append(unbnd)

    species.append('Initiator_C')
    species.append('Initiator_CC')
    
    for i in range (1,30):
        sim.addReaction(('ssDNAunboundSite_' + str(i),'M_DnaA_c'), ('ssDNAunboundSite_' + str(i+1),'ssDNABoundSite_' + str(i)),k_on)
        sim.addReaction(('ssDNAunboundSite_' + str(i+1), 'ssDNABoundSite_' + str(i)),('ssDNAunboundSite_' + str(i),'M_DnaA_c'),k_off)


    sim.addReaction(('ssDNAunboundSite_30', 'M_DnaA_c'),('ssDNABoundSite_30','ssDNA_Unbinding_30','Initiator_C','Initiator_CC'),k_on)

    
    #Chromosme 2
    nonOricSpec = ['High_Affinity_Site2','High_Affinity_Bound2','High_Affinity_Site_oriC2',
                   'High_Affinity_Bound_oriC2','Low_Affinity_Site2_1','Low_Affinity_Site2_2',
                   'Low_Affinity_Bound2_1','Low_Affinity_Bound2_2']
   
    sim.addReaction(('High_Affinity_Site2', 'M_DnaA_c'),'High_Affinity_Bound2',k_high)
    
    sim.addReaction(('High_Affinity_Site_oriC2', 'M_DnaA_c'),('High_Affinity_Bound_oriC2','Low_Affinity_Site2_1'),k_high)
    sim.addReaction(('Low_Affinity_Site2_1', 'M_DnaA_c'),('Low_Affinity_Bound2_1','Low_Affinity_Site2_2'),k_low)
    sim.addReaction(('Low_Affinity_Site2_2', 'M_DnaA_c'),('Low_Affinity_Bound2_2','ssDNAunboundSite2_1'),k_low)

    species = []

    for i in range(30):         #loop adds 30 terms for unbound sites
        term = 'ssDNAunboundSite2_'
        term = term + str(i+1)
        species.append(term)
    for j in range(30):         #adds 30 more terms of Bound sites
        bnd = 'ssDNABoundSite2_'
        bnd = bnd + str(j+1)
        species.append(bnd)
    for k in range(30):
        unbnd = 'ssDNA_Unbinding2_'
        unbnd = unbnd + str(k+1)
        species.append(unbnd)

    species.append('Initiator_C2')
    species.append('Initiator_CC2')

    # Binding and unbinding reactions are constructed so that DnaA in the middle of the filament cannot unbind and
    # DnaA only bind next to the last DnaA in the filament.

    # Add the binding reactions for each of the 30 binding events in filament formation.
    for i in range (1,30):
        sim.addReaction(('ssDNAunboundSite2_' + str(i),'M_DnaA_c'), ('ssDNAunboundSite2_' + str(i+1),'ssDNABoundSite2_' + str(i)),k_on_2)
        sim.addReaction(('ssDNAunboundSite2_' + str(i+1), 'ssDNABoundSite2_' + str(i)),('ssDNAunboundSite2_' + str(i),'M_DnaA_c'),k_off)

    repInitProds = ['ssDNABoundSite2_30','ssDNA_Unbinding2_30','Initiator_C','Initiator_CC']
    
    sim.addReaction(('ssDNAunboundSite2_30', 'M_DnaA_c'),tuple(repInitProds),k_on_2)

    
    #Chromosme 3
    nonOricSpec = ['High_Affinity_Site3','High_Affinity_Bound3','High_Affinity_Site_oriC3',
                   'High_Affinity_Bound_oriC3','Low_Affinity_Site3_1','Low_Affinity_Site3_2',
                   'Low_Affinity_Bound3_1','Low_Affinity_Bound3_2']
    
    sim.addReaction(('High_Affinity_Site3', 'M_DnaA_c'),'High_Affinity_Bound3',k_high)
    
    sim.addReaction(('High_Affinity_Site_oriC3', 'M_DnaA_c'),('High_Affinity_Bound_oriC3','Low_Affinity_Site3_1'),k_high)
    sim.addReaction(('Low_Affinity_Site3_1', 'M_DnaA_c'),('Low_Affinity_Bound3_1','Low_Affinity_Site3_2'),k_low)
    sim.addReaction(('Low_Affinity_Site3_2', 'M_DnaA_c'),('Low_Affinity_Bound3_2','ssDNAunboundSite3_1'),k_low)

    species = []

    for i in range(30):         #loop adds 30 terms for unbound sites
        term = 'ssDNAunboundSite3_'
        term = term + str(i+1)
        species.append(term)
    for j in range(30):         #adds 30 more terms of Bound sites
        bnd = 'ssDNABoundSite3_'
        bnd = bnd + str(j+1)
        species.append(bnd)
    for k in range(30):
        unbnd = 'ssDNA_Unbinding3_'
        unbnd = unbnd + str(k+1)
        species.append(unbnd)

    species.append('Initiator_C3')
    species.append('Initiator_CC3')
    
    # Binding and unbinding reactions are constructed so that DnaA in the middle of the filament cannot unbind and
    # DnaA only bind next to the last DnaA in the filament.

    # Add the binding reactions for each of the 30 binding events in filament formation.
    for i in range (1,30):
        sim.addReaction(('ssDNAunboundSite3_' + str(i),'M_DnaA_c'), ('ssDNAunboundSite3_' + str(i+1),'ssDNABoundSite3_' + str(i)),k_on_2)
        sim.addReaction(('ssDNAunboundSite3_' + str(i+1), 'ssDNABoundSite3_' + str(i)),('ssDNAunboundSite3_' + str(i),'M_DnaA_c'),k_off)

    repInitProds = ['ssDNABoundSite3_30','ssDNA_Unbinding3_30','Initiator_C','Initiator_CC']
    
    sim.addReaction(('ssDNAunboundSite3_30', 'M_DnaA_c'),tuple(repInitProds),k_on_2)
    
    
    # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
    for i in range (2,31):
        sim.addReaction(('ssDNA_Unbinding_' + str(i), 'ssDNABoundSite_' + str(i)),('ssDNA_Unbinding_' + str(i-1),'M_DnaA_c'),helicase_removal_rate)

    sim.addReaction(('ssDNA_Unbinding_1', 'ssDNABoundSite_1'),'M_DnaA_c',helicase_removal_rate)
    
    for i in range (2,31):
        sim.addReaction(('ssDNA_Unbinding2_' + str(i), 'ssDNABoundSite2_' + str(i)),('ssDNA_Unbinding2_' + str(i-1),'M_DnaA_c'),helicase_removal_rate)

    sim.addReaction(('ssDNA_Unbinding2_1', 'ssDNABoundSite2_1'),'M_DnaA_c',helicase_removal_rate)
    
    for i in range (2,31):
        sim.addReaction(('ssDNA_Unbinding3_' + str(i), 'ssDNABoundSite3_' + str(i)),('ssDNA_Unbinding3_' + str(i-1),'M_DnaA_c'),helicase_removal_rate)

    sim.addReaction(('ssDNA_Unbinding3_1', 'ssDNABoundSite3_1'),'M_DnaA_c',helicase_removal_rate)
    return 


def addReplicationTwo(sim,specDict,genome3A_DNA,ModelSpecies):    
    dna = genome3A_DNA[0]

    K0rep = 0.26e-3
    KDrep = 0.001
    kcatrep = 100
    
    datp = getConc(max(1,specDict['M_datp_c']),specDict)#0.009*2 # mM
    dttp = getConc(max(1,specDict['M_dttp_c']),specDict)#0.011*2 # mM
    dctp = getConc(max(1,specDict['M_dctp_c']),specDict)#0.006*2 # mM
    dgtp = getConc(max(1,specDict['M_dgtp_c']),specDict)#0.0035*2 # mM
    DNApol3 = 35*countToMiliMol # mM
    
    chromosome_C = ['chromosome_C']
    
    chromosome_CC = ['chromosome_CC']

    gene_list = []
    for i in range(len(dna.features)):
        if ('product' in dna.features[i].qualifiers.keys()):
            if dna.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
                gene_list.append(i)

    
    intergenic_list = []
    
    C_bp = 0
    CC_bp = 0

    position = 0
    
    CC_genes = []

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        if start < len(dna.seq)/2:
            geneSeq = Seq(str(dna.seq[position:end]))#, generic_dna)

            if start == 0:

                baseCount = defaultdict(int)
                for base in set(geneSeq):
                    baseCount[base] = geneSeq.count(base)


                n_tot = sum(list(baseCount.values()))

                C_bp = C_bp + n_tot

                NMono_A = baseCount["A"]

                NMono_C = baseCount["C"]

                NMono_G = baseCount["G"]

                NMono_T = baseCount["T"]

                NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]

                NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

                k_gene_dup = kcatrep/((1+K0rep/DNApol3)*(KDrep**2/datp/dttp) + NMonoSum + n_tot - 1)

                intergenic_region = locusTag+'_inter'

                geneID = locusTag + '_gene'

                if geneID not in ModelSpecies:
                    ModelSpecies.append(geneID)
                    spec = [geneID]

                intergenic = [intergenic_region]

                RepProd = [intergenic_region,geneID]

                for i in range(len(geneSeq)):
                    RepProd.append('ATP_DNArep')

                for i in range(NMono_A):
                    RepProd.append('dATP_DNArep')
                    RepProd.append('dTTP_DNArep')

                for i in range(NMono_T):
                    RepProd.append('dATP_DNArep')
                    RepProd.append('dTTP_DNArep')

                for i in range(NMono_G):
                    RepProd.append('dGTP_DNArep')
                    RepProd.append('dCTP_DNArep')
                for i in range(NMono_C):
                    RepProd.append('dCTP_DNArep')
                    RepProd.append('dGTP_DNArep')

                sim.addReaction(('chromosome_C','Initiator_C'), tuple(RepProd), k_gene_dup)

                intergenic_list.append(intergenic_region)

                position  = end

            else:

                baseCount = defaultdict(int)
                for base in set(geneSeq):
                    baseCount[base] = geneSeq.count(base)


                n_tot = sum(list(baseCount.values()))

                C_bp = C_bp + n_tot

                NMono_A = baseCount["A"]

                NMono_C = baseCount["C"]

                NMono_G = baseCount["G"]

                NMono_T = baseCount["T"]
                NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]

                NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

                k_gene_dup = kcatrep/ (NMonoSum + n_tot - 1)

                intergenic_region = locusTag+'_inter' 

                geneID = locusTag + '_gene'

                if geneID not in ModelSpecies:
                    ModelSpecies.append(geneID)
                    spec = [geneID]

                intergenic = [intergenic_region]

                RepProd = [intergenic_region,geneID]

                for i in range(len(geneSeq)):
                    RepProd.append('ATP_DNArep')

                for i in range(NMono_A):
                    RepProd.append('dATP_DNArep')
                    RepProd.append('dTTP_DNArep')

                for i in range(NMono_T):
                    RepProd.append('dATP_DNArep')
                    RepProd.append('dTTP_DNArep')

                for i in range(NMono_G):
                    RepProd.append('dGTP_DNArep')
                    RepProd.append('dCTP_DNArep')
                for i in range(NMono_C):
                    RepProd.append('dCTP_DNArep')
                    RepProd.append('dGTP_DNArep')

                sim.addReaction(intergenic_list[-1], tuple(RepProd), k_gene_dup)

                intergenic_list.append(intergenic_region)

                position = end


    position = 543086

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        if start > len(dna.seq)/2:
    
            CC_genes.append(gene)

    CC_genes.reverse()

    for gene in CC_genes:

        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        geneSeq = Seq(str(dna.seq[start:position]))#, generic_dna)

        if end == 543086:

            baseCount = defaultdict(int)
            for base in set(geneSeq):
                baseCount[base] = geneSeq.count(base)


            n_tot = sum(list(baseCount.values()))

            CC_bp = CC_bp + n_tot

            NMono_A = baseCount["A"]

            NMono_C = baseCount["C"]

            NMono_G = baseCount["G"]

            NMono_T = baseCount["T"]

            NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]

            NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

            k_gene_dup = kcatrep/((1+K0rep/DNApol3)*(KDrep**2/datp/dttp) + NMonoSum + n_tot - 1)

            intergenic_region = locusTag+'_inter'

            geneID = locusTag + '_gene'

            if geneID not in ModelSpecies:
                ModelSpecies.append(geneID)
                spec = [geneID]
            intergenic = [intergenic_region]

            RepProd = [intergenic_region,geneID]

            for i in range(len(geneSeq)):
                RepProd.append('ATP_DNArep')

            for i in range(NMono_A):
                RepProd.append('dATP_DNArep')
                RepProd.append('dTTP_DNArep')

            for i in range(NMono_T):
                RepProd.append('dATP_DNArep')
                RepProd.append('dTTP_DNArep')

            for i in range(NMono_G):
                RepProd.append('dGTP_DNArep')
                RepProd.append('dCTP_DNArep')

            sim.addReaction(('chromosome_CC','Initiator_CC'), tuple(RepProd), k_gene_dup)

            intergenic_list.append(intergenic_region)

            position  = start


        elif dna.features[gene].qualifiers['locus_tag'][0] == 'JCVISYN3A_0421':

            print('End of Replication')

            baseCount = defaultdict(int)
            for base in set(geneSeq):
                baseCount[base] = geneSeq.count(base)

            n_tot = sum(list(baseCount.values()))

            CC_bp = CC_bp + n_tot

            NMono_A = baseCount["A"]
            NMono_C = baseCount["C"]

            NMono_G = baseCount["G"]

            NMono_T = baseCount["T"]

            NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]

            NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

            k_gene_dup = kcatrep/ (NMonoSum + n_tot - 1)

            intergenic_region = locusTag+'_inter' 

            geneID = locusTag + '_gene'

            if geneID not in ModelSpecies:
                ModelSpecies.append(geneID)
                spec = [geneID]

            intergenic = [intergenic_region]

            gene_rep_end_products = ['High_Affinity_Site_oriC2','High_Affinity_Site_oriC3',
                                     'chromosome_C','chromosome_CC','chromosome_C','chromosome_CC',
                                     intergenic_region,geneID]

            for i in range(16):

                gene_rep_end_products.append('High_Affinity_Site2')

            for i in range(len(geneSeq)):
                gene_rep_end_products.append('ATP_DNArep')

            for i in range(NMono_A):
                gene_rep_end_products.append('dATP_DNArep')
                gene_rep_end_products.append('dTTP_DNArep')

            for i in range(NMono_T):
                gene_rep_end_products.append('dATP_DNArep')
                gene_rep_end_products.append('dTTP_DNArep')
            for i in range(NMono_G):
                gene_rep_end_products.append('dGTP_DNArep')
                gene_rep_end_products.append('dCTP_DNArep')

            for i in range(NMono_C):
                gene_rep_end_products.append('dCTP_DNArep')
                gene_rep_end_products.append('dGTP_DNArep')

            sim.addReaction(intergenic_list[-1], tuple(gene_rep_end_products), k_gene_dup)

            intergenic_list.append(intergenic_region)

            position = start


        else:


            baseCount = defaultdict(int)
            for base in set(geneSeq):
                    baseCount[base] = geneSeq.count(base)


            n_tot = sum(list(baseCount.values()))

            CC_bp = CC_bp + n_tot

            NMono_A = baseCount["A"]

            NMono_C = baseCount["C"]

            NMono_G = baseCount["G"]

            NMono_T = baseCount["T"]

            NMonoDict = [NMono_A,NMono_C,NMono_G,NMono_T]
            NMonoSum = NMono_A*KDrep/datp + NMono_C*KDrep/dctp + NMono_T*KDrep/dttp + NMono_G*KDrep/dgtp

            k_gene_dup = kcatrep/ (NMonoSum + n_tot - 1)

            intergenic_region = locusTag+'_inter' 

            geneID = locusTag + '_gene'

            if geneID not in ModelSpecies:
                ModelSpecies.append(geneID)
                spec = [geneID]

            intergenic = [intergenic_region]
            RepProd = [intergenic_region,geneID]

            for i in range(len(geneSeq)):
                RepProd.append('ATP_DNArep')

            for i in range(NMono_A):
                RepProd.append('dATP_DNArep')
                RepProd.append('dTTP_DNArep')

            for i in range(NMono_T):
                RepProd.append('dATP_DNArep')
                RepProd.append('dTTP_DNArep')

            for i in range(NMono_G):
                RepProd.append('dGTP_DNArep')
                RepProd.append('dCTP_DNArep')
            for i in range(NMono_C):
                RepProd.append('dCTP_DNArep')
                RepProd.append('dGTP_DNArep')

            sim.addReaction(intergenic_list[-1], tuple(RepProd), k_gene_dup)

            intergenic_list.append(intergenic_region)

            position = start

    return

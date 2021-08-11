#### Cell Geometry ####

from jLM.RegionBuilder import RegionBuilder
import jLM

import pandas as pd
import numpy as np
import random

def buildRegions(sim, N_edges, N_2, sim_center, ptn_ratio, dna_monomers, cyto_radius, riboFile, dnaFile, filename, pmap, PartIdxMap, partIdx):
    ext = sim.region("extracellular")
    mem = sim.region("membrane")
    cyt = sim.region("cytoplasm")
    ribo = sim.region("ribosomes")
    dna = sim.region("DNA")
    she = sim.region("outer_cytoplasm")

    build = RegionBuilder(sim)

    cytoplasm = build.ellipsoid(radius = cyto_radius, center = sim_center)
    cyto_dilation = build.dilate(cytoplasm, se = build.se26)
    cyto_shell = cyto_dilation & ~cytoplasm
    cyto_dilation = build.dilate(cyto_dilation, se = build.se26)
    membrane = cyto_dilation & ~cyto_shell & ~cytoplasm
    extracellular = ~cyto_dilation
    
#     dnaFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/08262020/s1c8_free_origin_rep00004_CG_coords.dat' #CG file
    # dnaFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/10062020/s1c8_free_origin_reps00051_00080_CG/s1c8_free_origin_rep00052_CG_coords.dat'
    # dnaFile = '../model_data/s1c15/s1c15_base_rep00001_CG_coords.dat'
#     dnaFile = '../model_data/s1c15/s1c15_base_CG_reps00001_00090/s1c15_base/CG/s1c15_base_rep00070_CG_coords.dat'
    
#     dnaFile = configs[rep][0]
#     print(dnaFile)
    
    dnaDF = pd.read_csv(dnaFile, header = None)
#     dnaDF

    geneBlocks = []
    genePoints = []
    genome_placement = []

    # 7348

    for index, row in dnaDF.iterrows():
        x = row[0] + N_2
        y = row[1] + N_2
        z = row[2] + N_2

        genome_placement.append([x,y,z])

    # genePoints.append(genome_placement[-1])

    for i in range(len(genome_placement)):

    #     if i<7350:
            gene_coord = genome_placement[i]
            genePoints.append(gene_coord)
            
            
    genes = np.full((N_edges, N_edges, N_edges), False)
    
    for gene in genePoints:
    #     print(gene)
        x = int(gene[0])
        y = int(gene[1])
        z = int(gene[2])

        genes[x,y,z] = True
        
        
    riboDF = pd.read_csv(riboFile, header = None)
    
    ribosome_radius = 1e-8/sim.latticeSpacing

    ribosome_spheres = []
    ribo_points_x = []
    ribo_points_y = []
    ribo_points_z = []

    its = 0
    ribo_points = []
    ribo_center_points = []

    for index, row in riboDF.iterrows():
    #     if index < 7348:
        x_int = row[1] + N_2
        y_int = row[2] + N_2
        z_int = row[3] + N_2
    
#         x_int = row[2]//2 + N_2
#         y_int = row[3]//2 + N_2
#         z_int = row[4]//2 + N_2

        center_point = [x_int,y_int,z_int]

        riboPoints = []

        xpoint1 = [x_int+1,y_int,z_int]
        xpoint2 = [x_int-1,y_int,z_int]
        ypoint1 = [x_int,y_int+1,z_int]
        ypoint2 = [x_int,y_int-1,z_int]
        zpoint1 = [x_int,y_int,z_int+1]
        zpoint2 = [x_int,y_int,z_int-1]

        riboPoints.append(center_point)
        riboPoints.append(xpoint1)
        riboPoints.append(ypoint1)
        riboPoints.append(zpoint1)
        riboPoints.append(xpoint2)
        riboPoints.append(ypoint2)
        riboPoints.append(zpoint2)

        ribo_points.append(center_point)
        ribo_center_points.append(center_point)

        ribo_points.append(xpoint1)
        ribo_points.append(ypoint1)
        ribo_points.append(zpoint1)
        ribo_points.append(xpoint2)
        ribo_points.append(ypoint2)
        ribo_points.append(zpoint2)

    #     ribosome = build.ellipsoid(radius = ribosome_radius, center = center_point)

    #     ribosome_spheres.append(ribosome)
    #     ribo_points_x.append(x_int)
    #     ribo_points_y.append(y_int)
    #     ribo_points_z.append(z_int)

    ribo_center_points = np.array(ribo_center_points,dtype=np.int)
    
    ribosomes = np.full((N_edges, N_edges, N_edges), False)
    for coord in ribo_points:
        x = int(coord[0])
        y = int(coord[1])
        z = int(coord[2])

        ribosomes[x,y,z] = True
        
    build.compose(
    (ext, extracellular),
    (cyt, cytoplasm),
    (she, cyto_shell),
    (ribo, ribosomes),
    (mem, membrane),
    (dna, genes))
    
    degradosome = sim.species('Degradosome')
    PartIdxMap['Degradosome'] = partIdx
    partIdx = partIdx + 1
    
    sim, occupied_mem_spaces, membrane_spaces = addDegParticles(sim, pmap, N_edges, ptn_ratio, cyto_shell, ribo_points, genePoints, degradosome)
    
    partIdx = addsecY(sim, occupied_mem_spaces, membrane_spaces, ptn_ratio, PartIdxMap, partIdx)
    
#     sim.showRegionStack(scl=8)
    
    print('Geometry constructed')
    
    return sim, genePoints, ribo_points, ribo_center_points, ext, mem, cyt, ribo, dna, she, cyto_shell, partIdx


def addDegParticles(sim, pmap, N_edges, ptn_ratio, cyto_shell, ribo_points, genePoints, degradosome):
    
    membrane_spaces = []

    for i in range(N_edges):
        for j in range(N_edges):
            for k in range(N_edges):

                if cyto_shell[i][j][k]:
                    membrane_spaces.append([i,j,k])
    
    deg_coords = []
    occupied_mem_spaces = []
    deg_num = 0

    for point in ribo_points:
        occupied_mem_spaces.append(point)

    for point in genePoints:
        occupied_mem_spaces.append(point)

    while deg_num < 120*ptn_ratio: #*700/500
        position = random.choice(membrane_spaces)
        if position not in occupied_mem_spaces:
            deg_num = deg_num + 1
            occupied_mem_spaces.append(position)
            deg_coords.append(position)

        else:
#             print(deg_num)
            continue

    for coord in deg_coords:

        x = int(coord[0])
        y = int(coord[1])
        z = int(coord[2])

        degradosome.placeParticle(x,y,z,1)

    print('Degradosomes placed ',deg_num)
    
    return sim, occupied_mem_spaces, membrane_spaces

def addsecY(sim, occupied_mem_spaces, membrane_spaces, ptn_ratio, PartIdxMap, partIdx):
    secy_coords = []
    # occupied_mem_spaces = []
    secy_num = 0

    # for point in ribo_points:
    #     occupied_mem_spaces.append(point)

    # for point in genePoints:
    #     occupied_mem_spaces.append(point)

    secy = sim.species('secY')
    PartIdxMap['secY'] = partIdx
    partIdx = partIdx + 1

#     secy.diffusionRate(mem,sim.diffusionZero)
#     secy.diffusionRate(ext,sim.diffusionZero)
#     secy.diffusionRate(ribo,sim.diffusionZero)
#     secy.diffusionRate(dna,sim.diffusionZero)
#     secy.diffusionRate(cyt,sim.diffusionZero)
#     secy.diffusionRate(she,sim.diffusionZero)

    while secy_num < 66*ptn_ratio:
        position = random.choice(membrane_spaces)
        if position not in occupied_mem_spaces:
            secy_num = secy_num + 1
            occupied_mem_spaces.append(position)
            secy_coords.append(position)

        else:
            print(secy_num)
            continue

    for coord in secy_coords:

        x = int(coord[0])
        y = int(coord[1])
        z = int(coord[2])

        secy.placeParticle(x,y,z,1)
        
    return partIdx

    # print(len(occupied_mem_spaces))


def readDNAoccupancies(dnaPartFile):
#     dnaPartFile = configs[rep][1]
#     dnaPartFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/08262020/s1c8_free_origin_rep00004_FG_nodes.dat' #CG occupancies
    # dnaPartFile = '../model_data/syn3A_DNA_4nm_lattice_FGtoCG/10062020/s1c8_free_origin_reps00051_00080_CG/s1c8_free_origin_rep00052_FG_nodes.dat'
    # dnaPartFile = '../model_data/s1c15/s1c15_base_rep00001_FG_nodes.dat'
#     dnaPartFile = '../model_data/s1c15/s1c15_base_CG_reps00001_00090/s1c15_base/CG/s1c15_base_rep00070_FG_nodes.dat'
    print(dnaPartFile)

    dnaPartDF = pd.read_csv(dnaPartFile, header = None)

    geneOccupancies = []

    for index, row in dnaPartDF.iterrows():

        occupancy = []

        for position in row:
            if position != -1:
                occupancy.append(position)

        geneOccupancies.append(occupancy)

    print('DNA occupancies read')
    
    return geneOccupancies


def mapDNA(gene_starts, dna_monomers):
    genes_added = 0
    position_added = 0

    DNA_map = []

    locus_added = []

    locus_finished = []

    locus_back = gene_starts[genes_added-1][0]
    locusTag = gene_starts[genes_added][0]
    start = gene_starts[genes_added][1]

    for i in range(dna_monomers-len(gene_starts)):

        if (position_added >= start - 11.9) and (locusTag not in locus_added):

            locus_added.append(locusTag)

            locusNum = locusTag.split('3A_')[1].lstrip('0')

            if locusTag == 'JCVISYN3A_0910':

                DNA_map.append('g_909')

                DNA_map.append('E_910')

                position_added = position_added + 11.9*2

    #                 print('Added ' + geneMetID)

            elif locusTag == 'JCVISYN3A_0001':

                locus_back = locusTag

                genes_added = genes_added + 1

                locusTag = gene_starts[genes_added][0]

                start = gene_starts[genes_added][1]

                DNA_map.append('g_1')

                position_added = position_added + 11.9

    #                 print('Added ' + geneMetID)

            else:

                locus_back = locusTag

                genes_added = genes_added + 1

                locusTag = gene_starts[genes_added][0]

                start = gene_starts[genes_added][1]

                direction = gene_starts[genes_added-1][2]

                locus_end = gene_starts[genes_added-2][0]

                direction_end = gene_starts[genes_added-2][2]

                locusNumEnd = locus_end.split('3A_')[1].lstrip('0')

                if direction_end == 1:

                    endMetID = 'E_' + locusNumEnd

                    DNA_map.append(endMetID)

                elif direction_end == -1:

                    endMetID = 'g_' + locusNumEnd

                    DNA_map.append(endMetID)

                if direction == 1:

                    geneMetID = 'g_' + locusNum

                    DNA_map.append(geneMetID)

                elif direction == -1:

                    geneMetID = 'E_' + locusNum

                    DNA_map.append(geneMetID)

                position_added = position_added + 11.9*2

    #                 print('Added ' + geneMetID)


        else:

            if (locusTag == 'JCVISYN3A_0910') and (locusTag in locus_added):

                locusNum = locusTag.split('3A_')[1].lstrip('0')

                intergeneMetID = 'C_' + locusNum

                position_added = position_added + 11.9

                DNA_map.append(intergeneMetID)

    #             print('Added ' + intergeneMetID)

            else:

                locusNum = locus_back.split('3A_')[1].lstrip('0')

                intergeneMetID = 'C_' + locusNum

                position_added = position_added + 11.9

                DNA_map.append(intergeneMetID)

    DNA_map.append('g_910')           
    print('DNA map written')
    
    return DNA_map


def addDNApart(sim, DNA_map, genePoints, geneOccupancies, ext, mem, cyt, ribo, dna, she, RDME_species_list, PartIdxMap, partIdx):
    
    geneEnds = {}
    geneStarts = {}
    
    for i in range(len(DNA_map)):

        for j in range(len(genePoints)):

            occupancy = geneOccupancies[j]

            if i+1 in occupancy:

                coord = genePoints[j]
                x = int(coord[0])
                y = int(coord[1])
                z = int(coord[2])

                geneMetID = DNA_map[i]
                
#                 print(i,geneMetID,occupancy,coord)
                
                if 'E_' in geneMetID:
                    geneEnds[geneMetID] = coord
                    RDME_species_list.append(geneMetID)
                    PartIdxMap[geneMetID] = partIdx
                    partIdx = partIdx + 1
                    
                if 'g_' in geneMetID:
                    geneStarts[geneMetID] = coord
                    RDME_species_list.append(geneMetID)
                    PartIdxMap[geneMetID] = partIdx
                    partIdx = partIdx + 1
                    
                if geneMetID not in RDME_species_list:
                    RDME_species_list.append(geneMetID)
                    PartIdxMap[geneMetID] = partIdx
                    partIdx = partIdx + 1

                geneSpecies = sim.species(geneMetID, annotation = geneMetID)

                geneSpecies.placeParticle(x,y,z,1)

                geneSpecies.diffusionRate(cyt,sim.diffusionZero)
                geneSpecies.diffusionRate(mem,sim.diffusionZero)
                geneSpecies.diffusionRate(ext,sim.diffusionZero)
                geneSpecies.diffusionRate(ribo,sim.diffusionZero)
                geneSpecies.diffusionRate(dna,sim.diffusionZero)

    #             print('Added ' + geneMetID,x,y,z)

                break
        
    print('DNA particles added') 
    
    return sim, geneEnds, geneStarts, partIdx
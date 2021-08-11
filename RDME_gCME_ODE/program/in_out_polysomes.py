"""
A file to deal with writing from the ODE Simulation back to the CME Simulation

Author: David Bianchi
"""

import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import importlib
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np

import pySTDLM.PostProcessing as pp

### CONSTANTS
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*3.14159*(r_cell)**3)*(1000) # for a spherical cell

# Lipid Groups and headgroup surface area compositions (values surface areas in nm^2)
# From various literature sources including: Jo et al Biophys. J. (2009), Bjorkborn et al Biophys J. (2010) ...
saDict = {
    'M_clpn_c':0.4,
    'M_chsterol_c':0.35, # 0.35, test for Chol. value smaller
    'M_sm_c':0.45,
    'M_pc_c':0.55,
    'M_pg_c':0.6,
    'M_galfur12dgr_c':0.6,
    'M_fa_c':0.5, # Should this be here??
    'M_12dgr_c':0.5, # Scale down, should cdp-dag be added???
    'M_pa_c':0.5,
    'M_cdpdag_c':0.5,
}

# Add code that was used to do this parsing
memProtList = ['JCVISYN3A_0005','JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 'JCVISYN3A_0034', 'JCVISYN3A_0060','JCVISYN3A_0095', 'JCVISYN3A_0113','JCVISYN3A_0114','JCVISYN3A_0116','JCVISYN3A_0117','JCVISYN3A_0132', 'JCVISYN3A_0143','JCVISYN3A_0146','JCVISYN3A_0164','JCVISYN3A_0165', 'JCVISYN3A_0166', 'JCVISYN3A_0167', 'JCVISYN3A_0168', 'JCVISYN3A_0169', 'JCVISYN3A_0195', 'JCVISYN3A_0196', 'JCVISYN3A_0197','JCVISYN3A_0235','JCVISYN3A_0239','JCVISYN3A_0248','JCVISYN3A_0249','JCVISYN3A_0296','JCVISYN3A_0304','JCVISYN3A_0314','JCVISYN3A_0317','JCVISYN3A_0326','JCVISYN3A_0332','JCVISYN3A_0338','JCVISYN3A_0345', 'JCVISYN3A_0346','JCVISYN3A_0371','JCVISYN3A_0372','JCVISYN3A_0379','JCVISYN3A_0388','JCVISYN3A_0398','JCVISYN3A_0399','JCVISYN3A_0411','JCVISYN3A_0425', 'JCVISYN3A_0426', 'JCVISYN3A_0427', 'JCVISYN3A_0428','JCVISYN3A_0439','JCVISYN3A_0440','JCVISYN3A_0478','JCVISYN3A_0481','JCVISYN3A_0505','JCVISYN3A_0516','JCVISYN3A_0601','JCVISYN3A_0639', 'JCVISYN3A_0641', 'JCVISYN3A_0642', 'JCVISYN3A_0643', 'JCVISYN3A_0652', 'JCVISYN3A_0685', 'JCVISYN3A_0686', 'JCVISYN3A_0691','JCVISYN3A_0696', 'JCVISYN3A_0706', 'JCVISYN3A_0707', 'JCVISYN3A_0708', 'JCVISYN3A_0774', 'JCVISYN3A_0777','JCVISYN3A_0778','JCVISYN3A_0779', 'JCVISYN3A_0787', 'JCVISYN3A_0789', 'JCVISYN3A_0790', 'JCVISYN3A_0791', 'JCVISYN3A_0792','JCVISYN3A_0795', 'JCVISYN3A_0797', 'JCVISYN3A_0822', 'JCVISYN3A_0827', 'JCVISYN3A_0830', 'JCVISYN3A_0835','JCVISYN3A_0836', 'JCVISYN3A_0839', 'JCVISYN3A_0852','JCVISYN3A_0870', 'JCVISYN3A_0872', 'JCVISYN3A_0876', 'JCVISYN3A_0878', 'JCVISYN3A_0879', 'JCVISYN3A_0881','JCVISYN3A_0908']


def getProtSA(pmap,memProtList):
    """
    Calculates the protein surface are contribution to the growing cell surface area by using a dictionary of known
    membrane proteins and using a predicted average membrane protein surfacea area of (35 nm^2).

    Parameters:
        pmap (LM pmap) - The Lattice Microbes particle map

        memProtList (list) - A list containing the Syn3A IDs of all membrane proteins

    Returns:
        saInt (integer) - The projected surface area contribution of membrane proteins (in nm^2)
    """

#     avgProtSA = 35.2 # nm^2, average protein surface area to produce expected 47% coverage for 7.3K membrane proteins
    avgProtSA = 24.75 # nm^2, average protein surface area to produce expected 47% coverage for 9.6K membrane proteins

    otherNamesDict = {'JCVISYN3A_0234':['crr','crr_P'],'JCVISYN3A_0779':['ptsg','ptsg_P'],'JCVISYN3A_0694':['ptsh','ptsh_P']} # ptsh, ptsi, crrr

    count = 0 # Count of number of membrane proteins
    for elementID in memProtList:
        # JCVISYN3A_0008
        # ptsi_P,MMSYN1_0233,0,0.95ptsh,MMSYN1_0694,1,0.05 ptsh_P,MMSYN1_0694,0,0.95crr,MMSYN1_0234,1,0.05crr_P,MMSYN1_0234,0,0.95ptsg,MMSYN1_0779,1,0.15ptsg_P,MMSYN1_0779,0,0.85
        if elementID in otherNamesDict.keys():
            elementIDList = otherNamesDict[elementID]
            for obj in elementIDList:
                count+= pmap[obj]
        else:
            elementID = 'M_PTN_' + elementID + '_c'
            count += pmap[elementID]
    upProtSA = int(count*avgProtSA) # nm^2
    #pmap['cellSA_Prot']=int(upProtSA)

    return upProtSA


def calcLipidToSA(pmap):
    """
    Calculates the Lipid headgroup contribution to the growing cell surface area
    
    Note: The cell "surface area" should only grow in relation to the outer cell surface area, not the inner edge
    of the bilayer. From previous work accounting for the thickness of the membrane bilayer being approximately 5 nm we can calculate the percentage
    of the total surface area on the outside of the membrane to be approximately 51.3% of the total membrane surface area
    
    Parameters:
        pmap (LM pmap) - the Lattice Microbes particle map
        
    Returns:
        saInt (int) - the "external" cell surface area in nm^2 rounded off to an integer.
    """
    
    saFlt = 0.0 # surface area float (nm^2)
    for key,val in saDict.items():
        saFlt += pmap[key]*val
        
    # Here we partition between inner and outer membrane surface area
    fracOutMem = 0.513 # 51.3% of the membrane surface is the outer layer of the membrane
    saFltOut = saFlt*fracOutMem
        
    saInt = int(round(saFltOut)) # Obtain the surface area as an integer nm^2 value for easy storage
   
    # Protein is assumed to take up 40% of membrane, 60% lipid or 2/3 of lipid value
    #protScale=0.667 # For 40% prot. surface area membrane
    #protScale=0.8868 # For 47% prot. surface area membrane
    #protScale=1.13 # For 53% prot. surface area membrane
    #protScale=0.8868 # For 47% prot. surface area membrane
    #protSA=protScale*saInt
    
    protSA = getProtSA(pmap,memProtList)

    # Assign the values
    pmap['CellSA_Lip']=saInt 
    pmap['CellSA']=saInt+protSA
    pmap['CellSA_Prot']=protSA

    # Output the values and recalculate readius and volume
    #print("Lipid SA: ",saInt)
    #print("Prot SA: ",protSA)
    #print("Mem SA: ",saInt+protSA)
    #print(" ")

    totalSAint = saInt+protSA
    cellRadius = ((totalSAint/4/np.pi)**(1/2))*1e-9
    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)

    #print("Cell Radius: ",cellRadius, " m")
    #print("Cell Volume: ",cellVolume, " L")

    return
    #return totalSAint

def calcCellVolume(pmap):
    
    # Update and Calculate the Lipid Surface area at the current time point
    #saInt = CalcLipidToSA(pmap)

    SurfaceArea = pmap['CellSA']
    
    cellRadius_calc = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
    cellRadius = min(cellRadius_calc,255e-9)
    #print('Radius',cellRadius)
    
    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
    #print('Volume',cellVolume)
    
    return cellVolume


def mMtoPart(conc,pmap):
    """
    Convert ODE concentrations to CME Particle Counts

    Arguments:
    conc (float): ODE chemical species concentration

    Returns:
    particle (int): particle integer number

    """
    
    cellVolume = calcCellVolume(pmap)
#     print('mmtopart',cellVolume)
    
    particle = int(round((conc/1000)*NA*cellVolume))
        
    return particle


def RDME_to_pmap(pmap, RDME_cts, rdme_sim, multiStatePtnDict, tRNAstateDict, RDME_species_list):
    

    for specID in RDME_species_list:
        
        try:
        
            new_count = RDME_cts['countBySpecies'][rdme_sim.species(specID)]

            pmap[specID] = new_count
            
        except:
            
            print(specID,' not in RDME')

            continue
            
                
    for key, dic in multiStatePtnDict.items():
        
        if dic['Type'] == 'cytoplasm':

            stateList = dic['States']

            new_count = RDME_cts['countBySpecies'][rdme_sim.species(key)]

            base_state_number = new_count 

            for stateID in stateList:

                if stateID != key:

                    base_state_number = base_state_number - pmap[stateID]

            pmap[key] = base_state_number
            
        elif dic['Type'] == 'ptsg':
            
            new_count = RDME_cts['countBySpecies'][rdme_sim.species('ptsg')]
            
            base_state_number = new_count - pmap['ptsg_P']
            
            pmap['ptsg'] = base_state_number
            
            
    for key, stateList in tRNAstateDict.items():
        
        new_count = RDME_cts['countBySpecies'][rdme_sim.species(key)]
        
        base_state_number = new_count 

        for stateID in stateList:

            if stateID != key:

                base_state_number = base_state_number - pmap[stateID]

        pmap[key] = base_state_number
        
        
    pmap['P_227'] = pmap['M_lpl_PdhC_c'] + pmap['M_dhlpl_PdhC_c'] + pmap['M_acdhlpl_PdhC_c']
        


def calc_RDME_costs(pmap, RDME_cts, rdme_sim, degDict, singleStatePtnDict, multiStatePtnDict):
    
    ### Updated costs from mRNA degradation ###
    for key, dic in degDict.items():
        
        new_count = RDME_cts['countBySpecies'][rdme_sim.species(key)]
        
        if new_count < pmap[key]:
            
            number_deg = pmap[key] - new_count
            
            for cost, val in dic.items():
                
                pmap[cost] = pmap[cost] + val*number_deg
                
#                 print(key,cost,pmap[cost])
    
    
    ### Update costs from translation of single-state proteins ###
    for key, dic in singleStatePtnDict.items():
        
        new_count = RDME_cts['countBySpecies'][rdme_sim.species(key)]
        
        if new_count > pmap[key]:
            
            number_made = new_count - pmap[key]
            
            for cost, val in dic.items():
                
                pmap[cost] = pmap[cost] + val*number_made
                
            print(key,'ATP_translat',pmap['ATP_translat'],pmap[key],new_count)
    
                
    ### Update costs from translation of multi-state proteins ###
    for key, dic in multiStatePtnDict.items():
                
        if dic['Type'] == 'membrane':
            
            stateList = dic['States']
            
            totalPtnPmap = 0
            totalPtnRdme = 0
            
            for stateID in stateList:
                
                totalPtnPmap = totalPtnPmap + pmap[stateID]
                
                new_count = RDME_cts['countBySpecies'][rdme_sim.species(stateID)]
                
                totalPtnRdme = totalPtnRdme + new_count
                
            
            if totalPtnRdme > totalPtnPmap:
                
                translat_costs = dic['Translat_costs']
                
                number_made = totalPtnRdme - totalPtnPmap
                
                for cost, val in translat_costs.items():
                    
                    pmap[cost] = pmap[cost] + val*number_made
                    
                print(key,'ATP_translat',pmap['ATP_translat'])
                    
            new_count_mem = RDME_cts['countBySpecies'][rdme_sim.species(key)]
            
            count_mem = pmap[key]
            
            if new_count_mem > count_mem:
                
                number_inserted = new_count_mem - count_mem
                
#                 print(key)
                
                transloc_per_protein = dic['Transloc_cost']
                
                pmap['ATP_transloc'] = pmap['ATP_transloc'] + transloc_per_protein*number_inserted
                    
            
        elif dic['Type'] == 'cytoplasm':
            
            stateList = dic['States']
            
            totalPtn = 0
            
            for stateID in stateList:
                
                totalPtn = totalPtn + pmap[stateID]
                
            new_count = RDME_cts['countBySpecies'][rdme_sim.species(key)]
            
            if new_count > totalPtn:
                
                print(new_count,totalPtn)
                
                for stateID in stateList:
                    print(stateID, pmap[stateID])
                
                translat_costs = dic['Translat_costs']
                
                number_made = new_count - totalPtn
                
                for cost, val in translat_costs.items():
                    
                    pmap[cost] = pmap[cost] + val*number_made
                    
                print(key,'ATP_translat',pmap['ATP_translat'])
                    
            base_state_number = new_count 
            
            for stateID in stateList:
                
                if stateID != key:
                    
                    base_state_number = base_state_number - pmap[stateID]
                    
            pmap[key] = base_state_number
         
        
        elif dic['Type'] == 'ptsg':
            
            stateList = dic['States']
            
            totalPtnPmap = 0
            totalPtnRdme = 0
            
            total_ptsg = RDME_cts['countBySpecies'][rdme_sim.species('ptsg')] + RDME_cts['countBySpecies'][rdme_sim.species('S_779')] + RDME_cts['countBySpecies'][rdme_sim.species('ptsg_C')]
            
            for stateID in stateList:
                
                totalPtnPmap = totalPtnPmap + pmap[stateID]
                
            
            if total_ptsg > totalPtnPmap:
                
                translat_costs = dic['Translat_costs']
                
                number_made = total_ptsg - totalPtnPmap
                
                for cost, val in translat_costs.items():
                    
                    pmap[cost] = pmap[cost] + val*number_made
                    
                    
            new_count_mem = RDME_cts['countBySpecies'][rdme_sim.species('ptsg')]
            
            count_mem = pmap['ptsg'] + pmap['ptsg_P']
            
            if new_count_mem > count_mem:
                
                number_inserted = new_count_mem - count_mem
                
                transloc_per_protein = dic['Transloc_cost']
                
                pmap['ATP_transloc'] = pmap['ATP_transloc'] + transloc_per_protein*number_inserted
                
                
    new_5S = RDME_cts['countBySpecies'][rdme_sim.species('M_rRNA_5S_c')]
    old_5S = pmap['M_rRNA_5S_c']
    
    if new_5S > old_5S:
        
        new_made = new_5S - old_5S

        print('New rRNA', new_5S, old_5S, new_made)
        
        pmap['ATP_trsc'] = pmap['ATP_trsc'] + (1451 + 1285 + 1202 + 900)*new_made
        pmap['ATP_rRNA'] = pmap['ATP_rRNA'] + 1451*new_made
        
        pmap['CTP_rRNA'] = pmap['CTP_rRNA'] + 1285*new_made
        pmap['GTP_rRNA'] = pmap['GTP_rRNA'] + 1202*new_made
        pmap['UTP_rRNA'] = pmap['UTP_rRNA'] + 900*new_made
        
        
# From Tyler Biopolymers
def deleteParticle(particles, x, y, z, pid):
    ps = np.array([p for p_ in particles[:,z,y,x,:] for p in p_ if p != pid])
    pps = 16 # Particles Per Site
#     pps = cfg.particlesPerSite # Don't know how Tyler added pps as argument to his config
    ps.resize(pps)
    ps = ps.reshape((particles.shape[0], particles.shape[4]))
    particles[:,z,y,x,:] = ps

def checkParticle(particles, x, y, z, pid):
#     print(x,y,z)
    ps = np.array([p for p_ in particles[:,z,y,x,:] for p in p_ if p == pid])
    if pid in ps:
#         print(ps)
        return True
    else:
        return False
    
# def checkParticle(particles, x, y, z, pid):
#     ps = np.array([p for p_ in particles[:,z,y,x,:] for p in p_ if p == pid])
#     if ps:
#         return True
# #         print(ps)


def updatePolysomes(pmap,rdme_sim,lattice,PartIdxMap,ordered_poly_ribo):
    
    particleLatt = lattice.getParticleLatticeView()
    
    for jcvi3AID in genomePtnLocDict:
        
        locusNum = jcvi3AID.split('3A_')[1].lstrip('0')
        
        poly_start = 'PS_' + locusNum
        
        if poly_start in pmap:
            
#             print('Checking polysomes')
            
            poly_start_idx = PartIdxMap[poly_start]
            poly_mid = 'PM_' + locusNum
            poly_mid_idx = PartIdxMap[poly_mid]
            poly_end = 'PE_' + locusNum
            poly_end_idx = PartIdxMap[poly_end]

            poly_start_done = 'PS_' + locusNum + '_d'
            poly_start_done_idx = PartIdxMap[poly_start_done]
            poly_mid_done = 'PM_' + locusNum + '_d'
            poly_mid_done_idx = PartIdxMap[poly_mid_done]

            start_ribosome = 'PolysomeStart'
            start_ribosome_idx = PartIdxMap[start_ribosome]
            mid_ribosome = 'PolysomeMid'
            mid_ribosome_idx = PartIdxMap[mid_ribosome]
            end_ribosome = 'PolysomeEnd'
            end_ribosome_idx = PartIdxMap[end_ribosome]
            
            RNA_ID = 'R_' + locusNum
            RNAidx = PartIdxMap[RNA_ID]
            
            if pmap[poly_start_done]>0 or pmap[poly_mid_done]>0:
                
                print('Updating polysomes for gene ',locusNum)
                
                for cluster in ordered_poly_ribo:
                    
#                     print(cluster)
                    
                    if len(cluster) == 2:
                        
                        start_ribo_coord = cluster[0][1]
                        end_ribo_coord = cluster[1][1]
                        
                        xe = int(end_ribo_coord[0])
                        ye = int(end_ribo_coord[1])
                        ze = int(end_ribo_coord[2])
                        
                        xs = int(start_ribo_coord[0])
                        ys = int(start_ribo_coord[1])
                        zs = int(start_ribo_coord[2])
                        
                        start_done = checkParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                        end_open = checkParticle(particleLatt, ze, ye, xe, end_ribosome_idx)
                        
                        if start_done and end_open:
                            
                            print(cluster)
                            
                            # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                            print('Polysome 2: moving mRNA from start')
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                            deleteParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                            lattice.addParticle(zs, ys, xs, start_ribosome_idx)
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                            
                            # Remove unoccupied end ribosome and replace it with an occupied one.
                            print('Polysome 2: binding mRNA to end')
                            print('End occ: ', lattice.getOccupancy(ze,ye,xe))
                            deleteParticle(particleLatt, ze, ye, xe, end_ribosome_idx)
                            print('End occ: ', lattice.getOccupancy(ze,ye,xe))
                            lattice.addParticle(ze, ye, xe, poly_end_idx)
                            print('End occ: ', lattice.getOccupancy(ze,ye,xe))
                            
                        elif start_done and not end_open:
                            
                            print(cluster)
                            
                            # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                            # Since end is occupied, place back unoccupied ribosome and mRNA
                            print('Polysome 2: End occupied, releasing mRNA')
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                            deleteParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                            lattice.addParticle(zs, ys, xs, start_ribosome_idx)
                            lattice.addParticle(zs, ys, xs, RNAidx)
                            print('Start occ: ', lattice.getOccupancy(zs,ys,xs))

                        
                    else:
                    
                        for riboP in range(len(cluster)-1):
                            
                            if riboP == 0:
                                
                                start_ribo_coord = cluster[0][1]
                                next_ribo_coord = cluster[1][1]

                                xn = int(next_ribo_coord[0])
                                yn = int(next_ribo_coord[1])
                                zn = int(next_ribo_coord[2])

                                xs = int(start_ribo_coord[0])
                                ys = int(start_ribo_coord[1])
                                zs = int(start_ribo_coord[2])

                                start_done = checkParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                                next_open = checkParticle(particleLatt, zn, yn, xn, mid_ribosome_idx)

                                if start_done and next_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    print('Polysome: moving mRNA from start')
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                                    deleteParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                                    lattice.addParticle(zs, ys, xs, start_ribosome_idx)
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))

                                    # Remove unoccupied next ribosome and replace it with an occupied one.
                                    print('Polysome: binding mRNA to next ribosome')
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))
                                    deleteParticle(particleLatt, zn, yn, xn, mid_ribosome_idx)
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))
                                    lattice.addParticle(zn, yn, xn, poly_mid_idx)
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))

                                elif start_done and not next_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    # Since next is occupied, place back unoccupied ribosome and mRNA
                                    print('Polysome: next ribosome occupied, releasing mRNA')
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                                    deleteParticle(particleLatt, zs, ys, xs, poly_start_done_idx)
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                                    lattice.addParticle(zs, ys, xs, start_ribosome_idx)
                                    lattice.addParticle(zs, ys, xs, RNAidx)
                                    print('Start occ: ', lattice.getOccupancy(zs,ys,xs))
                                    
                                    
                            elif (riboP == len(cluster)-2):
                                
                                current_ribo_coord = cluster[riboP][1]
                                end_ribo_coord = cluster[riboP+1][1]

                                xe = int(end_ribo_coord[0])
                                ye = int(end_ribo_coord[1])
                                ze = int(end_ribo_coord[2])

                                xc = int(current_ribo_coord[0])
                                yc = int(current_ribo_coord[1])
                                zc = int(current_ribo_coord[2])

                                current_done = checkParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                end_open = checkParticle(particleLatt, ze, ye, xe, end_ribosome_idx)

                                if current_done and end_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    print('Polysome: moving mRNA from current ribosome')
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    deleteParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    lattice.addParticle(zc, yc, xc, mid_ribosome_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))

                                    # Remove unoccupied end ribosome and replace it with an occupied one.
                                    print('Polysome: binding mRNA to end')
                                    print('End occ: ', lattice.getOccupancy(ze,ye,xe))
                                    deleteParticle(particleLatt, ze, ye, xe, end_ribosome_idx)
                                    print('End occ: ', lattice.getOccupancy(ze,ye,xe))
                                    lattice.addParticle(ze, ye, xe, poly_end_idx)
                                    print('End occ: ', lattice.getOccupancy(ze,ye,xe))

                                elif current_done and not end_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    # Since end is occupied, place back unoccupied ribosome and mRNA
                                    print('Polysome: End ribosome occupied, releasing mRNA')
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    deleteParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    lattice.addParticle(zc, yc, xc, mid_ribosome_idx)
                                    lattice.addParticle(zc, yc, xc, RNAidx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    
                                    
                            else:

                                current_ribo_coord = cluster[riboP][1]
                                next_ribo_coord = cluster[riboP+1][1]

                                xn = int(next_ribo_coord[0])
                                yn = int(next_ribo_coord[1])
                                zn = int(next_ribo_coord[2])

                                xc = int(current_ribo_coord[0])
                                yc = int(current_ribo_coord[1])
                                zc = int(current_ribo_coord[2])

                                current_done = checkParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                next_open = checkParticle(particleLatt, zn, yn, xn, mid_ribosome_idx)

                                if current_done and next_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    print('Polysome: moving mRNA from current ribosome')
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    deleteParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    lattice.addParticle(zc, yc, xc, mid_ribosome_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))

                                    # Remove unoccupied end ribosome and replace it with an occupied one.
                                    print('Polysome: binding mRNA to next ribosome')
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))
                                    deleteParticle(particleLatt, zn, yn, xn, mid_ribosome_idx)
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))
                                    lattice.addParticle(zn, yn, xn, poly_mid_idx)
                                    print('Next occ: ', lattice.getOccupancy(zn,yn,xn))

                                elif current_done and not next_open:
                                    
                                    print(cluster)

                                    # Remove ribosome that is done translating and replace it with an unoccupied ribosome
                                    # Since end is occupied, place back unoccupied ribosome and mRNA
                                    print('Polysome: next ribosome occupied, releasing mRNA')
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    deleteParticle(particleLatt, zc, yc, xc, poly_mid_done_idx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))
                                    lattice.addParticle(zc, yc, xc, mid_ribosome_idx)
                                    lattice.addParticle(zc, yc, xc, RNAidx)
                                    print('Curr occ: ', lattice.getOccupancy(zc,yc,xc))                        
                        
            
        else:
            
            continue


def placeNewRNA(pmap,rdme_sim,lattice,geneEnds,geneStarts,PartIdxMap,rtRNA_ID_dict):
    
    particleLatt = lattice.getParticleLatticeView()
    
    for jcvi3AID in genomePtnLocDict:
        
        locusNum = jcvi3AID.split('3A_')[1].lstrip('0')
        
        newRNA = 'New_mRNA_' + locusNum
        
        newRnaCnt = pmap[newRNA]
        
        if newRnaCnt>1:
            print(newRNA, ' greater than 1')
        
        if newRnaCnt>0:
            ## Get Location of gene end ##
            print('Trying to place mRNA',newRNA)
            
            gene_end = 'E_' + locusNum
            coord_end = geneEnds[gene_end]
            xe = int(coord_end[0])
            ye = int(coord_end[1])
            ze = int(coord_end[2])
            
            gene_start = 'g_' + locusNum
            coord_start = geneStarts[gene_start]
            xs = int(coord_start[0])
            ys = int(coord_start[1])
            zs = int(coord_start[2])
            
            if (lattice.getOccupancy(zs,ys,xs) < 15) and (lattice.getOccupancy(ze,ye,xe) < 14):
                # Add a product particle
#                 lattice.addParticle(x,y,z,5)
                print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))
                print(xe,ye,ze,'End Occ:',lattice.getOccupancy(ze,ye,xe))
#                 print(gene_end,lattice.getParticle(ze,ye,xe,PartIdxMap[gene_end]))
#                 RNApol = rdme_sim.species('RNApol')
#                 RNApol.placeParticle(xe,ye,ze,1)
#                 rnaID = 'RNA_' + locusNum
#                 RNA_part = rdme_sim.species(rnaID)
#                 RNA_part.placeParticle(xe,ye,ze,1)
#                 Trsc_Term = rdme_sim.species('Trsc_Term')
#                 Trsc_Term.placeParticle(xs,ys,zs,1)
                
                RNApolIdx = PartIdxMap['RNApol']
                lattice.addParticle(ze,ye,xe,RNApolIdx)
            
                RNAidx = PartIdxMap['R_' + locusNum]
                lattice.addParticle(ze,ye,xe,RNAidx)
                
                print(xe,ye,ze,'End Occ:',lattice.getOccupancy(ze,ye,xe))
                
                RNAPidx = PartIdxMap['RP_' + locusNum]
                deleteParticle(particleLatt, zs, ys, xs, RNAPidx)
                
                print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))
                
                geneIdx = PartIdxMap['g_' + locusNum]
                lattice.addParticle(zs,ys,xs,geneIdx)
                
                print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))

#                 TrscTermIdx = PartIdxMap['Trsc_Term']
#                 lattice.addParticle(zs,ys,xs,TrscTermIdx)

                pmap[newRNA] = 0
                
                print('Placed mRNA',newRNA)
                
            else:
                
                print('Site full, waiting to place mRNA until next communication time')


    for jcvi3AID in genomeRnaLocDict:
        
        locusNum = jcvi3AID.split('3A_')[1]
        
        if ('0067' not in locusNum) and ('0068' not in locusNum) and ('0069' not in locusNum) and ('0532' not in locusNum) and ('0533' not in locusNum) and ('0534' not in locusNum):
            
            locusNum = jcvi3AID.split('3A_')[1].lstrip('0')

            newRNA = 'New_RNA_' + locusNum

            newRnaCnt = pmap[newRNA]
            
#             print(pmap['RNAP_' + locusNum])

            if newRnaCnt>1:
                print(newRNA, ' greater than 1')

            if newRnaCnt>0:
                print(newRNA)
                ## Get Location of gene end ##
                print('Trying to place tRNA',newRNA)

                gene_end = 'E_' + locusNum
                coord_end = geneEnds[gene_end]
                xe = int(coord_end[0])
                ye = int(coord_end[1])
                ze = int(coord_end[2])

                gene_start = 'g_' + locusNum
                coord_start = geneStarts[gene_start]
                xs = int(coord_start[0])
                ys = int(coord_start[1])
                zs = int(coord_start[2])

                if (lattice.getOccupancy(zs,ys,xs) < 15) and (lattice.getOccupancy(ze,ye,xe) < 14):
                    # Add a product particle
    #                 lattice.addParticle(x,y,z,5)
                    print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))
                    print(xe,ye,ze,'End Occ:',lattice.getOccupancy(ze,ye,xe))
    #                 print(gene_end,lattice.getParticle(ze,ye,xe,PartIdxMap[gene_end]))
    #                 RNApol = rdme_sim.species('RNApol')
    #                 RNApol.placeParticle(xe,ye,ze,1)
    #                 rnaID = 'RNA_' + locusNum
    #                 RNA_part = rdme_sim.species(rnaID)
    #                 RNA_part.placeParticle(xe,ye,ze,1)
    #                 Trsc_Term = rdme_sim.species('Trsc_Term')
    #                 Trsc_Term.placeParticle(xs,ys,zs,1)

                    RNApolIdx = PartIdxMap['RNApol']
                    lattice.addParticle(ze,ye,xe,RNApolIdx)

                    rtRNA_ID = rtRNA_ID_dict[newRNA]

                    RNAidx = PartIdxMap[rtRNA_ID]
                    lattice.addParticle(ze,ye,xe,RNAidx)
                    
                    print(xe,ye,ze,'End Occ:',lattice.getOccupancy(ze,ye,xe))
                    
                    RNAPidx = PartIdxMap['RP_' + locusNum]
                    deleteParticle(particleLatt, zs, ys, xs, RNAPidx)
                    
                    print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))

                    geneIdx = PartIdxMap['g_' + locusNum]
                    lattice.addParticle(zs,ys,xs,geneIdx)
                    
                    print(xs,ys,zs,'Start Occ:',lattice.getOccupancy(zs,ys,xs))

#                     TrscTermIdx = PartIdxMap['Trsc_Term']
#                     lattice.addParticle(zs,ys,xs,TrscTermIdx)

                    pmap[newRNA] = 0

                    print('Placed RNA',newRNA)

                else:

                    print('Site full, waiting to place RNA until next communication time')
                    
#         else:
            
#             print(locusNum)

def writeOdeResults(pmap,model,res,rdme_sim,lattice,geneEnds,geneStarts,PartIdxMap):
    """
    Write results of ODE model simulation timestep back to the CME data structure (HDF5)

    Parameters:
    pmap (particle map): the CME particle map storing species counts data

    model (model obj.): The ODE simulation kinetic model object

    res (np ndarray): The final timestep of results from an ODE simulation

    Returns:
    None
    """

    # Get the list of metabolites
    mL = model.getMetList()

    ## Do mapping with enumeration over mL for index and name from metList

    # For loop iterating over the species and updating their counts with ODE results
    for ind in range(len(mL)):
#         print("Updating PC: ",mL[ind].getID())
        pmap[mL[ind].getID()] = mMtoPart(res[ind],pmap) # Assign updated species counts to particle map using species IDs
        
        
def updateGipCosts(pmap,model,res,rdme_sim,lattice,geneEnds,geneStarts,PartIdxMap):
    ### Calculate ATP Costs from RDME particles ###
    
    ATP_hydro_counters = ['ATP_translat','ATP_trsc','ATP_mRNAdeg','ATP_DNArep','ATP_transloc']
    
    NTP_counters = [['ATP_mRNA','M_atp_c'],['CTP_mRNA','M_ctp_c'],['UTP_mRNA','M_utp_c'],
                    ['GTP_mRNA','M_gtp_c'],['ATP_tRNA','M_atp_c'],['CTP_tRNA','M_ctp_c'],
                    ['UTP_tRNA','M_utp_c'],['GTP_tRNA','M_gtp_c'],['ATP_rRNA','M_atp_c'],
                    ['CTP_rRNA','M_ctp_c'],['UTP_rRNA','M_utp_c'],['GTP_rRNA','M_gtp_c']]
    
    NMP_counters = [['AMP_mRNAdeg','M_amp_c'],['UMP_mRNAdeg','M_ump_c'],
                    ['CMP_mRNAdeg','M_cmp_c'],['GMP_mRNAdeg','M_gmp_c']]
    
    dNTP_counters = [['dATP_DNArep','M_datp_c'],['dTTP_DNArep','M_dttp_c'],
                     ['dCTP_DNArep','M_dctp_c'],['dGTP_DNArep','M_dgtp_c']]
    
    aatrna_counters = [["FMET_cost","M_fmettrna_c","M_trnamet_c"]]#,
#                        ["ALA_cost","M_alatrna_c","M_trnaala_c"], ["ARG_cost","M_argtrna_c","M_trnaarg_c"],
#                        ["ASN_cost","M_asntrna_c","M_trnaasn_c"], ["ASP_cost","M_asptrna_c","M_trnaasp_c"],
#                        ["CYS_cost","M_cystrna_c","M_trnacys_c"], ["GLU_cost","M_glutrna_c","M_trnaglu_c"],
#                        ["GLN_cost","M_glntrna_c","M_trnagln_c"], ["GLY_cost","M_glytrna_c","M_trnagly_c"],
#                        ["HIS_cost","M_histrna_c","M_trnahis_c"], ["ILE_cost","M_iletrna_c","M_trnaile_c"],
#                        ["LEU_cost","M_leutrna_c","M_trnaleu_c"], ["LYS_cost","M_lystrna_c","M_trnalys_c"],
#                        ["MET_cost","M_mettrna_c","M_trnamet_c"], ["PHE_cost","M_phetrna_c","M_trnaphe_c"],
#                        ["PRO_cost","M_protrna_c","M_trnapro_c"], ["SER_cost","M_sertrna_c","M_trnaser_c"],
#                        ["THR_cost","M_thrtrna_c","M_trnathr_c"], ["TRP_cost","M_trptrna_c","M_trnatrp_c"],
#                        ["TYR_cost","M_tyrtrna_c","M_trnatyr_c"], ["VAL_cost","M_valtrna_c","M_trnaval_c"]]
                       
    
    for costID in ATP_hydro_counters:

#         print(pmap['M_atp_c'])
        if 'translat' in costID:
        
            costCnt = pmap[costID]
            print(costID,costCnt)
            gtpCnt = pmap['M_gtp_c']

            if costCnt > gtpCnt:

                pmap[costID] = costCnt - gtpCnt

                pmap['M_gdp_c'] = pmap['M_gdp_c'] + gtpCnt
                pmap['M_pi_c'] = pmap['M_pi_c'] + gtpCnt
                pmap['M_gtp_c'] = 0

            else:

                pmap['M_gtp_c'] = pmap['M_gtp_c'] - pmap[costID]
                pmap['M_gdp_c'] = pmap['M_gdp_c'] + pmap[costID]
                pmap['M_pi_c'] = pmap['M_pi_c'] + pmap[costID]

                pmap[costID] = 0
                
        
        else:
            
            costCnt = pmap[costID]
            print(costID,costCnt)
            atpCnt = pmap['M_atp_c']

            if costCnt > atpCnt:

                pmap[costID] = costCnt - atpCnt

                pmap['M_adp_c'] = pmap['M_adp_c'] + atpCnt
                pmap['M_pi_c'] = pmap['M_pi_c'] + atpCnt
                pmap['M_atp_c'] = 0

            else:

                pmap['M_atp_c'] = pmap['M_atp_c'] - pmap[costID]
                pmap['M_adp_c'] = pmap['M_adp_c'] + pmap[costID]
                pmap['M_pi_c'] = pmap['M_pi_c'] + pmap[costID]

                pmap[costID] = 0
    
    print('ATP: ',pmap['M_atp_c'])
    print('Pi: ',pmap['M_pi_c'])
    print('ADP: ',pmap['M_adp_c'])
    print('FDP: ',pmap['M_fdp_c'])
    
    #print('ATP: ',pmap['M_atp_c'])
    #print('Pi: ',pmap['M_pi_c'])
    #print('SA: ',pmap['CellSA'])
    #print('CHOL: ',pmap['M_chsterol_c'])
    #print('CLPN: ',pmap['M_clpn_c'])

    for cost in NTP_counters:
        
#         costID = cost[0]
#         metID  = cost[1]
        
#         pmap[metID] = pmap[metID] - pmap[costID]
#         pmap['M_ppi_c'] = pmap['M_ppi_c'] + pmap[costID]
#         pmap[costID] = 0
        
        costID = cost[0]
        metID  = cost[1]
        
        cost_count = pmap[costID]
        met_count = pmap[metID]
        
        if cost_count>met_count:
            
            print('Used more NTP than available: ' + metID)
            
            pmap[costID] = cost_count - met_count
            pmap[metID] = 0
            pmap['M_ppi_c'] = pmap['M_ppi_c'] + met_count
            
        else:
            
            pmap[metID] = pmap[metID] - pmap[costID]
            pmap['M_ppi_c'] = pmap['M_ppi_c'] + pmap[costID]
            pmap[costID] = 0
        
        
    for cost in dNTP_counters:
        
        costID = cost[0]
        metID  = cost[1]
        
        cost_count = pmap[costID]
        met_count = pmap[metID]
        
        if cost_count>met_count:
            
            print('Used more dNTP than available: ' + metID)
            
            pmap[costID] = cost_count - met_count
            pmap[metID] = 0
            pmap['M_ppi_c'] = pmap['M_ppi_c'] + met_count
            
        else:
            
            pmap[metID] = pmap[metID] - pmap[costID]
            pmap['M_ppi_c'] = pmap['M_ppi_c'] + pmap[costID]
            pmap[costID] = 0

        
    for recycle in NMP_counters:
        
        recycleID = recycle[0]
        metID     = recycle[1]
        
        pmap[metID] = pmap[metID] + pmap[recycleID]
        pmap[recycleID] = 0
        
        
    for cost in aatrna_counters:
        
        costID = cost[0]
        chargedID = cost[1]
        unchargedID= cost[2]

        cost_count = pmap[costID]
        charged_pool = pmap[chargedID]
        
        if cost_count>charged_pool:
            
            print('Used more charged trna than available: ',chargedID)
            
            pmap[costID] = cost_count - charged_pool
            pmap[unchargedID] = pmap[unchargedID] + charged_pool
            pmap[chargedID] = 0
            
        else:
            
            pmap[chargedID] = charged_pool - cost_count
            pmap[unchargedID] = pmap[unchargedID] + cost_count
            pmap[costID] = 0
            
   
    # Call the surface area calculation every write step, every comm. timestep
#     calcLipidToSA(pmap)

    return


def calcODEfluxes(time,solver,model,resStart,resFinal,totalTime,delt,simFolder):

    if (np.rint(time)/60).is_integer():
        minute = int(np.rint(time)/60)
        currentFluxes = solver.calcFlux(0, resStart )
#                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

        # Create list of reactions and fluxes
        fluxList = []
        for indx,rxn in enumerate(model.getRxnList()):
            fluxList.append( (rxn.getID(), currentFluxes[indx]) )

        fluxDF = pd.DataFrame(fluxList)

        fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'min_start.csv'

        fluxDF.to_csv(fluxFileName,header=False)

#             minute = int(int(time)/60)
        currentFluxes = solver.calcFlux(0, resFinal )
#                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

        # Create list of reactions and fluxes
        fluxList = []
        for indx,rxn in enumerate(model.getRxnList()):
            fluxList.append( (rxn.getID(), currentFluxes[indx]) )

        fluxDF = pd.DataFrame(fluxList)

        fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'min_end.csv'

        fluxDF.to_csv(fluxFileName,header=False)

        fluxFileName = simFolder+'fluxes/fluxDF_'+str(minute)+'min.csv'

        fluxDF.to_csv(fluxFileName,header=False)

        print('Saved fluxes at ' + str(minute) + ' minutes.')

    if time > totalTime-delt:
        print(time)
        finalFluxes = solver.calcFlux(0, resFinal )
#                         print("Is INF:", np.where( np.isinf(finalFluxes) ) )

        # Create list of reactions and fluxes
        fluxList = []
        for indx,rxn in enumerate(model.getRxnList()):
            fluxList.append( (rxn.getID(), finalFluxes[indx]) )

        fluxDF = pd.DataFrame(fluxList)

        fluxDF.to_csv(simFolder+'fluxDF_final.csv',header=False)

        print('Saved final fluxes.')


def updateGlobalCmeCnts(pmap,CSIMfilename):
    
    fHandle=pp.openLMFile(CSIMfilename)
    
    aaMetIDs = ["M_ala__L_c", "M_arg__L_c","M_glu__L_c","M_gln__L_c",
        "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_gly_c",
        "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c",
        "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c", "M_amet_c"]
    
    glugln = 'M_glutrnagln_c'
    glugln_enz = 'glutrnagln_enz'
    glugln_enz_atp = 'glutrnagln_enz_atp'
    glugln_enz_atp_aa = 'glutrnagln_enz_atp_gln'

    glugln_list = [glugln,glugln_enz,glugln_enz_atp,glugln_enz_atp_aa]
    
    for key, val in pmap.items(): 
    
        if 'New_' in key:
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
        if (key=='M_atp_c') or (key=='M_adp_c') or (key=='M_amp_c') or (key=='M_ppi_c') or (key=='M_pi_c'):
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
#         if ('UTP_' in key) or ('GTP_' in key) or ('CTP_' in key):
            
#             count_trace = pp.getSpecieTrace(fHandle, key)

#             pcount = count_trace[-1]
            
#             print(key,pcount)

#             pmap[key] = pcount
            
        if ('RP_' in key) or ('RP2_' in key):
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
        if ('_AA' in key) or ('ATP_' in key) or ('_ATP' in key) or ('_tRNA' in key) or ('P_' in key):
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
        if ('M_trna' in key) or ('trna_c' in key):
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
        if '_cost' in key:
            
            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
        
        if key in aaMetIDs:

            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
        
        if key in glugln_list:

            count_trace = pp.getSpecieTrace(fHandle, key)

            pcount = count_trace[-1]
            
#             print(key,pcount)

            pmap[key] = pcount
            
    pp.closeLMFile(fHandle)        
        

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
# cellVolume = CytoVolume

subvolume_vol = 1000*(8e-9)**3 # L

# print(cellVolume)

# Avogadro:
avgdr   = 6.022e23 # molec/mol
Avognum = avgdr

NaV = Avognum * subvolume_vol

countToMiliMol = 1000/(avgdr*CytoVolume)

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

ribo_init = 20*Ecoli_V*avgdr/60/6800

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


# def calcCellVolume(pmap):
    
#     SurfaceArea = pmap['CellSA']
    
#     cellRadius_calc = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
#     cellRadius = min(cellRadius_calc,255e-9)
    
#     cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
# #     print('Volume',cellVolume)
    
#     return cellVolume

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
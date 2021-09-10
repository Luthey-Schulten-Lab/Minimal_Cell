"""
A file to deal with writing from the ODE Simulation back to the CME Simulation

Author: David Bianchi
"""

import csv
import numpy as np
import pandas as pd

### CONSTANTS
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*3.14159*(r_cell)**3)*(1000) # for a spherical cell

# Lipid Groups and headgroup surface area compositions (values surface areas in nm^2)
# From various literature sources including: Jo et al Biophys. J. (2009), Bjorkborn et al Biophys J. (2010), Saenz et al PNAS (2012)
saDict = {
    'M_clpn_c':0.4,
    'M_chsterol_c':0.35, 
    'M_sm_c':0.45,
    'M_pc_c':0.55,
    'M_pg_c':0.6,
    'M_galfur12dgr_c':0.6,
    'M_12dgr_c':0.5, 
    'M_pa_c':0.5,
    'M_cdpdag_c':0.5,
}

# Add code that was used to do this parsing, recalculate Area/protein, remove non-needed ATPase subunits: delta, a, b
memProtList = ['JCVISYN3A_0005','JCVISYN3A_0008', 'JCVISYN3A_0009', 'JCVISYN3A_0010', 'JCVISYN3A_0011', 'JCVISYN3A_0030', 'JCVISYN3A_0034', 'JCVISYN3A_0060','JCVISYN3A_0095', 'JCVISYN3A_0113','JCVISYN3A_0114','JCVISYN3A_0116','JCVISYN3A_0117','JCVISYN3A_0132', 'JCVISYN3A_0143','JCVISYN3A_0146','JCVISYN3A_0164','JCVISYN3A_0165', 'JCVISYN3A_0166', 'JCVISYN3A_0167', 'JCVISYN3A_0168', 'JCVISYN3A_0169', 'JCVISYN3A_0195', 'JCVISYN3A_0196', 'JCVISYN3A_0197','JCVISYN3A_0235','JCVISYN3A_0239','JCVISYN3A_0248','JCVISYN3A_0249','JCVISYN3A_0296','JCVISYN3A_0304','JCVISYN3A_0314','JCVISYN3A_0317','JCVISYN3A_0326','JCVISYN3A_0332','JCVISYN3A_0338','JCVISYN3A_0345', 'JCVISYN3A_0346','JCVISYN3A_0371','JCVISYN3A_0372','JCVISYN3A_0379','JCVISYN3A_0388','JCVISYN3A_0398','JCVISYN3A_0399','JCVISYN3A_0411','JCVISYN3A_0425', 'JCVISYN3A_0426', 'JCVISYN3A_0427', 'JCVISYN3A_0428','JCVISYN3A_0439','JCVISYN3A_0440','JCVISYN3A_0478','JCVISYN3A_0481','JCVISYN3A_0505','JCVISYN3A_0516','JCVISYN3A_0601','JCVISYN3A_0639', 'JCVISYN3A_0641', 'JCVISYN3A_0642', 'JCVISYN3A_0643', 'JCVISYN3A_0652', 'JCVISYN3A_0685', 'JCVISYN3A_0686', 'JCVISYN3A_0691','JCVISYN3A_0696', 'JCVISYN3A_0706', 'JCVISYN3A_0707', 'JCVISYN3A_0708', 'JCVISYN3A_0774', 'JCVISYN3A_0777','JCVISYN3A_0778','JCVISYN3A_0779', 'JCVISYN3A_0787', 'JCVISYN3A_0789', 'JCVISYN3A_0790', 'JCVISYN3A_0791', 'JCVISYN3A_0792','JCVISYN3A_0795', 'JCVISYN3A_0797', 'JCVISYN3A_0822', 'JCVISYN3A_0827', 'JCVISYN3A_0830', 'JCVISYN3A_0835','JCVISYN3A_0836', 'JCVISYN3A_0839', 'JCVISYN3A_0852','JCVISYN3A_0870', 'JCVISYN3A_0872', 'JCVISYN3A_0876', 'JCVISYN3A_0878', 'JCVISYN3A_0879', 'JCVISYN3A_0881','JCVISYN3A_0908'] # _0793, _0794, _0796 = ATPase not in transmembrane


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

    avgProtSA = 28.0 # nm^2, average protein surface area to produce expected 54% coverage for 9.6K membrane proteins

    # Proteins with special naming conventions in the model
    otherNamesDict = {'JCVISYN3A_0779':['ptsg','ptsg_P']}

    count = 0 # Count of number of membrane proteins
    for elementID in memProtList:
        if elementID in otherNamesDict.keys():
            elementIDList = otherNamesDict[elementID]
            for obj in elementIDList:
                count+= pmap[obj]
        else:
            elementID = 'M_PTN_' + elementID + '_c'
            count += pmap[elementID]
    #print("Count is: ", count)
    upProtSA = int(count*avgProtSA) # nm^2

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
    
    protSA = getProtSA(pmap,memProtList)

    # Assign the values
    pmap['CellSA_Lip']=saInt 
    pmap['CellSA']=saInt+protSA
    pmap['CellSA_Prot']=protSA


    totalSAint = saInt+protSA
    cellRadius = ((totalSAint/4/np.pi)**(1/2))*1e-9
    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)

    return

def calcCellVolume(pmap):
    
    SurfaceArea = pmap['CellSA']
    
    cellRadius = ((SurfaceArea/4/np.pi)**(1/2))*1e-9
    #print('Radius',cellRadius)
    
    cellVolume = ((4/3)*np.pi*(cellRadius)**3)*(1000)
    #print('Volume',cellVolume)
    
    # Put a ceiling on cellV growth, stop the cell volume growth when the volume has doubled
    if (cellVolume > 6.70e-17):
        #print("stopping volume growth at 6.70e-17 L")
        cellVolume = 6.70e-17
        pmap['CellV'] = int(670)

    # Have a recorded of the cell volume times 10^19 L
    else:
        pmap['CellV'] = int(round(cellVolume*1e19))
    #print("Cell Volume is: ",cellVolume)

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
    
    particle = int(round((conc/1000)*NA*cellVolume))
        
    return particle


def writeResults(pmap,model,res,time,procid):
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
        if (mL[ind].getID() == 'CellSA') or (mL[ind].getID() == 'CellSA_Prot') or (mL[ind].getID() == 'CellSA_Lip'):
            continue
        else:
            pmap[mL[ind].getID()] = mMtoPart(res[ind],pmap) # Assign updated species counts to particle map using species IDs
        
        
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

    
    for costID in ATP_hydro_counters:

        if 'translat' in costID:

            costCnt = pmap[costID]
            atpCnt = pmap['M_gtp_c']

            if costCnt > atpCnt:

                pmap[costID] = costCnt - atpCnt

                pmap['M_gdp_c'] = pmap['M_gdp_c'] + atpCnt
                pmap['M_pi_c'] = pmap['M_pi_c'] + atpCnt
                pmap['M_gtp_c'] = 0

            else:

                pmap['M_gtp_c'] = pmap['M_gtp_c'] - pmap[costID]
                pmap['M_gdp_c'] = pmap['M_gdp_c'] + pmap[costID]
                pmap['M_pi_c'] = pmap['M_pi_c'] + pmap[costID]

                pmap[costID] = 0

        else:
            costCnt = pmap[costID]
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

    

    for cost in NTP_counters:
        
        
        costID = cost[0]
        metID  = cost[1]
        
        cost_count = pmap[costID]
        met_count = pmap[metID]
        
        if cost_count>met_count:
            
            #print('Used more NTP than available: ' + metID)
            
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
            
            #print('Used more dNTP than available: ' + metID)
            
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
            
            #print('Used more charged trna than available: ',chargedID)
            
            pmap[costID] = cost_count - charged_pool
            pmap[unchargedID] = pmap[unchargedID] + charged_pool
            pmap[chargedID] = 0
            
        else:
            
            pmap[chargedID] = charged_pool - cost_count
            pmap[unchargedID] = pmap[unchargedID] + cost_count
            pmap[costID] = 0

    # Call the surface area calculation every write step, every comm. timestep
    calcLipidToSA(pmap)

    print("Recalculated Lipid SA")

    # Out Met Csvs
    if (int(time)/60).is_integer():
        minute = int(int(time)/60)
        # print("Calling updated write to csvs at minute: ",minute)

    return

def outMetCsvs(pmap,minute,procID):
    """
    Write metabolite csvs at each timestep

    Parameters:
        pmap (particle map): the CME particle map storing species counts data

        minute (int): the minute of cell cycle simulation time

        procID (int): the process ID (necessary if running parallel/batch simulations)

    Returns:
        None
    """

    specIDs = []
    newCounts = []
    
    if int(minute) == -1:
        print("The -1 LOOP WAS ENTERED!!")
        for met in pmap.particleMap.keys():
            if (pmap.particleMap[met] == 'CellSA') or (pmap.particleMap[met] == 'CellSA_Prot') or (pmap.particleMap[met] == 'CellSA_Lip') or (pmap.particleMap[met] == 'CellV'):
                specIDs.append(met)
                newCounts.append(pmap.particleMap[met])
            else:
                specIDs.append(met)
                newCounts.append(pmap.particleMap[met])    

    else:
        for met in pmap.keys():
            if (pmap[met] == 'CellSA') or (pmap[met] == 'CellSA_Prot') or (pmap[met] == 'CellSA_Lip') or (pmap[met] == 'CellV'):
                specIDs.append(met)
                newCounts.append(pmap[met])
            else:
                specIDs.append(met)
                newCounts.append(pmap[met])

    if int(minute) == -1:
        minute = 1

    metsDF = pd.DataFrame()
    metsDF['Time'] = specIDs
    print("Species IDs len is: ",len(specIDs))
    metsDF[np.rint(int(minute))] = newCounts
    print("DF size on time: ",minute, " is: ", metsDF.size)
    print(metsDF)

    metFileName = '../simulations/rep-'+str(procID)+'.csv'
    if (int(minute) == 0):
        metsDF.to_csv(metFileName,index=False,mode='w+')
    else:
        print("writing df")
        metsDF.to_csv(metFileName,index=False,mode='a')

    return

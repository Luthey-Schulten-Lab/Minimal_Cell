"""
A file to define Sink Reactions

Author: David Bianchi
"""

import numpy as np
import pandas as pd
from in_out import calcCellVolume

### Constants
NA = 6.022e23 # Avogadro's
r_cell = 200e-9 # 200 nm radius, 400 nm diameter
V = ((4/3)*np.pi*(r_cell)**3)*(1000) # for a spherical cell

countToMiliMol = 1000/(NA*V)

def partTomM(particles,pmap):
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


def addSrcSinkRxns(model,pmap):
    """
    Adds source and sink reactions

    Parameters:
        model (odecell.model) - the simulation model object

        pmap - the Lattice Microbes particle map

    Returns:

        None
    """    


    SrcSinkRxnsList = []

#     rxnIndx = model.addReaction("pseudo_ATP_hydrolysis","zeroOrderOnOff","Pseudo ATP hydrolysis")
#     model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#     model.addProduct(rxnIndx,"Prod1","M_adp_c")
#     # model.addProduct(rxnIndx,"Prod2","M_pi_c")
#     hydro_rate = 0.686 - 0.047 - 0.024 - 0.46 - 0.08 - 0.046 # mM/s - ZT update #1.33 - 0.64 #mM/s
#     #print(hydro_rate)
#     # 0.80 mM/s (9.7 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": hydro_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)

   # rxnIndx = model.addReaction("ATP_AMP","zeroOrderOnOff","ATP to AMP")
   # model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
   # model.addProduct(rxnIndx,"Prod1","M_amp_c")
   # model.addProduct(rxnIndx,"Prod2","M_pi_c")
   # atp_amp_rate = 0.149 #mM/s
   # paramDict = {"K": atp_amp_rate}
   # for key,val in paramDict.items():
   #     model.addParameter(rxnIndx, key, val)

   # SrcSinkRxnsList.append(rxnIndx)

    #rxnIndx = model.addReaction("RNAdeg_AMP","zeroOrderOnOff","RNA deg AMP")
    #model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
    #model.addProduct(rxnIndx,"Prod1","M_amp_c")
    #model.addProduct(rxnIndx,"Prod2","M_pi_c")
    #RNAdeg_amp_rate = 0.021 #mM/s
    #paramDict = {"K": RNAdeg_amp_rate}
    #for key,val in paramDict.items():
        #model.addParameter(rxnIndx, key, val)

    #SrcSinkRxnsList.append(rxnIndx)

#     rxnIndx = model.addReaction("pi_c_src","zeroOrderOnOff","phosphate production")
#     model.addProduct(rxnIndx,"Prod1","M_pi_c")
#     pi_c_rate = 0.71 - 0.047 - 0.048 - 0.46 - 0.08 - 0.046 # mM/s - Zane update #1.38 - 0.64 #mM/s
#     # 0.85 mM/s (10.24 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": pi_c_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)
    
#     rxnIndx = model.addReaction("nad_c_src","zeroOrderOnOff","phosphate production")
#     model.addProduct(rxnIndx,"Prod1","M_nad_c")
#     nad_c_rate = 0.000175 # mM/s - Zane update #1.38 - 0.64 #mM/s
#     # 0.85 mM/s (10.24 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": nad_c_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)
    
#     rxnIndx = model.addReaction("NACabc","zeroOrderOnOff","nac production")
#     model.addProduct(rxnIndx,"Prod1","M_nac_c")
#     model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#     model.addProduct(rxnIndx,"Prod2","M_adp_c")
#     model.addProduct(rxnIndx,"Prod3","M_pi_c")
#     nad_c_rate = 0.000175 # mM/s - Zane update #1.38 - 0.64 #mM/s
#     # 0.85 mM/s (10.24 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": nad_c_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)
    
    
#     rxnIndx = model.addReaction("RIBFLVabc","zeroOrderOnOff","phosphate production")
#     model.addProduct(rxnIndx,"Prod1","M_ribflv_c")
#     model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#     model.addProduct(rxnIndx,"Prod2","M_adp_c")
#     model.addProduct(rxnIndx,"Prod3","M_pi_c")
#     nad_c_rate = 0.00000175 # mM/s - Zane update #1.38 - 0.64 #mM/s
#     # 0.85 mM/s (10.24 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": nad_c_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)
    
    
#     rxnIndx = model.addReaction("5FTHFabc","zeroOrderOnOff","phosphate production")
#     model.addProduct(rxnIndx,"Prod1","M_5fthf_c")
#     model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#     model.addProduct(rxnIndx,"Prod2","M_adp_c")
#     model.addProduct(rxnIndx,"Prod3","M_pi_c")
#     nad_c_rate = 0.000035 # mM/s - Zane update #1.38 - 0.64 #mM/s
#     # 0.85 mM/s (10.24 FBA) + 0.4 mM/s translation + 0.08 mM/s transcription + 0.05 mM/s RNA deg
#     paramDict = {"K": nad_c_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)

#     SrcSinkRxnsList.append(rxnIndx)
    

    for rxnIndx in SrcSinkRxnsList:
        if "onoff" in model.getReaction(rxnIndx).getKeys():
            model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

    sourceSinkRxns = []
    
#     RNA_sources = False
    
#     if RNA_sources:

#         rxnIndx = model.addReaction("atp_RNA_sink","zeroOrderOnOff","rna pol")
#         model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#         atp_c_rate = 175/6300 #mM/s
#         # print(atp_c_rate)
#         paramDict = {"K": atp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("utp_RNA_sink","zeroOrderOnOff","rna pol")
#         model.addSubstrate(rxnIndx,"Sub1","M_utp_c")
#         utp_c_rate = 175/6300 #mM/s
#         paramDict = {"K": utp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("gtp_RNA_sink","zeroOrderOnOff","rna pol")
#         model.addSubstrate(rxnIndx,"Sub1","M_gtp_c")
#         gtp_c_rate = 75/6300 #mM/s
#         paramDict = {"K": gtp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("ctp_RNA_sink","zeroOrderOnOff","rna pol")
#         model.addSubstrate(rxnIndx,"Sub1","M_ctp_c")
#         ctp_c_rate = 75/6300 #mM/s
#         paramDict = {"K": ctp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("amp_RNA_src","zeroOrderOnOff","rna deg")
#         model.addProduct(rxnIndx,"Prod1","M_amp_c")
#         amp_c_rate = 100/6300 #+ 0.14 #mM/s
#         paramDict = {"K": amp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("ump_RNA_src","zeroOrderOnOff","rna deg")
#         model.addProduct(rxnIndx,"Prod1","M_ump_c")
#         ump_c_rate = 100/6300 #mM/s
#         paramDict = {"K": ump_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("gmp_RNA_src","zeroOrderOnOff","rna deg")
#         model.addProduct(rxnIndx,"Prod1","M_gmp_c")
#         gmp_c_rate = 45/6300 #mM/s
#         paramDict = {"K": gmp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("cmp_RNA_src","zeroOrderOnOff","rna deg")
#         model.addProduct(rxnIndx,"Prod1","M_cmp_c")
#         cmp_c_rate = 45/6300 #+ 0.0025 #mM/s
#         paramDict = {"K": cmp_c_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#         rxnIndx = model.addReaction("ppi_RNA_src","zeroOrderOnOff","rna deg")
#         model.addProduct(rxnIndx,"Prod1","M_ppi_c")
#         ppi_c_rna_rate = 175/6300 + 175/6300 + 75/6300 + 75/6300 #0.05 #mM/s
#         paramDict = {"K": ppi_c_rna_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)
#         model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

#         SrcSinkRxnsList.append(rxnIndx)

#    rxnIndx = model.addReaction("ura_transport_supp","zeroOrderOnOff","Adding to ura transport")
#    model.addProduct(rxnIndx,"Prod1","M_ump_c")
#    ura_c_rate = 0.005 #mM/s
#    paramDict = {"K": ura_c_rate}
#    for key,val in paramDict.items():
#        model.addParameter(rxnIndx, key, val)
#    model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )
#
#    SrcSinkRxnsList.append(rxnIndx)

    #model.addMetabolite("M_uri_c","uridine",1.69)

    #rxnIndx = model.addReaction("URIabc","zeroOrderOnOff","URI active transport")
    #model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
    # model.addSubstrate(rxnIndx,"Sub2","M_uri_e")
    #model.addProduct(rxnIndx,"Prod1","M_adp_c")
    #model.addProduct(rxnIndx,"Prod2","M_pi_c")
    #model.addProduct(rxnIndx,"Prod3","M_uri_c")
    #uriabc_rate = 0.017 #mM/s
    #paramDict = {"K": uriabc_rate}
    #for key,val in paramDict.items():
        #model.addParameter(rxnIndx, key, val)
    
    #SrcSinkRxnsList.append(rxnIndx)

#     #print("Adding PIabc reaction, and appending to srcsink")
#     rxnIndx = model.addReaction("PIabc","zeroOrderOnOff","PI active transport")
#     model.addSubstrate(rxnIndx,"Sub1","M_atp_c")
#     # model.addSubstrate(rxnIndx,"Sub2","M_uri_e")
#     model.addProduct(rxnIndx,"Prod1","M_adp_c")
#     model.addProduct(rxnIndx,"Prod2","M_pi_c",stoich=2)
#     piabc_rate = 0.2 #0.024*3.0 #mM/s
#     paramDict = {"K": piabc_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)
#     model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )   
#     #print("rxnIndx: ", rxnIndx)
 
#     SrcSinkRxnsList.append(rxnIndx)


    ### Amino Acid charging and metabolism
#     aaMets = ["M_ala__L_c", "M_arg__L_c", 
#     "M_asn__L_c", "M_asp__L_c", "M_cys__L_c", "M_glu__L_c", "M_gln__L_c", "M_gly_c", 
#     "M_his__L_c", "M_ile__L_c", "M_leu__L_c", "M_lys__L_c", "M_met__L_c", "M_phe__L_c", 
#     "M_pro__L_c", "M_ser__L_c", "M_thr__L_c", "M_trp__L_c", "M_tyr__L_c", "M_val__L_c",]

#     aaTRNAMets = ["M_alatrna_c", "M_argtrna_c", 
#         "M_asntrna_c", "M_asptrna_c", "M_cystrna_c", "M_glutrna_c", "M_glntrna_c", "M_glytrna_c", 
#         "M_histrna_c", "M_iletrna_c", "M_leutrna_c", "M_lystrna_c", "M_mettrna_c", "M_phetrna_c", 
#         "M_protrna_c", "M_sertrna_c", "M_thrtrna_c", "M_trptrna_c", "M_tyrtrna_c", "M_valtrna_c"]

#     aaTRNAFreeMets = ["M_trnaala_c", "M_trnaarg_c", 
#         "M_trnaasn_c", "M_trnaasp_c", "M_trnacys_c", "M_trnaglu_c", "M_trnagln_c", "M_trnagly_c", 
#         "M_trnahis_c", "M_trnaile_c", "M_trnaleu_c", "M_trnalys_c", "M_trnamet_c", "M_trnaphe_c", 
#         "M_trnapro_c", "M_trnaser_c", "M_trnathr_c", "M_trnatrp_c", "M_trnatyr_c", "M_trnaval_c"]

#     AA_sources = False
    
#     if AA_sources:
        
#         aaTransport = [['R_ARGt2r', 'M_arg__L_c', 0.0047],
#          ['R_ASPt2pr', 'M_asp__L_c', 0.009],
#          ['R_GLYt2r', 'M_gly_c', 0.0081],
#          ['R_ISOt2r', 'M_ile__L_c', 0.014],
#          ['R_ALAt2r', 'M_ala__L_c', 0.0090],
#          ['R_ASNt2r', 'M_asn__L_c', 0.0105],
#          ['R_LEUt2r', 'M_leu__L_c', 0.014],
#          ['R_HISt2r', 'M_his__L_c', 0.0020],
#          ['R_LYSt2r', 'M_lys__L_c', 0.0151],
#          ['R_PROt2r', 'M_pro__L_c', 0.0045],
#          ['R_PHEt2r', 'M_phe__L_c', 0.0064],
#          ['R_THRt2r', 'M_thr__L_c', 0.0086],
#          ['R_TRPt2r', 'M_trp__L_c', 0.0012],
#          ['R_TYRt2r', 'M_tyr__L_c', 0.0051],
#          ['R_VALt2r', 'M_val__L_c', 0.0095],
#          ['R_SERt2r', 'M_ser__L_c', 0.0089],
#          ['R_METt2r', 'M_met__L_c', 0.0031],
#          ['R_CYSt2r', 'M_cys__L_c', 0.00084],
#          ['R_GLUt2pr', 'M_glu__L_c', 0.0035], #0.0023],
#          ['R_GLNt2r', 'M_gln__L_c', 0.013]] #0.014]]

#     # ['R_GLNt2r', 'M_gln__L_c', 0.0084]['R_GLUt2pr', 'M_glu__L_c', 0.0020]

#     # aaTransport

#         for rxn in aaTransport:
#             #print(rxn[0])
#             rxnID = rxn[0] + "_zero"
#             metID = rxn[1]
#             rate = rxn[2]

#             if ('glu' not in metID) and ('gln' not in metID):
#                 model.addMetabolite(metID, metID, partTomM(pmap[metID],pmap))

#             # Add each of th eintracellular amino acids
#             #if ((metID != "M_glu__L_c") and (metID != "M_gln__L_c")):
#                 #model.addMetabolite(metID,metID,int(round(countToMiliMol*pmap[metID])))    
#             rxnIndx = model.addReaction(rxnID,"zeroOrderOnOff",metID + ' Transport')
#             # model.addSubstrate(rxnIndx,"Sub1","M_adp_c")
#             model.addProduct(rxnIndx,"Prod1",metID)
#             paramDict = {"K": rate}
#             for key,val in paramDict.items():
#                 model.addParameter(rxnIndx, key, val)

#             SrcSinkRxnsList.append(rxnIndx)
    
#     if AA_sources:     
#         rxnIndx = model.addReaction("Translation","zeroOrderOnOff","translation tRNA src sinks")
#         translat_rate = 0.0003856899829120103 #mM/s
#         model.addSubstrate(rxnIndx,"Sub1","M_alatrna_c",stoich=-23.0)
#         model.addProduct(rxnIndx,"Prod1","M_trnaala_c",stoich=23.0)
#         model.addSubstrate(rxnIndx,"Sub2","M_argtrna_c",stoich=-12.0)
#         model.addProduct(rxnIndx,"Prod2","M_trnaarg_c",stoich=12.0)
#         model.addSubstrate(rxnIndx,"Sub3","M_asntrna_c",stoich=-27.0)
#         model.addProduct(rxnIndx,"Prod3","M_trnaasn_c",stoich=27.0)
#         model.addSubstrate(rxnIndx,"Sub4","M_asptrna_c",stoich=-23.0)
#         model.addProduct(rxnIndx,"Prod4","M_trnaasp_c",stoich=23.0)
#         model.addSubstrate(rxnIndx,"Sub5","M_cystrna_c",stoich=-2.0)
#         model.addProduct(rxnIndx,"Prod5","M_trnacys_c",stoich=2.0)
#         model.addSubstrate(rxnIndx,"Sub6","M_fmettrna_c",stoich=-1.0)
#         model.addSubstrate(rxnIndx,"Sub7","M_glntrna_c",stoich=-15.0)
#         model.addProduct(rxnIndx,"Prod7","M_trnagln_c",stoich=15.0)
#         model.addSubstrate(rxnIndx,"Sub8","M_glutrna_c",stoich=-27.0)
#         model.addProduct(rxnIndx,"Prod8","M_trnaglu_c",stoich=27.0)
#         model.addSubstrate(rxnIndx,"Sub9","M_glytrna_c",stoich=-21.0)
#         model.addProduct(rxnIndx,"Prod9","M_trnagly_c",stoich=21.0)
#         model.addSubstrate(rxnIndx,"Sub10","M_histrna_c",stoich=-5.0)
#         model.addProduct(rxnIndx,"Prod10","M_trnahis_c",stoich=5.0)
#         model.addSubstrate(rxnIndx,"Sub11","M_iletrna_c",stoich=-36.0)
#         model.addProduct(rxnIndx,"Prod11","M_trnaile_c",stoich=36.0)
#         model.addSubstrate(rxnIndx,"Sub12","M_leutrna_c",stoich=-36.0)
#         model.addProduct(rxnIndx,"Prod12","M_trnaleu_c",stoich=36.0)
#         model.addSubstrate(rxnIndx,"Sub13","M_lystrna_c",stoich=-39.0)
#         model.addProduct(rxnIndx,"Prod13","M_trnalys_c",stoich=39.0)
#         model.addSubstrate(rxnIndx,"Sub14","M_mettrna_c",stoich=-7.0)
#         # model.addSubstrate(rxnIndx,"Sub14","M_mettrna_c",stoich=-8.0)
#         model.addProduct(rxnIndx,"Prod14","M_trnamet_c",stoich=8.0)
#         model.addSubstrate(rxnIndx,"Sub15","M_phetrna_c",stoich=-16.0)
#         model.addProduct(rxnIndx,"Prod15","M_trnaphe_c",stoich=16.0)
#         model.addSubstrate(rxnIndx,"Sub16","M_protrna_c",stoich=-11.0)
#         model.addProduct(rxnIndx,"Prod16","M_trnapro_c",stoich=11.0)
#         model.addSubstrate(rxnIndx,"Sub17","M_sertrna_c",stoich=-23.0)
#         model.addProduct(rxnIndx,"Prod17","M_trnaser_c",stoich=23.0)
#         model.addSubstrate(rxnIndx,"Sub18","M_thrtrna_c",stoich=-22.0)
#         model.addProduct(rxnIndx,"Prod18","M_trnathr_c",stoich=22.0)
#         model.addSubstrate(rxnIndx,"Sub19","M_trptrna_c",stoich=-3.0)
#         model.addProduct(rxnIndx,"Prod19","M_trnatrp_c",stoich=3.0)
#         model.addSubstrate(rxnIndx,"Sub20","M_tyrtrna_c",stoich=-13.0)
#         model.addProduct(rxnIndx,"Prod20","M_trnatyr_c",stoich=13.0)
#         model.addSubstrate(rxnIndx,"Sub21","M_valtrna_c",stoich=-24.0)
#         model.addProduct(rxnIndx,"Prod21","M_trnaval_c",stoich=24.0)
#         paramDict = {"K": translat_rate}
#         for key,val in paramDict.items():
#             model.addParameter(rxnIndx, key, val)

#         print("Adding tRNAs")
#         print("rxnIndx ",rxnIndx)
#         SrcSinkRxnsList.append(rxnIndx)    


#     rxnIndx = model.addReaction("Fmet_cofactors","zeroOrderOnOff",'Fmet cofactor source sinks')
#     fmet_rate = 0.00035
# #     model.addSubstrate(rxnIndx,"Sub1","M_thfglu3_c")
#     model.addProduct(rxnIndx,"Prod1","M_10fthfglu3_c")
#     paramDict = {"K": fmet_rate}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)
    
#     SrcSinkRxnsList.append(rxnIndx) 
    
#     rxnIndx = model.addReaction("Fmet_cofactors_2","zeroOrderOnOff",'Fmet cofactor source sinks')
#     fmet_rate_2 = 0.00025
#     model.addSubstrate(rxnIndx,"Sub1","M_thfglu3_c")
# #     model.addProduct(rxnIndx,"Prod1","M_10fthfglu3_c")
#     paramDict = {"K": fmet_rate_2}
#     for key,val in paramDict.items():
#         model.addParameter(rxnIndx, key, val)
    
#     SrcSinkRxnsList.append(rxnIndx) 
    
    
    ### CELL GROWTH VIA SURFACE AREA ###
    SurfaceArea = pmap['CellSA_Prot']
    model.addMetabolite('CellSA_Prot','Prot Surface Area',partTomM(SurfaceArea,pmap))
#     print('CellSA',SurfaceArea)

    SAbool = False # We don't want to grow the SA in this way rn, instead we will do this in the in_out.py stage
    if (SAbool == True):
    #if SurfaceArea<780000:
        
        rxnIndx = model.addReaction("SA_growth","zeroOrderOnOff",'Surface area growth source sinks')
        sa_rate = 0.002
        sa_rate_2 = partTomM(500000,pmap)/105/60
        #     print('SArate',sa_rate,sa_rate_2)
        model.addProduct(rxnIndx,"Prod1","CellSA")
        paramDict = {"K": sa_rate}
        for key,val in paramDict.items():
            model.addParameter(rxnIndx, key, val)

        SrcSinkRxnsList.append(rxnIndx)
    

    #if (SAbool ==False):
        #rxnIndx = model.addReaction("SA_growth_Prot","zeroOrderOnOff",'Surface area growth source sinks')
        #sa_rate = (partTomM(int(round(321329*0.513)),pmap)/(105*60))*0.01
        #model.addProduct(rxnIndx,"Prod1","CellSA_Prot")
        #paramDict = {"K": sa_rate}
        #for key,val in paramDict.items():
            #model.addParameter(rxnIndx, key, val)
        
        #SrcSinkRxnsList.append(rxnIndx)

    #rxnIndx = model.addReaction("gln_transport_supp","zeroOrderOnOff","Adding to gln transport")
    #model.addProduct(rxnIndx,"Prod1","M_gln__L_c")
    #gln_rate = 0.0025 #mM/s
    #paramDict = {"K": gln_rate}
    #for key,val in paramDict.items():
        #model.addParameter(rxnIndx, key, val)
    #model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

    #SrcSinkRxnsList.append(rxnIndx)

    for rxnIndx in SrcSinkRxnsList:
        if "onoff" in model.getReaction(rxnIndx).getKeys():
            model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )

    ## NOT really necessary in CME-ODE format
    SrcSinkRxnsIDsList = [ model.getReaction(rxnIndx).getID() for rxnIndx in SrcSinkRxnsList ]

    return

    
    def buildSinkRxn(substrate, met, rate):
    #     rate_sec = rate*convesionNTP
        # First add the moved metabolite as a new metabolite in the model, for example metabolite metID goes to metabolite metID2
#         model.addMetabolite("{}_sink_c".format(substrate),"{}_Sink".format(substrate),0)

        # add your reaction to the model
        rxnIndx = model.addReaction("{}_SINK".format(substrate),"zeroOrderOnOff","{}_Sink_Rxn".format(substrate))

        # add the metabolite you want to move to the reaction as a substrate
        model.addSubstrate(rxnIndx, "Sub1", met)

        # add the moved metabolite as a new metabolite as a product to the reaction
#         model.addProduct(rxnIndx, "Prod1", "{}_sink_c".format(substrate))

        # add the rate you want to move the metabolite as a parameter in mM/s
        model.addParameter(rxnIndx, "K", rate)
        model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )
    #     sourceSinkRxns.append(rxnIndx)
        return 'Sink reaction for {} ({}) built with a zero order rate of {} mM/s '.format(met,substrate,rate)
    def buildSourceRxn(substrate, met, rate):
    #     rate_sec = rate*convesionNTP
        # First add the moved metabolite as a new metabolite in the model, for example metabolite metID goes to metabolite metID2
    #     model.addMetabolite("{}_source_c".format(substrate),"{}_Source".format(substrate),10)

        # add your reaction to the model
        rxnIndx = model.addReaction("{}_SOURCE".format(substrate),"zeroOrderOnOff","{}_Source_Rxn".format(substrate))

        # add the metabolite you want to move to the reaction as a substrate
    #     model.addSubstrate(rxnIndx, "Sub1", )

        # add the moved metabolite as a new metabolite as a product to the reaction
        model.addProduct(rxnIndx, "Prod1", met)

        # add the rate you want to move the metabolite as a parameter in mM/s
        model.addParameter(rxnIndx, "K", rate)
        model.addParameter(rxnIndx, 'onoff',1,lb=0, ub=1, unit="mM", parName="Debug On/Off switch" )
    #     sourceSinkRxns.append(rxnIndx)
#         return 'Source reaction for {} ({}) built with a zero order rate of {} mM/s '.format(met,substrate,rate)
        print('Source reaction for {} ({}) built with a zero order rate of {} mM/s '.format(met,substrate,rate))

    gDW = 10.2e-15
    r = 2.0*(10**-7) #meters 
    volume = 4.0 / 3.0 * r**3 * np.pi* 1000 #Liters
    #FBA in mmol/(gDW * h) to  Sim in mM/s
    convFBAtoSim = gDW / 3600 / volume

    # sourceSink = pd.read_csv('../../model_data/nucleo_sourceSink_params.csv',header=0)
#     sourceSink = pd.read_csv('../../../model_data/TB_nuc/nucleo_sourceSink_params_noCentral.csv',header=0)
    sourceSink = pd.read_csv('../model_data/nucleo_sourceSink_params_noLip_noCentral_noAA.csv',header=0)
#     sourceSink.head()

    for ind,row in sourceSink.iterrows(): 
        if row[4] and row[3]=='Sink':
            buildSinkRxn(row[0],row[1],row[2]*convFBAtoSim)
            sourceSinkRxns.append('{}_SINK'.format(row[0]))
    #         print('sink')
        if row[4] and row[3]=='Source':
    #         print('source')
            buildSourceRxn(row[0],row[1],row[2]*convFBAtoSim)
            sourceSinkRxns.append('{}_SOURCE'.format(row[0]))
#     sourceSink
    # sourceSink[sourceSink['Sink / Source']=='Sink']
#     return model

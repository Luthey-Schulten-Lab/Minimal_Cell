Author: David Bianchi
Date: 9/29/21

A Description of the model data input files for multi-scale CME-ODE Simulations of *in-silico* JCVI-syn3A cells is given below:

Metabolic Module Parameter Files
=============================================

Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv - Parameter File for Central and Amino Acid Metabolism Modules, with parameters resulting from parameter balancing and literature parameter value collation process.

lipid_NoH2O_balanced_model.tsv - Parameter File for Lipid Metabolism Module, with parameters resulting from parameter balancing and literature parameter value collation process.

Nucleotide_Kinetic_Parameters.tsv - Parameter File for Nucleotide Metabolism Module, with parameters resulting from parameter balancing and literature parameter value collation process.

GlobalParameters_Zane-TB-DB.csv - Global Simulation Parameters including external growth medium metabolite concentrations *etc.*

transport_NoH2O_Zane-TB-DB.tsv - Parameter File for Cellular Transport Reactions which are handled via different rate laws that intracellular metabolic reactions and result in a signficant energy expense in JCVI-syn3A relative to other bacterial cells.



Gene Expression Condition Files
=============================================

RNApol_proteins.csv -

mRNA_counts.csv -

proteomics.xlsx -

riboPtnInfo.csv -

rrna_metabolites_1.csv -

rrna_metabolites_2.csv -

trna_metabolites.csv -

trna_metabolites_synthase.csv -




Reaction IDs Files
=============================================

lipid_rxns_list.txt -

nucleo_rxns_list.txt -




Metabolic Model and Genome Files
=============================================

syn2.gb -

syn3A.gb -

manual_GPR_conversion.csv -


./FBA
==================================================


Syn3A_annotation_compilation.xlsx  -

iMB155.json -




Miscellaneous Files
==================================================

membrane_protein_metabolites.csv -

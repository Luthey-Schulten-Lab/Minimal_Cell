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

RNApol_proteins.csv - List of the RNA polymerase gene loci present in JCVI-syn3A. The protein counts of these genes are used for the molecular machine that "catalyzes" transcription reactions in the gene expression module.

mRNA_counts.csv - Average values of mRNA counts for each gene in JCVI-syn3A resulting from stochastic simulations. These averages are used to generate Poisson Distributions of mRNA counts for each, from which initial mRNA counts for each gene are sample in multi-scale CME-ODE/RDME-gCME-ODE simulations of JCVI-syn3A cells. 

proteomics.xlsx - The proteomics values for JCVI-syn3A genes reported in *Breuer et al eLife (2019)*. These values are used as initial conditions for simulations.

rrna_metabolites_1.csv - List of the first rRNA operon present in JCVI-syn3A, used in the transcription module of JCVI-syn3A stochastic gene expression simulations.

rrna_metabolites_2.csv - List of the second rRNA operon present in JCVI-syn3A, used in the transcription module of JCVI-syn3A stochastic gene expression simulations.

trna_metabolites.csv - List of the tRNA gene loci present in JCVI-syn3A, used in the stochastic gene expression (transcription/translation) modules of JCVI-syn3A simulations.

trna_metabolites_synthase.csv -List of the tRNA synthase gene loci present in JCVI-syn3A, used in the stochastic gene expression (transcription/translation) modules of JCVI-syn3A simulations.




Reaction IDs Files
=============================================

lipid_rxns_list.txt - List of non-uptake (*e.g.* sphingomyelin) lipid metabolic module reaction IDs used in JCVI-syn3A simulations.

nucleo_rxns_list.txt - List of non-uptake (*e.g.* nucleoside uptake) nucleotide metabolic module reaction IDs used in JCVI-syn3A simulations.




Metabolic Model and Genome Files
=============================================

syn2.gb - The JCVI-syn2.0 genome file, used to connect IDs reported for proteomics in (Breuer et al eLife 2019) to the JCVI-syn3A gene locus tag IDs used in these simulations.

syn3A.gb - The JCVI-syn3A genome file, used in simulation species (genes, RNA, *etc.) definition in JCVI-syn3A simulations.

manual_GPR_conversion.csv - A file used for outlier genes (*.i.e* those without reported Reaction to Gene Locus tag definitions in the metabolic module presented in Breuer et al eLife 2019), to connect IDs reported for proteomics in (Breuer et al eLife 2019) to the JCVI-syn3A gene locus tag IDs used in these simulations.


./FBA
==================================================


Syn3A_annotation_compilation.xlsx  - File containing annotations and protein product descriptions for all genes in JCVI-syn3A (modified from the version given in Breuer et al eLife 2019). 

iMB155.json - The json metabolic model file for the JCVI-syn3A model presented in Breuer et al eLife 2019.




Miscellaneous Files
==================================================

membrane_protein_metabolites.csv - A list of all proteins in JCVI-syn3A predicted to have transmembrane helix components and that will contribute to growing the cell membrane surface in the *in silico* growth model. Initial conditions for genes with different states (*i.e.* the phosphorylation states of PtsG) are given. 

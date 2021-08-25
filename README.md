# coevolution-compatibility
Companion software for "Coevolutionary methods enable robust design of modular repressors by reestablishing intra-protein interactions"

AUTHORS
Xian-Li Jiang, Rey P. Dimas, Clement T. Y. Chan, Faruck Morcos

Description: This README file describes the codes and data accompanying the above publication.

Files:

DCAparameters.m: MATLAB code, used to extract coevolutionary parameters, eij and hi, from homologous multiple sequence alignments (MSAs) of interested protein family. 

CompatibilityScore.m: MATLAB code, compute compatibility score for any aligned sequence, using parameters obtained from DCAparameter.m.

StructueFitness.m: MATLAB code, compute structural fitness score for any aligned sequence using parameters obtained from DCAparameter.m.

MSA_LacI_homologs.fasta: Fasta file, contains homologous multiple sequence alignments (MSAs) of protein members from a family (LacI family for this paper). It is used as input file for coevolutionary analysis. 

LacI family residue pairs with top 1500 DI values: Txt file used for compatibility score computation. It contains strongly coevolved residue pairs between DNA recognition module and environmental sensing module computed from same MSAs file used for DCAparameter.m code. Their direct information values are computed by DCA method (dca.m) and are among the top 1500. The original dca code can be found in dca.rice.edu.

1lbg_monomer.txt: Txt file used for structural fitness score computation. It contains contacting residue pairs obtained from LacI crystal structure (PDB:1lbg).

Aligned sequence of hybrid mutants.fasta: Fasta file, contains aligned sequences for 8 hybrid repressors and 11 mutants (selected by this paper) for each hybrid. The 8 original aligned sequence can be used to generate any mutants of interest, such as all possible single point or double point mutations.


System requirements and Installation:
All the codes and MATLAB data file are ran and open by the software MATLAB (tested on MATLAB R2019b and newer), which can be downloaded from: https://www.mathworks.com/products/matlab.html.
Bioinformatics Toolbox needs to be installed in MATLAB.


Steps:


1. MSAs and DCAparameters.m:
	Generate coevolutionary parameters (familycouplings and local fields, denoted as eij and hi in the paper) using the DCAparameters.m code. 
	Example command: [lf,eij]=DCAparameters('MSA_LacI_homologs.fasta',1)
	The output variables lf and eij can be save to mat file for future use.

2. MSAs and dca.m:
	Calculate the pairwise direct information (DI) values by using the dca.m code in the following link: http://dca.rice.edu/portal/dca/
	Example command: dca('MSA_LacI_homologs.fasta','DI_values_each_residue_pairs.DI')
	Sort the pairs with inclusion of only inter-module pairs, select top 1500 pairs.
	Alternatively, you can use the file 'LacI family residue pairs with top 1500 DI values.txt' directly.

3. Calculate the compatibility score by CompatibilityScore.m: 
	LacI family residue pairs with top 1500 DI values.txt (provided or calculated from step 2)
	Coevolution parameters (eij and lf in this demo) calculated from step 1
	Aligned sequence of hybrid mutants.fasta

	Example command: CS=CompatibilityScore('Aligned sequence of hybrid mutants.fasta',eij,lf,1)
	Expected output: An array of CS scores for all the sequences in the fasta file.

4. Calculate the structural fitness score by StructueFitness.m: 
	1lbg_monomer.txt (provided)
	Coevolution parameters (eij and lf in the demo) calculated from step 1
	Aligned sequence of hybrid mutants.fasta

	Example command: SF=StructureFitness('Aligned sequence of hybrid mutants.fasta',eij,lf,1)
	Expected output: An array of SF scores for all the sequences in the fasta file.

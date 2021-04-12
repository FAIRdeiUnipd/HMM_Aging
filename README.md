Repository for matlab data and codes to generate the figures from the paper "Unveil whole-brain dynamics in normal aging through Hidden Markov Models" by Moretto et al.

Code:

	Code/MAIN_code_Moretto_NeuroImage.m
	main code to run the analyses and generate the figures of the paper

Results folder is organized in subfolders, each of which describes  an analysis

	Results/Chronnectome/Metrics.mat
	Results/Cvs/CV.mat
	Results/Graph_metrics/Graph_metrics.mat
	Results/HMM_6states_rep5_conf2/fe_rep.mat
	Results/HMM_6states_rep5_conf2/HMM_rep_5_conf_2.mat
	Results/HMM_6states_rep5_conf2/option_conf_2.mat
	Results/Modularity/Modularity_50000rep.mat
	Results/Modularity/Results_modularity_50000rep.mat
	names_IC_RSNs.mat
	NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz

Figures folder contains the figures of the Manuscript and of Supplementary Material

Required toolboxes/utilities:

	HMM-MAR: https://github.com/OHBA-analysis/HMM-MAR
	BCT Version 2019-03-03:	https://sites.google.com/site/bctnet/
Matlab Functions:

	NIfTI toolbox: https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
	fdr_bh: https://it.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh

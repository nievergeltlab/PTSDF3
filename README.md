# PTSDF3
This is the public code base for PGC-PTSD Freeze 3 (Nat Genet, Accepted). Each folder contains the code for a given analysis. Details of operation can be found in code.
  
Meta-analysis summary statistics will be made available upon publication here: https://pgc.unc.edu/for-researchers/download-results/. 

Access to study level summary statics and individual level data is facilitated here: https://pgc.unc.edu/for-researchers/data-access-committee/data-access-portal/


###Notes on the meta-analysis
•	The GWAS meta-analysis included studies where PTSD was measured as either case/control or PTSD symptoms
o	The beta estimates from these analysis have distinct meanings and are therefore on different measure scales
	i.e. log odds of PTSD diagnosis / increase in symptom scores
	thus not directly combinable in traditional meta-analysis, which creates an average beta value
•	Instead,  we performed a sample-size weighted, z-score based meta-analysis. 
o	This does not require the assumption of the same measure scale.
o	Sample size weighted meta-analysis however does not produce beta values - a beta generated from disparate measure types would not have any direct meaning 
•	Many post-GWAS still require beta values
o	We propose to generate these pseudo-beta values based on the formalue provided to you
o	These betas provide
	a relative difference in effect sizes between SNPs - i.e. that larger betas mean relatively stronger effects than smaller betas
	effect directions are preserved - positive is risk increasing, negative is risk decreasing
o	Because relative effect size and direction are preserved, have been adequate for most analyses
	In Mendelian randomization, however, this gives no direct quantification of causal effect magnitude - significant MR can only be interpreted as 'risk increasing' or 'risk decreasing'



### COJO
Conditional analysis of SNPs within significant genes (from gene-set analyses)

### Finemapping
Fine-mapping of risk loci using Polyfun

### Forest plots
Forest plots for significant loci

### GO term plot
Plot significant GO term results

### GWAS and Meta-analysis
Format phenotype data for GWAS  
Perform GWAS of datasets using PLINK and BOLT (where applicable)  
Reformat summary statistics for METAL  
Meta analysis in METAL  

### Gene-mapping Venn Diagram
For supplementary figure of # of genes mapped by each gene-mapping method

### Gene-set analyses
Conduct gene-set analysis of brain datasets

### LDSC
Heritability analyses

### Local rg
LAVA local heritability

### Locuszoom
Regional association plots

### male-female rg
Genetic correlations between men and women

### Manhattan plots
Generate Manhattan plots 

### MiXeR
Mixer (univ and biv) analyses

### PRS
Polygenic risk scoring plots

### Prior hit lookup
Determine if hits overlap previously reported hits

### Prioritization
Gene prioritization analyses

### Study Bubble PLots
Make bubble plots from figure 1

### rg with psychiatric
Calculate rg with other psychiatric disorders



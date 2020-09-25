# TFA-evaluation inputs from Ma&Brent 2020

These files are adapted from the Supplemental Files of Ma&Brent 2020 to show replicate results shown in Fig2A

Example run:

python evaluateTFA.py -a Paper_Inputs/S5_unlabeled_50-TF_TFA_values_optimized_with_ChIP-CC-TFKO_on_ZEV.csv -m Paper_Inputs/S5_TF_mRNA_in_ZEV.csv -t Paper_Inputs/S5_TF_labels.csv -s Paper_Inputs/sample_labels.csv -p 2

Expected output:

    TFA values will be evaluated according to file of expected activity pattern
    Median rank percentile of perturbed TFs:
	    86.0 	p-value:  1.55e-09
    Percent of TFA changes in expected direction:
	    0.66 	p-value:  1.64e-02
    Percent positive correlation between TFA and mRNA:
	    0.74 	p-value:  4.68e-05
    Percent positive correlation between TFA and mRNA (bootstrapped median):
	    0.72


-a

    a csv file of TF activity values, arranged as a TF x samples matrix
    these TFA values were inferred by using the ChIP-CC network, sign constrained sign constrained and optimized with TFKO dataset
 
-m

    a csv file of TF gene expression values from the ZEV dataset that matches the arrangement of the TFA file
    this file is used to calculate the positive correlation metric
    the percent of TFs with a positive correlation is calculated 
      p-value is calculated as probability to get that percentage with # coin flips = # TFs
    the median percentage of TFs with a positive correlation across 1000 bootstrap samples is calculated

-t

    a file that labels the rows of the activity matrix with TF names
    
-s

    a file that labels the columns of the activity matrix
    if a sample has an expected change in TFA, it should be named with that TF to match the TF file
    wherever the sample and tf labels match, the TFA value will be used for calculating the correct direction and median rank metrics
    
-p

    the direction of TFA perturbation
    0 indicates low activity
    2 indicates high activity
    

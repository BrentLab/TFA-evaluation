# TFA-evaluation

Code that evaluates inferred TFA values according to three metrics described in Ma&Brent 2020

1. Fraction of TFs who are predicted to change activity in the *correct direction*
2. *Median rank* percentile of perturbed TF in each sample with a known perturbation
3. Fraction of TFs whose inferred TFA values and gene expression levels have *positive correlation*

Package requirements to run:

1. python 2.7 or python 3.6
2. numpy 1.18
3. scipy 1.4

Example run 1:

python evaluateTFA.py -a Inputs/TFA_inferred_from_ZEV.csv -m Inputs/measured_TF_mRNA_in_ZEV.csv -e Inputs/ChIP_ZEV_TFA_expectations_180samples.csv

Expected output:

    TFA values will be evaluated according to file of expected activity pattern
    Median rank percentile of perturbed TFs:
	    86.0 	p-value:  1.55e-09
    Percent of TFA changes in expected direction:
	    0.78 	p-value:  4.51e-05
    Percent positive correlation between TFA and mRNA:
	    0.78 	p-value:  4.51e-05
    Percent positive correlation between TFA and mRNA (bootstrapped median):
	    0.72


-a

    a csv file of TF activity values, arranged as a TF x samples matrix
    this is required input
-m

    a csv file of TF gene expression values that matches the arrangement of the TFA file
    this file is used to calculate the positive correlation metric
    the percent of TFs with a positive correlation is calculated 
      p-value is calculated as probability to get that percentage with # coin flips = # TFs
    the median percentage of TFs with a positive correlation across 1000 bootstrap samples is calculated
-e

    a csv file of TF activity expectations that matches the arrangement of the TFA file
    this file is used to calculate the median rank and correct direction metrics
    1 indicates no expectations
    0 indicates TFA is expected to be low
    2 indicates TFA is expected to be high
    the correct direction metric compares all 0/2 TFA values against the TF's activity in the reference sample
      reference sample is the last column by default, but can be chosen with the -w flag
      p-value is calculated as probability to get that percentage with # coin flips = # TFs
    the median rank percentile metric compares the z-scores all 0/2 TFA values against the z-scores of all other TFs in the same sample
      p-value is calculated as probability to get half the TFs above the median rank percentile given the random chance of achieving that percentile (1-percentile)
    
Example run 2:

python evaluateTFA.py -a Inputs/TFA_inferred_from_ZEV.csv -m Inputs/measured_TF_mRNA_in_ZEV.csv -t Inputs/tf_list.csv -s Inputs/sample_list.csv -w -1 -p 2 -b 1

Expected output:

    TFA values will be evaluated according to sample labeling of perturbed TFs
    Perturbation direction is increased activity
    Median rank percentile of perturbed TFs:
	    86.0 	p-value:  1.55e-09
    Percent of TFA changes in expected direction:
	    0.78 	p-value:  4.51e-05
    [====================]100% creating null distribution
    Percent positive correlation between TFA and mRNA:
	    0.78 	p-value:  4.51e-05
    Percent positive correlation between TFA and mRNA (bootstrapped median):
	    0.72 	p-value:  0.00e+00


-a

    a csv file of TF activity values, arranged as a TF x samples matrix
    this is required input
-m

    a csv file of TF gene expression values that matches the arrangement of the TFA file
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
    
-w

    the column index of the reference or WT sample
    default is -1 for the last column
    
-p

    the direction of TFA perturbation
    0 indicates low activity
    2 indicates high activity
    
-b

    yes (1) or no (0) to calculate the significance of the median positive correlation across 1000 bootstraps
    default is 0 since it takes up to an hour to run
    this will create a null distribution of 1000 results from randomly pairing activity levels of one TF with the gene expression levels of another


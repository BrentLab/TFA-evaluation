# Network Construction

Code and inputs needed to reconstruct the networks in Ma&Brent 2020
Special thanks for coding contribution from Lisa Liao

Package requirements to run:

1. python 3.6
2. numpy 1.18
3. scipy 1.4
4. pandas 1.0

Included input files:

cisBP.csv

    list of yeast TFs pulled from the cisBP database
    used to identify samples as TF perturbation samples

YeastGeneNames.csv

    list of yeast genes, one per row, with systematic name as first column and common names added as needed
    names were pulled from SGD: yeastgenome.org

harbison_ChIP_score_matrix_trimmed.csv

    a genes x TFs matrix of ChIP binding scores from Harbison, 2004

PWM_max_pVal_matrix_trimmed.csv

    a genes x TFs matrix of max -log10 p-values for a TF's binding site being found in the gene's promoter
    TF binding motifs were taken from the ScerTF database: stormo.wustl.edu/ScerTF

kemmerenKOexpressionMatrixTrimmed.csv

    a genes x samples matrix of log2FC gene expression data measured using microarray
    this is a trimmed version of the TFKO dataset published by the Holstege Lab:
    http://deleteome.holstegelab.nl/data/downloads/deleteome_all_mutants_controls.txt
    only includes samples of perturbed genes listed by cisBP as a TF gene and have a perturbation sample in the ZEV dataset
    a WT sample of 0 gene expression values is included at the end to represent the theoretical reference sample for all the experimental sample

zev15minExpressionMatrixTrimmed.csv

    a genes x samples matrix of log2FC gene expression data measured using microarray
    this is a trimmed version of the ZEV dataset published by the McIsaac Lab:
    https://storage.googleapis.com/calico-website-pin-public-bucket/datasets/pin_tall_expression_data.zip
    only includes 15min timepoint samples of perturbed genes listed by cisBP as a TF gene and have a perturbation sample in the ZEV dataset
    a WT sample of 0 gene expression values is included at the end to represent the theoretical reference sample for all the experimental sample
    

Example run 1:

python construct_network.py -n 50 -d 0 -e input/harbison_ChIP_score_matrix_trimmed.csv -lr 1250 -mt 2

-n

    the number of TFs desired in network
    this is required input

-d

    the direction of perturbation for TFs in samples that are labeled with that TF
    0 indicates knock-out or knock-down
    2 indicates induced or over-expression

-e

    a csv file of scores arranged as a genes x TFs scores for likelihood of a TF-target gene edge
    high, absolute values will be treated as high scores/support for a TF-target gene edge
    this is required input

-lr

    an integer value indicating the lowest ranking score to consider adding to the network

-mt

    an integer value indicating the minimum number of target genes per TF


Expected printed output:

    Edge score within rank 1250 : 4.37220825205
    Starting network construction...
    Found 50 TFs. Highest rank 419 , Score: 6.7254219200900005
    Removing unidentifiability...
    49 TFs in network
    Found 50 TFs. Highest rank 450 , Score: 6.5479757651199995
    Removing unidentifiability...
    50 TFs in network
    Total number of edges: 1104
    Final network: 50 TFs, 778 genes
    End network construction.
    Writing network to files...
    Done!
    Writing network to files...
    Done!

Expected file output:

hk/hk_binary_cs_corr_signed.csv
hk/hk_binary_cs_signed.csv

    the binary ChIP-CC and ChIP-PC networks that use the TFKO dataset to define sign constraints

hk/hk_binary_tfa_all.csv
hk/hk_binary_tfa_perturbed.csv
hk/hk_binary_tfa_unperturbed.csv

    the binary TFA matrices with a 0 value to indicate low activity expected for the TF in its perturbed TFKO sample
    perturbed is the subset of samples that perturb TFs in the network
    unperturbed is the subset of samples that perturb TFs not in the network
    a column of 1's is added to all files as the reference/wildtype sample where TFs are not perturbed

hk/hk_expression_all.csv
hk/hk_expression_perturbed.csv
hk/hk_expression_unperturbed.csv

    the expression datasets subsetted from TFKO to use for optimization of TFA and CS values
    perturbed is the subset of samples that perturb TFs in the network
    unperturbed is the subset of samples that perturb TFs not in the network
    a column of 0's is added to all files as the reference/wildtype sample

hk/hk_gene_list.csv
hk/hk_tf_list.csv
hk/hk_sample_list_all.csv
hk/hk_sample_list_perturbed.csv
hk/hk_sample_list_unperturbed.csv

    the gene, tf, and sample labels for the cs, tfa, and expression matrices

zev/zev_binary_cs_corr_signed.csv
zev/zev_binary_cs_signed.csv

    the binary ChIP-CC and ChIP-PC networks that use the ZEV dataset to define sign constraints

zev/zev_binary_tfa_all.csv
zev/zev_binary_tfa_perturbed.csv
zev/zev_binary_tfa_unperturbed.csv

    the binary TFA matrices with a 2 value to indicate high activity expected for the TF in its perturbed ZEV sample
    perturbed is the subset of samples that perturb TFs in the network
    unperturbed is the subset of samples that perturb TFs not in the network
    a column of 1's is added to all files as the reference/wildtype sample where TFs are not perturbed

zev/zev_expression_all.csv
zev/zev_expression_perturbed.csv
zev/zev_expression_unperturbed.csv

    the expression datasets subsetted from ZEV to use for optimization of TFA and CS values
    perturbed is the subset of samples that perturb TFs in the network
    unperturbed is the subset of samples that perturb TFs not in the network
    a column of 0's is added to all files as the reference/wildtype sample

zev/zev_gene_list.csv
zev/zev_tf_list.csv
zev/zev_sample_list_all.csv
zev/zev_sample_list_perturbed.csv
zev/zev_sample_list_unperturbed.csv

    the gene, tf, and sample labels for the cs, tfa, and expression matrices


Example run 2:





Edge score within rank 1250 : 1.3135128999999999
Starting network construction...
Found 50 TFs. Highest rank 567 , Score: 1.6255416999999999
Removing unidentifiability...
29 TFs in network
Cannot find 50 TFs with rank threshold of 1250


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

Calls for networks evaluated in Fig2A of Ma&Brent2020:

ChIP-CC and ChIP-PC using ZEV data as prior knowledge

    python construct_network.py -n 50 -d 2 -e input/harbison_ChIP_score_matrix_trimmed.csv -lr 1250 -mt 2

DE-PC using ZEV data as prior knowledge
    python construct_network.py -n 50 -d 2 -e input/harbison_ChIP_score_matrix_trimmed.csv -lr 1250 -mt 2









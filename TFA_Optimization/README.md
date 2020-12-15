# TFA Inference

Code that infers TFA values as described in Ma&Brent 2020

Package requirements to run:

1. gurobi
2. python 2.7 or python 3.6
3. numpy

Files in Paper datasets directory:

glucoseInflux_0p02percent_0p2percent_expressionData.csv

    a csv file of gene expression data published by Ronen and Botstein, 2006
    the first series, labeled by timepoints, is gene expression data after glucose influx of 0.02%
    the second series, labeled by timepoints, is gene expression data after glucose influx of 0.2%
    duplicate gene expression measurements were removed by dropping the row(s) with less total abs value

glucoseInflux_2percent_expressionData.csv

    a csv file of gene expression data published by Apweiler et al, 2012
    labeled by timepoints, samples are gene expression data after glucose influx of 2%

TFdoubleKOexpressionData.csv

    a csv file of gene expression data published by Sameith et al, 2015
    gene expression data of steady state cells with two TF-encoding genes deleted

zamanExpressionDataRescaled.csv
zamanSampleKey.csv

    a csv file of gene expression data published by Zaman et al, 2009
    includes time series of various conditions and genotypes, with a key file for looking up sample labels

zev15minExpressionData.csv
zevAllExpressionData.csv

    csv files of gene expression data published by Hackett et al, 2020
    includes multiple inductions of TF-encoding genes
	the first file with 15min timepoint
	the second file is grouped by timepoints 2.5, 5, 10, 20, 30, 45, 60, and 90min


Other notes:

1. code assumes the existence of subdirectories logFiles and modelFiles for logging results during optimization
2. Gurobi offers free academic licenses for their optimizer at gurobi.com > Academia > Academic Program and Licenses
3. data in Paper datasets directory can be subsetted to run TFA inference as reported in the paper

Example run:

python TFAinference.py -e hkMatrixE_optimize.csv -c hkSignedBinaryCS.csv -i 20 -v hkMatrixE_holdOut.csv -t Example

-e

    a csv file of target gene expression data, arranged as a gene x samples matrix
    this is required input

-c

    a csv file of TF-gene edges, arranged as a TF x gene matrix
    0  indicates no edge
    1  indicates an activating edge
    -1 indicates a repressing edge
    this is required input

-i

    max number of iterations before stopping
    default is 100

-v

    a csv file of target gene expression data
    current CS values will be used at each iteration to fit TFA values to this set and the variance explained calculated
    the variance explained on this second expression dataset can be used to find a better stopping point than an arbitrary number of iterations

-t

    a string to tag the output files with, useful if doing multiple random starts


Expected printed output:

**the numbers won't be exactly the same because a new random set of starting CS values is generated with each execution**
  
    Variance explained and error from random CS values:  -2369.709 3293052.697 
    cross check:  0.0 -496.659 3194524.619

    iteration  1 

    Variance explained and error at iteration  1 :  0.526 658.662
    cross check:  1.0 0.305 4458.43

    iteration  2 

    Variance explained and error at iteration  2 :  0.555 618.007
    cross check:  0.99 0.333 4282.723

    iteration  3 
    ...
    ...

    
Expected file outputs:

logFiles/ExampleCSiteration0.csv - ExampleCSiteration19.csv
logFiles/ExampleTFAiteration0.csv - ExampleTFAiteration19.csv

    these are the CS/TFA values inferred across 20 iterations

logFiles/ExampleCS.log
logFiles/ExampleTFA.log

    this is the gurobi output when optimizing CS values / TFA values / TFA values for the optional -v input
    logged output from each iteration is appended to the end of this file

modelFiles/ExampleCS.lp
modelFiles/ExampleTFA.lp
modelFiles/ExampleCrossCheck.lp

    this is the model file created for gurobi to optimize CS values / TFA values / TFA values for the optional -v input
    at each iteration, the entire file is updated to match the current model



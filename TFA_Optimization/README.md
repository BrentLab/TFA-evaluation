# TFA Inference

Code that infers TFA values as described in Ma&Brent 2020

Package requirements to run:

1. gurobi
2. python 2.7 or python 3.6
3. numpy

Other notes:

1. code assumes the existence of subdirectories logFiles and modelFiles for logging results during optimization
2. Gurobi offers free academic licenses for their optimizer at gurobi.com > Academia > Academic Program and Licenses

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



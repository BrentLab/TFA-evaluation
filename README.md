# TFA-evaluation

Code that evaluates inferred TFA values according to three metrics described in

1. Fraction of TFs who are predicted to change activity in the correct direction
2. Median rank percentile of perturbed TF in each sample with a known perturbation
3. Fraction of TFs whose inferred TFA values are positively correlated with their mRNA levels

Example run:

python evaluateTFA.py -a <inferred TFA matrix> -e < "binary" TFA matrix of expected behavior>  -m <TF mRNA expression data matrix> 

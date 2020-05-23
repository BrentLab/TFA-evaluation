"""
Evaluated inferred TFA values
"""

import argparse
import numpy
import scipy.stats
import math
import random

"""
path to a CSV file of TFA values is required input [-a]
(optional) column index of reference/WT, default: -1 [-w]

the values can then be evaluated based on two sets of additional inputs
1. 
	path to a CSV file of activity values expectated relative to reference/WT [-e]
		should have same dimensions as TFA input
		this path is required if checking for a mix of perturbation directions
                        e.g. if some samples should have increased activity while others should have decreased
	
2.
	path to a file of TF labels for the rows of TFA values [-t]
	path to a file of sample labels for the columns of TFA values [-s]
		sample labels should be the TFs expected to be perturbed
	(optional) direction of all perturbations relative to reference/WT, default: 3 [-p]
	
"""
def parseInputs():
	continueFlag = False

	parser = argparse.ArgumentParser(description='Evaluate TFA values')
	parser.add_argument('--activity', '-a', type=str, 
					dest='tfaFile', action='store', required=True,
                	help='path to csv file of activity matrix')
	parser.add_argument('--WT', '-w', '--wt', type=int, 
					dest='wtIndex', action='store', default=-1,
                	help='column index in activity that holds the WT or reference sample')
	parser.add_argument('--expectation', '-e', type=str, 
					dest='perturbFile', action='store',
                	help='csv file of expected activity behavior relative to WT\
                	\n -1: decrease\t 0: unknown\t 1: no expectation\t 2: increase')
	parser.add_argument('--TFs', '-t', '--tfs', type=str, 
					dest='tfList', action='store',
                	help='path to file of TFs labeling the rows of activity input')
	parser.add_argument('--samples', '-s', type=str, 
					dest='sampleList', action='store',
                	help='file of TFs labeling the columns of activity input')
	parser.add_argument('--perturbation', '-p', type=int, 
					dest='perturbFlag', action='store',
                	default=3, choices = [0, 2, 3],
                	help='direction of expected activity behavior relative to WT\
                	\n 0: decrease\t 2: increase\t: 3: unknown')
	parser.add_argument('--mRNA', '-m', '--mrna', type=str,
					dest='mrnaFile', action='store',
			help='path to csv file of measured TF mRNA')


	args = parser.parse_args()
	#fileInputs = tfaValues, expectedValues, tfList, sampleList, mrnaValues
	fileInputs = [0, 0, 0, 0, 0]

	#check for required TFA data
	try:
		tfaValues = numpy.genfromtxt(args.tfaFile, delimiter=',')
		fileInputs[0]=tfaValues
	except IOError as e:
		print("Error:\t", e)
		quit()

	#check for data needed to go route 1
	if args.perturbFile != None:
		try:
			expectedValues = numpy.genfromtxt(args.perturbFile, delimiter=',')
			if expectedValues.shape == tfaValues.shape:
				continueFlag = True
				fileInputs[1]=expectedValues
				print("TFA values will be evaluated according to file of expected activity pattern")
			else:
				print("Dimensions of activity values ", tfaValues.shape, 
					" do not match expected activity behavior ", expectedValues.shape,
					"\n Will attempt to use alternate method of evaluation")
		except IOError as e:
			print("Error:\t", e, "\n Will attempt to use alternate method of evaluation")

	#check for data needed to go route 2
	if not continueFlag and args.tfList != None and args.sampleList != None:
		try:
			with open(args.tfList, 'r') as tfFile:
				tfs = [x.strip('"') for x in tfFile.read().split()]
			with open(args.sampleList, 'r') as sampleFile:
				samples = [x.strip('"').split('_') for x in sampleFile.read().split()]
			if (len(tfs),len(samples))==tfaValues.shape:
                                continueFlag = True
                                fileInputs[2]=tfs
                                fileInputs[3]=samples
                                print("TFA values will be evaluated according to sample labeling of perturbed TFs")
                                if args.perturbFlag == 3:
                                        print("Perturbation direction is unknown")
                                elif args.perturbFlag == 0:
                                        print("Perturbation direction is decreased activity")
                                else:
                                        print("Perturbation direction is increased activity")
			else:
				print("Dimensions of activity values ", tfaValues.shape,
					" do not match tf and sample labels ", (len(tfs),len(samples)))
		except IOError as e:
			print("Error:\t", e)

	#check for mRNA data
	if args.mrnaFile != None:
                try:
                        mrnaValues = numpy.genfromtxt(args.mrnaFile, delimiter=',')
                        if tfaValues.shape==mrnaValues.shape:
                                fileInputs[4] = mrnaValues
                        else:
                                print("Dimensions of activity values ", tfaValues.shape,
					" do not match mRNA values ", mrnaValues.shape)
                                print("Correlation metric will not be calculated")
                except IOError as e:
                        print("Error:\t", e)
                        

	if not continueFlag:
		print("\nCan not evaluate TFA values with given inputs.")
		print("User must give either a csv file of expected activity relative to the reference/WT sample")
		print("\twhich matches the dimensions of the activity values [-e]")
		print("OR files that label the activity values with TFs for the rows [-t]")
		print("\tand label the activity values with expected perturbed TFs for the columns [-s]")
		quit()


	return [args, fileInputs]


"""
returns a matrix of tf x sample values
if the sample label includes a tf, that value should be perturbDirection
else it should be 1
options for perturb Direction should be 0, 2, 3
"""
def makeExpectationMatrix(tfLabels, sampleLabels, perturbDirection):
	expectedValues = numpy.empty([len(tfLabels), len(sampleLabels)])
	for i in range(len(tfLabels)):
		for j in range(len(sampleLabels)):
			if tfLabels[i] in sampleLabels[j]:
				expectedValues[i][j]=perturbDirection
			else:
				expectedValues[i][j]=1
	return expectedValues


"""
checks tfaValues for expected directed behavior
for each 0 or 2 value in expectedValues, check that index in tfaValues
	0: should be less than value at same row index, column wtIndex
	2: should be greater than value at same row index, column wtIndex
"""
def countDirected(tfaValues, expectedValues, wtIndex):
	numComparisons = 0
	numCorrect = 0
	for i in range(len(tfaValues)):
		for j in range(len(tfaValues[0])):
			if expectedValues[i][j] == 0:
				numComparisons += 1
				if tfaValues[i][wtIndex] > tfaValues[i][j]:
					numCorrect += 1
			elif expectedValues[i][j] == 2:
				numComparisons += 1
				if tfaValues[i][wtIndex] < tfaValues[i][j]:
					numCorrect += 1
	return [numComparisons, numCorrect]


"""
uses expected values to identify TFs in each sample that are expected to rank highly

1: convert to z-scores or z-scores after log if needed
2: rank expected TFs compared to TFs w/o expectations
3: return list of rank percentiles
"""
def rankTFs(tfaValues, expectedValues):
        #convert tfaValues to z-scores of logged tfaValues
        standardizedTFAmatrix = tfaValues.copy()
        for i in range(len(tfaValues)):
                logTFA=[math.log(x+0.0001,2) for x in tfaValues[i]]
                standardizedTFAmatrix[i] = [(x-numpy.average(logTFA))/numpy.std(logTFA)
                                            for x in logTFA]
        #calculate rankPercentiles
        rankPercentiles = []
        for i in range(len(expectedValues[0])):
                downIndex = numpy.where(expectedValues[:,i]==0)
                upIndex = numpy.where(expectedValues[:,i]==2)
                unknownIndex = numpy.where(expectedValues[:,i]==3)
                for j in range(len(downIndex[0])):
                        sample = standardizedTFAmatrix[:,i]*-1
                        sampleRank = scipy.stats.percentileofscore(sample, sample[downIndex[0][j]])
                        rankPercentiles.append(sampleRank)
                for j in range(len(upIndex[0])):
                        sample = standardizedTFAmatrix[:,i]
                        sampleRank = scipy.stats.percentileofscore(sample, sample[upIndex[0][j]])
                        rankPercentiles.append(sampleRank)
                for j in range(len(unknownIndex[0])):
                        sample = [abs(x) for x in standardizedTFAmatrix[:,i]]
                        sampleRank = scipy.stats.percentileofscore(sample, sample[unknownIndex[0][j]])
                        rankPercentiles.append(sampleRank)
                                
        return rankPercentiles

"""
checks for positive correlations b/w inferred TFA values and measured mRNA values of the TFs

will also do 1000 bootstrap samplings of the samples/columns to control for bias from outlier samples

returns the number of positive correlations
        the median number of positive correlations among the bootstraps
"""
def checkCorrelation(tfaValues, mrnaValues):
        numPositive = 0
        for i in range(len(tfaValues)):
                corr, pVal = scipy.stats.pearsonr(tfaValues[i], mrnaValues[i])
                if corr > 0:
                        numPositive += 1
        
        bootstrappedResults = []
        for i in range(1000):
                bootstrapSampling = random.choices(range(len(tfaValues[0])), k=len(tfaValues[0]))
                bootstrappedTFA = tfaValues[:,bootstrapSampling]
                bootstrappedRNA = mrnaValues[:,bootstrapSampling]
                positiveCorrelation = 0
                for i in range(len(tfaValues)):
                        corr, pVal = scipy.stats.pearsonr(bootstrappedTFA[i], bootstrappedRNA[i])
                        if corr > 0:
                                positiveCorrelation += 1
                bootstrappedResults.append(positiveCorrelation)

        return [numPositive, numpy.median(bootstrappedResults)]

def main():
        [args, fileInputs] = parseInputs()
        #fileInputs = tfaValues, expectedValues, tfList, sampleList, mrnaValues

        #check data given, merge into route 1 if needed
        tfaValues = fileInputs[0]
        if numpy.shape(fileInputs[1]) == ():
                expectedValues = makeExpectationMatrix(fileInputs[2], fileInputs[3], args.perturbFlag)
        else:
                expectedValues = fileInputs[1]

        rankPercentiles = rankTFs(tfaValues, expectedValues)
        rankSignificance = scipy.stats.binom_test(len(rankPercentiles)/2, len(rankPercentiles),
                                                  1-numpy.median(rankPercentiles)/100.0, 'greater')
        print("Median rank percentile of perturbed TFs:\n\t", numpy.median(rankPercentiles),
                      "\tp-value: ", "{:.2e}".format(rankSignificance))
        
        [numCompared, numCorrect] = countDirected(tfaValues, expectedValues, args.wtIndex)
        if numCompared > 0:
                directedSignificance = scipy.stats.binom_test(numCorrect, numCompared, 0.5, 'greater')
                print("Percent of TFA changes in expected direction:\n\t", float(numCorrect)/numCompared,
                      "\tp-value: ", "{:.2e}".format(directedSignificance))
        else:
                print("No directed expectations to evaluate.")
                
        if numpy.shape(fileInputs[4]) != ():
                mrnaValues = fileInputs[4]
                [numPositive, medianBootstrapped] = checkCorrelation(tfaValues, mrnaValues)
                correlationSignificance = scipy.stats.binom_test(numPositive, len(tfaValues), 0.5, 'greater')
                print("Percent positive correlation between TFA and mRNA:\n\t",
                      float(numPositive)/float(len(tfaValues)),
                      "\tp-value: ", "{:.2e}".format(correlationSignificance))
                print("Percent positive correlation between TFA and mRNA (bootstrapped median):\n\t",
                      float(medianBootstrapped)/float(len(tfaValues)))





main()



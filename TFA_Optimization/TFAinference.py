#!/usr/bin/python

#learn TFA matrix values from 
#	symbolic TFA matrix, 
#	CS matrix, 
#	and measured log expression values


from gurobipy import *
import numpy
import random
import argparse
import time


def parseInputs():
  parser = argparse.ArgumentParser(description='Infer CS and TFA values')
  parser.add_argument('--expression', '-e', type=str, 
		      dest='expressionFile', action='store', required=True,
                      help='path to csv file of expression data')
  parser.add_argument('--crossSet', '-v', type=str, 
		      dest='crossFile', action='store', 
                      help='path to csv file of expression data for cross checking overfitting')
  parser.add_argument('--WT', '-w', '--wt', type=int, 
		      dest='wtIndex', action='store', default=-1,
                      help='column index for the WT or reference sample')
  parser.add_argument('--binaryCS', '-c', type=str, 
		      dest='csFile', action='store',
                      help='csv file of edges between TFs (columns) and target genes (rows)\
                            \n -1: repressive edge\t 0: no edge\t 1: activating edge')
  parser.add_argument('--startingCS', '-s', type=str, 
		      dest='startFile', action='store',
                      help='csv file of edge strengths between TFs (columns) and target genes (rows)\
                            \n optimization will start with these values')
  parser.add_argument('--binaryTFA', '-a', type=str, 
                      dest='tfaFile', action='store',
                      help='csv file of known activity behavior for TFs (rows) in samples (columns)\
                            \n 0: no activity\t 1: no prior expectation\t 2: greater activity than WT/reference')
  #parser.add_argument('--freezeCS', '-f', type=int, 
  #			dest='freezeFlag', action='store',
  #              	default=0, choices = [0,1],
  #              	help='meant for second round fitting where CS values are held constant\
  #                            \n 0: optimize all \t1: only optimize TFA and baselines')
  parser.add_argument('--constrainMOR', '-m', type=int, 
			dest='morFlag', action='store',
                	default=0, choices = [0,1],
                	help='constrain mode of regulation by using signs from binaryCS\
                	\n 0: use signs to constrain CS values\t 1: don\'t constrain CS values')
  parser.add_argument('--maxIterations', '-i', type=int, 
		      dest='maxIter', action='store', default=100,
                      help='maximum number of alternative iterations of optimization')
  parser.add_argument('--fileTag', '-t', type=str, 
		      dest='fileTag', action='store', default='',
                      help='string to help label log and output files')

  args = parser.parse_args()
  #check for required expression data
  try:
    expressionValues = numpy.genfromtxt(args.expressionFile, delimiter=',')
  except IOError as e:
    print("Error:\t", e)
    quit()
  #check for binaryCS and/or startingCS file
  if args.startFile != None:
    try:
      startCS = numpy.genfromtxt(args.startFile, delimiter=',', dtype=float)
    except IOError as e:
      print("Error:\t", e)
      quit()
  elif args.csFile != None:
    try:
      binaryCS = numpy.genfromtxt(args.csFile, delimiter=',', dtype=float)
      startCS = makeRandomStart(binaryCS, args.morFlag==0)
    except IOError as e:
      print("Error:\t", e)
      quit()
  else:
    print("Either a binaryCS file or a startingCS file is required")
    quit()
  #check for binaryTFA data
  if args.tfaFile != None:
    try:
      binaryTFA = numpy.genfromtxt(args.tfaFile, delimiter=',')
    except IOError as e:
      print("Error:\t", e)
      quit()
  else:
    binaryTFA = numpy.ones((len(startCS[0])-1,len(expressionValues[0])))
    
  if len(expressionValues) != len(startCS):
    print("Error: number of genes (",len(expressionValues),"rows) in expression data does not match (",len(startCS),"rows) CS data")
    quit()
  if len(expressionValues[0]) != len(binaryTFA[0]):
    print("Error: number of samples (",len(expressionValues[0]),"columns) in expression data does not match (",len(binaryTFA[0]),"columns) TFA data")
    quit()
  if len(startCS[0])-1 != len(binaryTFA):
    print("Error: number of TFs ("+str(len(startCS[0])-1)+"columns) in CS data does not match ("+str(len(binaryTFA))+"rows) TFA data")
    quit()

  if abs(args.wtIndex) > len(expressionValues[0]):
    print("Error: index of reference sample outside the bounds of given expression data")
    quit()

  if args.morFlag == 0:
    print("will constrain the signs/mode of regulation for CS values")

  return [args, expressionValues, startCS, binaryTFA]


  

"""
Input:
  a numpy array of binaryCS values
  a boolean for whether to keep cs sign constraint
Output:
  a numpy array to use as a random start CS

makes a random start point for cs matrix and saves as a csv file
"""
def makeRandomStart(binaryCS, signFlag):
  startingCS = numpy.hstack((numpy.copy(binaryCS), numpy.ones((len(binaryCS),1))))
  for i in range(len(startingCS)):
    for j in range(len(binaryCS[i])):
      randVal = random.uniform(-10,10)
      if signFlag and abs(binaryCS[i][j])==1:
        startingCS[i,j]*=abs(randVal)
      else:
        startingCS[i,j]*=randVal
    startingCS[i,-1] = random.uniform(-10,10)
  return startingCS
  

"""
Input:
  two data matrices
  the first is assumed to be true data
  the second is assumed to be predicted/learned/fitted data
Output:
  list of var explained and SSE

calculates total error and variance explained between predicted and true data
"""
def calcError(data, dataLearned):
  numerator = 0
  denominator = 0
  for i in range(len(dataLearned)):
    sst = numpy.var(data[i])*len(data[i])
    sse = sum((data[i]-dataLearned[i])**2)
    numerator += sse
    denominator += sst
  return [1-(numerator/denominator), numerator]








"""
Input:
  binary activity matrix A
  latest control strength matrix C, including baseline
  expression matrix data
  freezeFlag to indicate whether doing second step fitting (or else constrain avg TFA)
  fileTag for use w/logging gurobi output
Output:
  False if learning failed
  learned A matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""
def learnTFA(A, lastC, data, freezeFlag, fileTag):
  numGenes = len(data)
  numSamples = len(data[0])
  numTFs = len(A)

  # initialize gurobi model
  model = Model()
  model.setParam('LogToConsole', False)
  model.setParam('LogFile', 'logFiles/'+fileTag+".log")

  #threads for parallelization
  model.setParam('Threads', 1)

  # Add tfa variables to the model
  varsMatrix = []   # holds the activity matrix, with pointers to coeff where relevant
  for i in range(numTFs):
    constraintCounter = 0   # counts the number of coeff in a row
    varsMatrix.append([])   # start a new row in the activity matrix
    constraint = LinExpr()    # initialize the constraint that each row's avg coeff value is 1
    for j in range(numSamples):
      if A[i][j]==0:
        varsMatrix[i].append(0)
      else:
        v = model.addVar(lb=0.0001, vtype=GRB.CONTINUOUS, name='A['+str(i)+','+str(j)+']')
        varsMatrix[i].append(v)
        constraint += v
        constraintCounter += 1
    if not freezeFlag:  # add the scaling constraint if not second round fitting
      model.addConstr(constraint/constraintCounter, GRB.EQUAL, 1.0, "c"+str(i))
    model.update()

  # Add overexpression perturb constraints to the model
  for i in range(numTFs):
    for j in range(numSamples):
      if A[i][j]==2:
        model.addConstr(varsMatrix[i][j], GRB.GREATER_EQUAL, 1.1*varsMatrix[i][-1])
  model.update()

  # Populate objective
  obj = QuadExpr()
  for i in range(numGenes):
    for j in range(numSamples):
      geneExpr = LinExpr()
      geneExpr += lastC[i][numTFs]
      for k in range(numTFs):
        if type(varsMatrix[k][j])==Var:
          geneExpr += lastC[i][k]*varsMatrix[k][j]
      geneError = data[i][j] - geneExpr
      obj += geneError * geneError
  model.setObjective(obj)
  model.update()

  # Solve
  try:
  	model.optimize()
  except:
  	return False

  # Write model to a file
  model.write('modelFiles/'+fileTag+'.lp')

  # check that optimization succeeded
  if model.status != GRB.Status.OPTIMAL:
    return False

  #convert back to matrix
  Atemp = []
  for i in range(numTFs):
    Atemp.append([])
    for j in range(numSamples):
      if A[i][j] == 0:
        Atemp[i].append(0)
      else:
        Atemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])

  return Atemp



"""
Input:
  latest activity matrix A
  binary cs matrix C, or latest cs matrix if validating
  expression matrix data
  morFlag to indicate whether to constrain mode of regulation (sign of CS values)
  freezeFlag to indicate whether doing second step fitting (or else constrain avg TFA)
  fileTag for use w/logging gurobi output
Output:
  False if learning failed
  learned CS matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""
def learnCS(lastA, C, data, morFlag, fileTag):
  numGenes = len(data)
  numSamples = len(data[0])
  numTFs = len(lastA)

  # Initialize the model
  model = Model()
  model.setParam('LogToConsole', False)
  model.setParam('LogFile', "logFiles/"+fileTag+".log")
 
  model.setParam('Threads', 1)
 
  # Add cs variables to the model
  varsMatrix = []   # holds the cs matrix, with pointers to coeff where relevant
  lassoConstraint = LinExpr()   # Intialize the LASSO constraint
  for i in range(numGenes):
    varsMatrix.append([])   # add a row to the cs matrix
    for j in range(numTFs):
        #if freezeFlag==1:
        # varsMatrix[i].append(C[i][j])
        #else:
        if C[i][j]==0:  #no influence
          varsMatrix[i].append(0)
        elif not morFlag: #an influence to be learned, no sign to constrain
          v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='C['+str(i)+','+str(j)+']')
          varsMatrix[i].append(v)
        else:
          if C[i][j] > 0: #an influence to be learned, activating
            v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='C['+str(i)+','+str(j)+']')
            varsMatrix[i].append(v)
          else: #an influence to be learned, repressing
            v = model.addVar(lb=-GRB.INFINITY, ub=-0.0001, vtype=GRB.CONTINUOUS, name='C['+str(i)+','+str(j)+']')
            varsMatrix[i].append(v)
    #learning baseline expression
    v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='C['+str(i)+','+str(numTFs)+']')
    varsMatrix[i].append(v)

  model.update()

  # Populate objective
  obj = QuadExpr()
  for i in range(numGenes):
    for j in range(numSamples):
      geneExpr = LinExpr()
      geneExpr += varsMatrix[i][numTFs]
      for k in range(numTFs):
          geneExpr += varsMatrix[i][k]*lastA[k][j]
      geneError = data[i][j] - geneExpr
      obj += geneError * geneError

  model.setObjective(obj)
  model.update()

  # Solve
  try:
    model.optimize()
  except:
    return False

  # Write model to a file
  model.write('modelFiles/'+fileTag+'.lp')

  # check that optimization succeeded
  if model.status != GRB.Status.OPTIMAL:
    return False

  #convert back to matrix
  """
  if freezeFlag==1:
    Ctemp = []
    for i in range(numGenes):
      Ctemp.append([])
      for j in range(numTFs+1):
        if j<numTFs:
          Ctemp[i].append(C[i][j])
        else:
          Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
  else:
  """
  Ctemp = []
  for i in range(numGenes):
      Ctemp.append([])
      for j in range(numTFs+1):
        if j<numTFs and C[i][j] == 0:
          Ctemp[i].append(0)
        else:
          Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
    
  return Ctemp


def main():
  start = time.time()
  [args, expressionValues, startCS, binaryTFA] = parseInputs()
  crossFlag = False
  if args.crossFile != None:
    try:
      crossExpression = numpy.genfromtxt(args.crossFile, delimiter=',', dtype=float)
      crossBinaryTFA = numpy.ones((len(startCS[0])-1,len(crossExpression[0])))
      crossFlag = True
    except IOError as e:
      print("Error:\t", e)
      quit()
  numGenes = len(startCS)
  numSamples = len(expressionValues[0])
  numTFs = len(binaryTFA)

  Atemp = learnTFA(binaryTFA, startCS, expressionValues, False, args.fileTag+"TFA")
  if Atemp == False:
    print("Could not learn the activity matrix")
    return

  numpy.savetxt("logFiles/"+args.fileTag+"CSiteration0.csv", startCS, delimiter=',', fmt='%1.5f')
  numpy.savetxt("logFiles/"+args.fileTag+"TFAiteration0.csv", Atemp, delimiter=',', fmt='%1.5f')

  predictedExpression = numpy.matmul(startCS, numpy.append(Atemp, numpy.ones((1,len(Atemp[0]))),axis=0))
  [ve, sse] = calcError(expressionValues, predictedExpression)

  print("Variance explained and error from random CS values: ", round(ve, 3), round(sse,3))

  if crossFlag:
      crossA = learnTFA(crossBinaryTFA, startCS, crossExpression, True, args.fileTag+"CrossCheck")
      if crossA == False:
        print("Could not cross check")
      else:
        predictedExpression = numpy.matmul(startCS, numpy.append(crossA, numpy.ones((1,len(crossA[0]))),axis=0))
        [ve, sse] = calcError(crossExpression, predictedExpression)
        print("cross check: ", round(numpy.mean(crossA[0]),2), round(ve, 3), round(sse,3))
  
  for itr in range(1,args.maxIter):
    print("\niteration ", itr, "\n")

    Ctemp = learnCS(Atemp, startCS, expressionValues, args.morFlag==0, args.fileTag+"CS")
    if Ctemp == False:
      print("Could not learn the control strength matrix")
      return
    
    Atemp = learnTFA(binaryTFA, Ctemp, expressionValues, False, args.fileTag+"TFA")
    if Atemp == False:
      print("Could not learn the activity matrix")
      return

    numpy.savetxt("logFiles/"+args.fileTag+"CSiteration"+str(itr)+".csv", Ctemp, delimiter=',', fmt='%1.5f')
    numpy.savetxt("logFiles/"+args.fileTag+"TFAiteration"+str(itr)+".csv", Atemp, delimiter=',', fmt='%1.5f')

    predictedExpression = numpy.matmul(Ctemp, numpy.append(Atemp, numpy.ones((1,len(Atemp[0]))),axis=0))
    [ve, sse] = calcError(expressionValues, predictedExpression)

    print("Variance explained and error at iteration ", itr, ": ", round(ve, 3), round(sse,3))

    if crossFlag:
      crossA = learnTFA(crossBinaryTFA, Ctemp, crossExpression, True, args.fileTag+"CrossCheck")
      if crossA == False:
        print("Could not cross check")
      else:
        predictedExpression = numpy.matmul(Ctemp, numpy.append(crossA, numpy.ones((1,len(crossA[0]))),axis=0))
        [ve, sse] = calcError(crossExpression, predictedExpression)
        print("cross check: ", round(numpy.mean(crossA[0]),2), round(ve, 3), round(sse,3))

  end = time.time()
  print(round(end-start))
main()











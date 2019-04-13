
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math
import pyBigWig
import numpy as np
from pysam import Samfile

# Input
ext = int(sys.argv[1])
featureSummitFileName = sys.argv[2]
bamNames = sys.argv[3].split(",")
bamCounts = [float(e) for e in sys.argv[4].split(",")]
bamList = sys.argv[5].split(",")
outputFileName = sys.argv[6]

###################################################################################################
# Functions
###################################################################################################

def fetchSignal(bamFile, region):
  return float(sum(1 for _ in bamFile.fetch(region[0],region[1],region[2])))

def fetchSignalBw(bwFile, region):
  valuesVec = bwFile.values(region[0], region[1], region[2])
  for i in range(0, len(valuesVec)):
    if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): valuesVec[i] = 0.0
  return sum(valuesVec)

###################################################################################################
# Intersection table
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Fetching regions
featureSummitFile = open(featureSummitFileName,"r")
regionList = []
for line in featureSummitFile:
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  region = [ll[0], int(ll[1])-ext, int(ll[2])+ext]
  if(int(region[1]) < 0): continue
  regionList.append(region)
featureSummitFile.close()

# Creating table
matrix = []
for i in range(0,len(bamList)):
  inputBamFileName = bamList[i]
  correctFactor = bamCounts[i]/1000000
  if(".bam" in inputBamFileName):
    bamFile = Samfile(inputBamFileName,"rb")
    vec = []
    for region in regionList:
      try: bamSignal = fetchSignal(bamFile, region) / correctFactor
      except Exception: bamSignal = 0
      vec.append(bamSignal)
  else:
    bamFile = pyBigWig.open(inputBamFileName)
    vec = []
    for region in regionList:
      #try: 
      bamSignal = fetchSignalBw(bamFile, region) / correctFactor
      #except Exception: bamSignal = 0
      vec.append(bamSignal) 
  matrix.append(vec)
  bamFile.close()
outputFile = open(outputFileName,"w")
outputFile.write("\t".join(bamNames)+"\n")
for j in range(0,len(matrix[0])):
  vec = []
  for i in range(0,len(matrix)):
    try: vec.append(str(matrix[i][j]))
    except Exception: vec.append("NA")
  outputFile.write("\t".join(vec)+"\n")
outputFile.close()



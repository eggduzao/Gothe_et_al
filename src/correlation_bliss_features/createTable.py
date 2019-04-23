
###################################################################################################
# Input
###################################################################################################

# Import
from __future__ import print_function
import os
import sys
import math
import pyBigWig
import numpy as np
from pysam import Samfile

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

def create_table(half_ext, feature_summit_file_name, bam_names, bam_counts, bam_list, output_file_name):

  # Initialization
  outLoc = "/".join(output_file_name.split("\t")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Allowed chromosomes
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Fetching regions
  featureSummitFile = open(feature_summit_file_name,"r")
  regionList = []
  for line in featureSummitFile:
    ll = line.strip().split("\t")
    if(ll[0] not in chrList): continue
    region = [ll[0], int(ll[1])-half_ext, int(ll[2])+half_ext]
    if(int(region[1]) < 0): continue
    regionList.append(region)
  featureSummitFile.close()

  # Creating table
  matrix = []
  for i in range(0,len(bam_list)):
    inputBamFileName = bam_list[i]
    correctFactor = int(bam_counts[i])/1000000
    extension = inputBamFileName.split(".")[-1]
    if(extension == "bam"):
      bamFile = Samfile(inputBamFileName,"rb")
      vec = []
      for region in regionList:
        try: bamSignal = fetchSignal(bamFile, region) / correctFactor
        except Exception: bamSignal = 0
        vec.append(bamSignal)
    elif(extension == "bw" or extension == "bigwig"):
      bamFile = pyBigWig.open(inputBamFileName)
      vec = []
      for region in regionList:
        try: bamSignal = fetchSignalBw(bamFile, region) / correctFactor
        except Exception: bamSignal = 0
      vec.append(bamSignal)
    else: print("The tool supports only BAM or BIGWIG files.")
    matrix.append(vec)
    bamFile.close()
  outputFile = open(output_file_name,"w")
  outputFile.write("\t".join(bam_names)+"\n")
  for j in range(0,len(matrix[0])):
    vec = []
    for i in range(0,len(matrix)):
      try: vec.append(str(matrix[i][j]))
      except Exception: vec.append("NA")
    outputFile.write("\t".join(vec)+"\n")
  outputFile.close()

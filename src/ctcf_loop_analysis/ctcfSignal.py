
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

###################################################################################################
# Functions
###################################################################################################

def fetchTotalSignalBam(bamFile, region):
  totalSignal = 0.0
  for read in bamFile.fetch(region[0], region[1], region[2]): totalSignal += 1.0
  return totalSignal

def fetchSignalBam(bamFile, region, ext):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  for read in bamFile.fetch(region[0], region[1], region[2]):
    vecReadLoc = read.reference_start - region[1]
    for i in range(max(vecReadLoc - ext, 0), min(vecReadLoc + ext, len(returnVec))): returnVec[i] += 1.0
  return returnVec

def ctcf_signal(): pass

###################################################################################################
# Reading genes
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Alias dictionary
aliasDict = dict()
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g] = value
aliasFile.close()

# Gene dictionary
geneDict = dict()
geneFile = open(geneFileName,"rU")
for line in geneFile:
  ll = line.strip().split("\t")
  try: geneDict[aliasDict[ll[3]]] = ll
  except Exception: continue
geneFile.close()

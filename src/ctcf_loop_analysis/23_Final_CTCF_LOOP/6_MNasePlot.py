
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
higherSide = sys.argv[1]
ctcfExt = int(sys.argv[2])
ctcfFileName = sys.argv[3] # motifs
etoBamFileName = sys.argv[4]
ctcfBamFileName = sys.argv[5] # ChIP-seq
mnaseBwFileName = sys.argv[6]
outputEtoFileName = sys.argv[7]
outputCtcfFileName = sys.argv[8]
outputMNaseFileName = sys.argv[9]

# Initialization
outLoc = "/".join(outputEtoFileName.split("/")[:-1])+"/"
command = "mkdir -p "+outLoc
os.system(command)
bamExt = 5

###################################################################################################
# Functions
###################################################################################################

def fetchSignalBam(bamFile, region, ext):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  for read in bamFile.fetch(region[0], region[1], region[2]):
    vecReadLoc = read.reference_start - region[1]
    for i in range(max(vecReadLoc - ext, 0), min(vecReadLoc + ext, len(returnVec))): returnVec[i] += 1.0
  return returnVec

def fetchSignalBw(bwFile, region):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  valuesVec = bwFile.values(region[0], region[1], region[2])
  for i in range(0, len(valuesVec)):
    if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): returnVec[i] = 0.0
    else: returnVec[i] = valuesVec[i]
  return returnVec

###################################################################################################
# Creating table
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Initialization
etoBamFile = Samfile(etoBamFileName, "rb")
ctcfBamFile = Samfile(ctcfBamFileName, "rb")
mnaseBwFile = pyBigWig.open(mnaseBwFileName)

# Fetching the bam signal in all categories
# CTCF_CHR, CTCF_P1, CTCF_P2, CTCF_STR, [SIGNAL...]
ctcfFile = open(ctcfFileName, "rU")
outputEtoFile = open(outputEtoFileName,"w")
outputCtcfFile = open(outputCtcfFileName,"w")
outputMNaseFile = open(outputMNaseFileName,"w")
for line in ctcfFile:

  # Initialization
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  ctcfChrom = ll[0]; ctcfP1 = int(ll[1]); ctcfP2 = int(ll[2]); ctcfStrand = ll[5]
  mid = (ctcfP1+ctcfP2)/2
  vector = [ctcfChrom, ctcfP1, ctcfP2, ctcfStrand]

  # Fetching DSB signal
  region11 = [ctcfChrom, mid-ctcfExt, mid-200]; region12 = [ctcfChrom, mid-200, mid-100]; region13 = [ctcfChrom, mid-100, mid]
  region21 = [ctcfChrom, mid, mid+100]; region22 = [ctcfChrom, mid+100, mid+200]; region23 = [ctcfChrom, mid+200, mid+ctcfExt]; 
  if(higherSide == "left"):
    r11p = 0.6
    r12p = 0.8
    r13p = 1.0
    r21p = 0.5
    r22p = 0.6
    r23p = 0.6
  elif(higherSide == "right"):
    r11p = 0.6
    r12p = 0.6
    r13p = 0.5
    r21p = 1.0
    r22p = 0.8
    r23p = 0.6
  elif(higherSide == "rest"):
    r11p = 0.3
    r12p = 0.3
    r13p = 0.2
    r21p = 0.2
    r22p = 0.3
    r23p = 0.3
  else:
    r11p = 0.6
    r12p = 0.6
    r13p = 0.6
    r21p = 0.6
    r22p = 0.6
    r23p = 0.6
  signal1 = [e*r11p for e in fetchSignalBam(etoBamFile, region11, bamExt)] + [e*r12p for e in fetchSignalBam(etoBamFile, region12, bamExt)] + [e*r13p for e in fetchSignalBam(etoBamFile, region13, bamExt)]
  signal2 = [e*r21p for e in fetchSignalBam(etoBamFile, region21, bamExt)] + [e*r22p for e in fetchSignalBam(etoBamFile, region22, bamExt)] + [e*r23p for e in fetchSignalBam(etoBamFile, region23, bamExt)]
  dsbSignal = signal1 + signal2

  # CTCF / MNase signal
  region = [ctcfChrom, mid-ctcfExt, mid+ctcfExt]
  ctcfSignal = fetchSignalBam(ctcfBamFile, region, 50)
  mnaseSignal = fetchSignalBw(mnaseBwFile, region)

  # Writing vector
  outputEtoFile.write("\t".join([str(e) for e in vector + dsbSignal])+"\n")
  outputCtcfFile.write("\t".join([str(e) for e in vector + ctcfSignal])+"\n")
  outputMNaseFile.write("\t".join([str(e) for e in vector + mnaseSignal])+"\n")

# Closing all files
ctcfFile.close()
outputEtoFile.close()
outputCtcfFile.close()
outputMNaseFile.close()
etoBamFile.close()
ctcfBamFile.close()
mnaseBwFile.close()



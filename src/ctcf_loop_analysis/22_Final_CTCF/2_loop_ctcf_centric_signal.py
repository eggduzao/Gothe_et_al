
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
ctcfExt = int(sys.argv[1])
loopFileName = sys.argv[2]
signalFileName = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
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

###################################################################################################
# Execution
###################################################################################################

# Initialization
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]
signalFile = Samfile(signalFileName, "rb")
outputFile = open(outputFileName, "w")

# Fetching CTCF sites
# CTCF_CHR, CTCF_P1, CTCF_P2, LOOP_SCORE, [SIGNAL...]
loopFile = open(loopFileName, "rU")
loopFile.readline()
for line in loopFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = "chr"+ll[0]
  if(chrom not in chrList): continue
  ctcfx1 = ll[20]; ctcfx2 = ll[21]; ctcfxo = ll[23]
  ctcfy1 = ll[25]; ctcfy2 = ll[26]; ctcfyo = ll[28]
  loopScore = ll[7]

  # First motif
  if(ctcfxo != "NA"):
    mid = (int(ctcfx1)+int(ctcfx2)) / 2
    region = [chrom, mid-ctcfExt, mid+ctcfExt]
    signal = fetchSignalBam(signalFile, region, bamExt)
    vector = [chrom, ctcfx1, ctcfx2, loopScore] + signal
    outputFile.write("\t".join([str(e) for e in vector])+"\n")

  # Second motif
  if(ctcfyo != "NA"):
    mid = (int(ctcfy1)+int(ctcfy2)) / 2
    region = [chrom, mid-ctcfExt, mid+ctcfExt]
    signal = fetchSignalBam(signalFile, region, bamExt)[::-1]
    vector = [chrom, ctcfy1, ctcfy2, loopScore] + signal
    outputFile.write("\t".join([str(e) for e in vector])+"\n")

# Termination
loopFile.close()
signalFile.close()
outputFile.close()



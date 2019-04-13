
# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# INPUT
###################################################################################################

# Input 
ctcfBamFileName = sys.argv[1]
ctcfBedFileName = sys.argv[2]
featureFileName = sys.argv[3]
tempLoc = sys.argv[4]
outputPrefix = sys.argv[5]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def getBestMotifList(bamFile, region):
  motifList = []
  bestScore = -9999
  motifFetch = bamFile.fetch(region[0], region[1], region[2])
  for read in motifFetch:
    score = float(read.qname.split(":")[-1])
    if(score < bestScore): continue
    elif(score == bestScore):
      strand = "+"
      if(read.is_reverse): strand = "-"
      motifList.append([region[0], str(read.pos), str(read.aend), str(score), strand])
      bestScore = score
    else:
      strand = "+"
      if(read.is_reverse): strand = "-"
      motifList = [[region[0], str(read.pos), str(read.aend), str(score), strand]]
      bestScore = score
  return motifList

###################################################################################################
# EXECUTION
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Input files
ctcfBamFile = Samfile(ctcfBamFileName, "rb")

# Output files
outFileMFRTFRName = outputPrefix+"_MFR_TFR.bed"
outFileMFTFRName = outputPrefix+"_MF_TFR.bed"
outFileMRTFRName = outputPrefix+"_MR_TFR.bed"
outFileMFRTFName = outputPrefix+"_MFR_TF.bed"
outFileMFRTRName = outputPrefix+"_MFR_TR.bed"
outFileMFTFMRTRName = outputPrefix+"_MFTF_MRTR.bed"
outFileMFTRMRTFName = outputPrefix+"_MFTR_MRTF.bed"
outFileMFRTFR = open(outFileMFRTFRName, "w")
outFileMFTFR = open(outFileMFTFRName, "w")
outFileMRTFR = open(outFileMRTFRName, "w")
outFileMFRTF = open(outFileMFRTFName, "w")
outFileMFRTR = open(outFileMFRTRName, "w")
outFileMFTFMRTR = open(outFileMFTFMRTRName, "w")
outFileMFTRMRTF = open(outFileMFTRMRTFName, "w")

# Iterating on features
featureFile = open(featureFileName, "rU")
for line in featureFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; name = ll[3]; score = ll[4]; strand = ll[5]
  if(chrom not in chrList): continue

  # Best motif list
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])

  # Output files writing
  for motif in bestMotifList:
    mStrand = motif[4]
    if(strand == "+" and mStrand == "+"):
      outFileMFRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFRTF.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTFMRTR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
    elif(strand == "+" and mStrand == "-"):
      outFileMFRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFRTF.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTRMRTF.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
    elif(strand == "-" and mStrand == "+"):
      outFileMFRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFRTR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTRMRTF.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
    elif(strand == "-" and mStrand == "-"):
      outFileMFRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMRTFR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFRTR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")
      outFileMFTFMRTR.write("\t".join([motif[0], motif[1], motif[2], name, motif[3], motif[4]])+"\n")

# Closing files
featureFile.close()
outFileMFRTFR.close()
outFileMFTFR.close()
outFileMRTFR.close()
outFileMFRTF.close()
outFileMFRTR.close()
outFileMFTFMRTR.close()
outFileMFTRMRTF.close()
ctcfBamFile.close()

# Rest output file
outFileRestMFName = outputPrefix+"_restMF.bed"
outFileRestMRName = outputPrefix+"_restMR.bed"

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+outFileMFRTFRName+" > "+tempSortFileName
os.system(command)

# Intersection
tempIntFileName = tempLoc+"tempIntFileName.bed"
command = "intersectBed -wa -v -a "+ctcfBedFileName+" -b "+tempSortFileName+" > "+tempIntFileName
os.system(command)

# Fetching rest files
command = "awk '$6 == \"+\"' "+tempIntFileName+" > "+outFileRestMFName
os.system(command)
command = "awk '$6 == \"-\"' "+tempIntFileName+" > "+outFileRestMRName
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)




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
loopFileName = sys.argv[3]
tempLoc = sys.argv[4]
outputPrefix = sys.argv[5]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

###################################################################################################
# FUNCTIONS
###################################################################################################

def getMotifList(bamFile, region):
  motifList = []
  motifFetch = bamFile.fetch(region[0], region[1], region[2])
  for read in motifFetch:
    strand = "+"
    if(read.is_reverse): strand = "-"
    motifList.append([region[0], str(read.pos), str(read.aend), read.qname.split(":")[0], read.qname.split(":")[1], strand])
  return motifList

###################################################################################################
# EXECUTION
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Input files
ctcfBamFile = Samfile(ctcfBamFileName, "rb")

# Output files
outFileTIOName = outputPrefix+"_TIO.bed"
outFileTOName = outputPrefix+"_TO.bed"
outFileTIName = outputPrefix+"_TI.bed"
outFileTIO = open(outFileTIOName, "w")
outFileTO = open(outFileTOName, "w")
outFileTI = open(outFileTIName, "w")

# Iterating on features
loopFile = open(loopFileName, "rU")
for line in loopFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom1 = ll[0]; p11 = ll[1]; p12 = ll[2]
  chrom2 = ll[3]; p21 = ll[4]; p22 = ll[5]
  if(chrom1 not in chrList): continue

  # Motif list
  motifList1 = getMotifList(ctcfBamFile, [chrom1, int(p11), int(p12)])
  motifList2 = getMotifList(ctcfBamFile, [chrom2, int(p21), int(p22)])

  # Output files writing
  for motif in motifList1:
    outFileTIO.write("\t".join(motif)+"\n")
    outFileTO.write("\t".join(motif)+"\n")
  for motif in motifList2:
    outFileTIO.write("\t".join(motif)+"\n")
    outFileTI.write("\t".join(motif)+"\n")

# Closing files
loopFile.close()
outFileTIO.close()
outFileTO.close()
outFileTI.close()
ctcfBamFile.close()

# Rest output file
outFileRestName = outputPrefix+"_rest.bed"

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+outFileTIOName+" > "+tempSortFileName
os.system(command)

# Intersection
command = "intersectBed -wa -v -a "+ctcfBedFileName+" -b "+tempSortFileName+" > "+outFileRestName
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)



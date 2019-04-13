
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
loopBedFileName = sys.argv[3]
activeGeneBedFileName = sys.argv[4]
inactiveGeneBedFileName = sys.argv[5]
activePromBedFileName = sys.argv[6]
inactivePromBedFileName = sys.argv[7]
intergenicBedFileName = sys.argv[8]
tempLoc = sys.argv[9]
outputPrefix = sys.argv[10]

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
    rr = read.qname.split(":")
    name = ":".join(rr[:-1])
    score = float(rr[-1])
    if(score < bestScore): continue
    elif(score == bestScore):
      strand = "+"
      if(read.is_reverse): strand = "-"
      motifList.append([region[0], str(read.pos), str(read.aend), name, str(score), strand])
      bestScore = score
    else:
      strand = "+"
      if(read.is_reverse): strand = "-"
      motifList = [[region[0], str(read.pos), str(read.aend), name, str(score), strand]]
      bestScore = score
  return motifList

###################################################################################################
# EXECUTION
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Input bam files
ctcfBamFile = Samfile(ctcfBamFileName, "rb")

###################################################################################################
# LOOPS
###################################################################################################

# Output files
allCtcfInLoopFileName = outputPrefix+"allCtcfInLoop.bed"
allCtcfOutLoopFileName = outputPrefix+"allCtcfOutLoop.bed"
allCtcfInLoopFile = open(allCtcfInLoopFileName, "w")
allCtcfOutLoopFile = open(allCtcfOutLoopFileName, "w")

# Iterating on features
featureFile = open(loopBedFileName, "rU")
for line in featureFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]
  if(chrom not in chrList): continue
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])
  for motif in bestMotifList: allCtcfInLoopFile.write("\t".join(motif)+"\n")

# Closing files
featureFile.close()
allCtcfInLoopFile.close()
allCtcfOutLoopFile.close()

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+allCtcfInLoopFileName+" > "+tempSortFileName
os.system(command)

# Intersection
command = "intersectBed -wa -v -a "+ctcfBedFileName+" -b "+tempSortFileName+" > "+allCtcfOutLoopFileName
os.system(command)

###################################################################################################
# GENE
###################################################################################################

# Output files
allCtcfAllGeneFileName = outputPrefix+"allCtcfAllGene.bed"
allCtcfPlusGeneFileName = outputPrefix+"allCtcfPlusGene.bed"
allCtcfMinusGeneFileName = outputPrefix+"allCtcfMinusGene.bed"
allCtcfAllGeneFile = open(allCtcfAllGeneFileName, "w")
allCtcfPlusGeneFile = open(allCtcfPlusGeneFileName, "w")
allCtcfMinusGeneFile = open(allCtcfMinusGeneFileName, "w")

# Iterating on features
featureFile = open(activeGeneBedFileName, "rU")
for line in featureFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; name = ll[3]; score = ll[4]; strand = ll[5]
  if(chrom not in chrList): continue
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])
  for motif in bestMotifList: allCtcfPlusGeneFile.write("\t".join([motif[0], motif[1], motif[2], name, motif[4], motif[5]])+"\n")
featureFile.close()
featureFile = open(inactiveGeneBedFileName, "rU")
for line in featureFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; name = ll[3]; score = ll[4]; strand = ll[5]
  if(chrom not in chrList): continue
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])
  for motif in bestMotifList: allCtcfMinusGeneFile.write("\t".join([motif[0], motif[1], motif[2], name, motif[4], motif[5]])+"\n")

# Closing files
featureFile.close()
allCtcfAllGeneFile.close()
allCtcfPlusGeneFile.close()
allCtcfMinusGeneFile.close()

# Cat
tempCatFileName = tempLoc+"tempCatFileName.bed"
command = "cat "+allCtcfPlusGeneFileName+" "+allCtcfMinusGeneFileName+" > "+tempCatFileName
os.system(command)

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+tempCatFileName+" > "+tempSortFileName
os.system(command)

# Merge
command = "mergeBed -c 4,5,6 -o first,mean,first -i "+tempSortFileName+" > "+allCtcfAllGeneFileName
os.system(command)

###################################################################################################
# PROMOTER
###################################################################################################

# Output files
allCtcfAllPromFileName = outputPrefix+"allCtcfAllProm.bed"
allCtcfPlusPromFileName = outputPrefix+"allCtcfPlusProm.bed"
allCtcfMinusPromFileName = outputPrefix+"allCtcfMinusProm.bed"
allCtcfAllPromFile = open(allCtcfAllPromFileName, "w")
allCtcfPlusPromFile = open(allCtcfPlusPromFileName, "w")
allCtcfMinusPromFile = open(allCtcfMinusPromFileName, "w")

# Iterating on features
featureFile = open(activePromBedFileName, "rU")
for line in featureFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; name = ll[3]; score = ll[4]; strand = ll[5]
  if(chrom not in chrList): continue
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])
  for motif in bestMotifList: allCtcfPlusPromFile.write("\t".join([motif[0], motif[1], motif[2], name, motif[4], motif[5]])+"\n")
featureFile.close()
featureFile = open(inactivePromBedFileName, "rU")
for line in featureFile:
  ll = line.strip().split("\t")
  chrom = ll[0]; p1 = ll[1]; p2 = ll[2]; name = ll[3]; score = ll[4]; strand = ll[5]
  if(chrom not in chrList): continue
  bestMotifList = getBestMotifList(ctcfBamFile, [chrom, int(p1), int(p2)])
  for motif in bestMotifList: allCtcfMinusPromFile.write("\t".join([motif[0], motif[1], motif[2], name, motif[4], motif[5]])+"\n")

# Closing files
featureFile.close()
allCtcfAllPromFile.close()
allCtcfPlusPromFile.close()
allCtcfMinusPromFile.close()
ctcfBamFile.close()

# Cat
tempCatFileName = tempLoc+"tempCatFileName.bed"
command = "cat "+allCtcfPlusPromFileName+" "+allCtcfMinusPromFileName+" > "+tempCatFileName
os.system(command)

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+tempCatFileName+" > "+tempSortFileName
os.system(command)

# Merge
command = "mergeBed -c 4,5,6 -o first,mean,first -i "+tempSortFileName+" > "+allCtcfAllPromFileName
os.system(command)

###################################################################################################
# LOOP + GENE
###################################################################################################

# Output files
allCtcfInLoopAllGeneFileName = outputPrefix+"allCtcfInLoopAllGene.bed"
command = "intersectBed -wa -u -a "+allCtcfAllGeneFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopAllGeneFileName
os.system(command)

allCtcfInLoopPlusGeneFileName = outputPrefix+"allCtcfInLoopPlusGene.bed"
command = "intersectBed -wa -u -a "+allCtcfPlusGeneFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopPlusGeneFileName
os.system(command)

allCtcfInLoopMinusGeneFileName = outputPrefix+"allCtcfInLoopMinusGene.bed"
command = "intersectBed -wa -u -a "+allCtcfMinusGeneFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopMinusGeneFileName
os.system(command)

allCtcfOutLoopAllGeneFileName = outputPrefix+"allCtcfOutLoopAllGene.bed"
command = "intersectBed -wa -u -a "+allCtcfAllGeneFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopAllGeneFileName
os.system(command)

allCtcfOutLoopPlusGeneFileName = outputPrefix+"allCtcfOutLoopPlusGene.bed"
command = "intersectBed -wa -u -a "+allCtcfPlusGeneFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopPlusGeneFileName
os.system(command)

allCtcfOutLoopMinusGeneFileName = outputPrefix+"allCtcfOutLoopMinusGene.bed"
command = "intersectBed -wa -u -a "+allCtcfMinusGeneFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopMinusGeneFileName
os.system(command)

###################################################################################################
# LOOP + PROM
###################################################################################################

# Output files
allCtcfInLoopAllPromFileName = outputPrefix+"allCtcfInLoopAllProm.bed"
command = "intersectBed -wa -u -a "+allCtcfAllPromFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopAllPromFileName
os.system(command)

allCtcfInLoopPlusPromFileName = outputPrefix+"allCtcfInLoopPlusProm.bed"
command = "intersectBed -wa -u -a "+allCtcfPlusPromFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopPlusPromFileName
os.system(command)

allCtcfInLoopMinusPromFileName = outputPrefix+"allCtcfInLoopMinusProm.bed"
command = "intersectBed -wa -u -a "+allCtcfMinusPromFileName+" -b "+allCtcfInLoopFileName+" > "+allCtcfInLoopMinusPromFileName
os.system(command)

allCtcfOutLoopAllPromFileName = outputPrefix+"allCtcfOutLoopAllProm.bed"
command = "intersectBed -wa -u -a "+allCtcfAllPromFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopAllPromFileName
os.system(command)

allCtcfOutLoopPlusPromFileName = outputPrefix+"allCtcfOutLoopPlusProm.bed"
command = "intersectBed -wa -u -a "+allCtcfPlusPromFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopPlusPromFileName
os.system(command)

allCtcfOutLoopMinusPromFileName = outputPrefix+"allCtcfOutLoopMinusProm.bed"
command = "intersectBed -wa -u -a "+allCtcfMinusPromFileName+" -b "+allCtcfOutLoopFileName+" > "+allCtcfOutLoopMinusPromFileName
os.system(command)

###################################################################################################
# ALL CTCF
###################################################################################################

# Output files
allCtcfFileName = outputPrefix+"allCtcf.bed"

# Cat
tempCatFileName = tempLoc+"tempCatFileName.bed"
command = "cat "+" ".join([allCtcfAllGeneFileName, allCtcfAllPromFileName, allCtcfInLoopFileName, allCtcfOutLoopFileName])+" > "+tempCatFileName
os.system(command)

# Sorting
tempSortFileName = tempLoc+"tempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+tempCatFileName+" > "+tempSortFileName
os.system(command)

# Merge
command = "mergeBed -c 4,5,6 -o first,mean,first -i "+tempSortFileName+" > "+allCtcfFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)



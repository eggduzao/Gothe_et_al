
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
nBins = int(sys.argv[1])
tssExt = int(sys.argv[2])
bamCount = int(sys.argv[3])
percentileList = sys.argv[4].split(",")
aliasFileName = sys.argv[5]
genesFileName = sys.argv[6]
featurePeakFileName = sys.argv[7]
bamFileName = sys.argv[8]
tempLocation = sys.argv[9]
outputFileName = sys.argv[10]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
outLoc = "/".join(outputFileName.split("/")[:-1])+"/"
command = "mkdir -p "+outLoc
os.system(command)
rpm = bamCount/1000000.

###################################################################################################
# Functions
###################################################################################################

def fetchTotalSignalBam(bamFile, region):
  totalSignal = 0.0
  for read in bamFile.fetch(region[0], region[1], region[2]): totalSignal += 1.0
  return totalSignal

def fetchSignalBw(bwFile, region, nBins, reverse = False):
  returnVec = [0.0] * nBins
  regionLen = (region[2] - region[1])
  bpPerBin = int( regionLen / nBins )
  valuesVec = bwFile.values(region[0], region[1], region[2])
  for i in range(0, len(valuesVec)):
    if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): valuesVec[i] = 0.0
  for i in range(0, regionLen, bpPerBin):
    returnVec[min(i/bpPerBin,len(returnVec)-1)] = sum(valuesVec[(i):(min(i+bpPerBin,len(valuesVec)-1))])
  return returnVec

def fetchSignalBam(bamFile, region, nBins, reverse = False):
  returnVec = [0.0] * nBins
  regionLen = (region[2] - region[1])
  bpPerBin = int( regionLen / nBins )
  correctionFactor = 200. / bpPerBin
  for read in bamFile.fetch(region[0], region[1], region[2]):
    returnVec[max(min(int((read.reference_start-region[1])/bpPerBin),nBins-1),0)] += (1.0 * correctionFactor)
  if(reverse): return returnVec[::-1]
  else: return returnVec

###################################################################################################
# Reading genes, lists and dictionaries
###################################################################################################

# Valid chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Reading alias
aliasDict = dict() # alias -> gene_symbol
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g] = value
aliasFile.close()

# Reading Tss
regionDict = dict() # gene_symbol -> [chr, tss-tssExt, tss, tss+tssExt, tts-tssExt, tts, tts+tssExt, gene, score, strand]
genesFile = open(genesFileName,"rU")
for line in genesFile:

  # Initialization
  ll = line.strip().split("\t")
  if(ll[2] not in chrList): continue
  try: gene = aliasDict[ll[12]]
  except Exception: continue
  try:
    eduardoehmuitolindo = regionDict[gene]
    continue
  except Exception: pass
  chrom = ll[2]; txStart = ll[4]; txEnd = ll[5]; score = ll[11]; strand = ll[3]

  # Region
  if(strand == "+"):
    p1 = int(txStart) - tssExt; p2 = int(txStart); p3 = int(txStart) + tssExt
    e1 = int(txEnd) - tssExt; e2 = int(txEnd); e3 = int(txEnd) + tssExt
  else:
    p1 = int(txEnd) + tssExt; p2 = int(txEnd); p3 = int(txEnd) - tssExt
    e1 = int(txStart) + tssExt; e2 = int(txStart); e3 = int(txStart) - tssExt

  # Region
  regionDict[gene] = [chrom, p1, p2, p3, e1, e2, e3, gene, score, strand]

genesFile.close()

###################################################################################################
# Expression stratification
###################################################################################################

# Fetching expression (feature) for all genes
featureList = []
featureDict = dict()
featurePeakFile = open(featurePeakFileName,"rU")
featurePeakFile.readline()
for line in featurePeakFile:
  ll = line.strip().split("\t")
  featureGene = ll[0]
  try:
    gene = aliasDict[featureGene]
    region = regionDict[gene]
  except Exception: continue
  featurevalue = float(ll[1]) / (max(region[2],region[5]) - min(region[2],region[5]))
  featureDict[gene] = featurevalue
  featureList.append(featurevalue)
featurePeakFile.close()

# Percentile values dictionary
percValueDict = dict()
featureListNp = np.array(featureList)
for percentile in percentileList: percValueDict[percentile] = np.percentile(featureListNp, int(percentile))

# Percentile dictionary
percentileDict = dict()
featureDictKeys = sorted(featureDict.keys())
for gene in featureDictKeys:
  for percentile in percentileList:
    percValue = percValueDict[percentile]
    featureValue = featureDict[gene]
    if(featureValue >= percValue):
      percentileDict[gene] = percentile
      break

###################################################################################################
# Creating table
###################################################################################################

# Initialization
featureKeysLen = len(featureDictKeys)
bamFile = Samfile(bamFileName, "rb")

# Fetching the bam signal in all categories
# GENE, GRO_VALUE, PERCENTILE, SIGNAL1, SIGNAL2, .....
outputFile = open(outputFileName,"w")
featureListNp = np.array(featureList)
prevPercValue = max(featureList)
for gene in featureDictKeys:

  # Initialization
  region = regionDict[gene]
  rev = True
  if(region[-1] == "+"): rev = False
  if(region[1] < 0): continue

  if(rev): totalSignal = fetchTotalSignalBam(bamFile, [region[0], region[5], region[2]])
  else: totalSignal = fetchTotalSignalBam(bamFile, [region[0], region[2], region[5]])
  vector = [gene, featureDict[gene] * totalSignal]

  # Verifying meta-gene
  if(rev and region[3]-region[4] < 2*nBins): continue
  if(not rev and region[4]-region[3] < 2*nBins): continue

  # Fetching percentile
  perc = percentileDict[gene]
  vector.append(perc)

  # Fetching signal
  if(rev):
    signal1 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[2], region[1]], nBins, reverse = rev)]
    signal2 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[3], region[2]], nBins, reverse = rev)]
    signal3 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[4], region[3]], 2*nBins, reverse = rev)]
    signal4 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[5], region[4]], nBins, reverse = rev)]
    signal5 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[6], region[5]], nBins, reverse = rev)]
  else:
    signal1 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[1], region[2]], nBins, reverse = rev)]
    signal2 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[2], region[3]], nBins, reverse = rev)]
    signal3 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[3], region[4]], 2*nBins, reverse = rev)]
    signal4 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[4], region[5]], nBins, reverse = rev)]
    signal5 = [e/rpm for e in fetchSignalBam(bamFile, [region[0], region[5], region[6]], nBins, reverse = rev)]

  # Updating vector
  vector = vector + signal1 + signal2 + signal3 + signal4 + signal5

  # Writing vector
  outputFile.write("\t".join([str(e) for e in vector])+"\n")

# Closing all files
outputFile.close()
bamFile.close()

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)



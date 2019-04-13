
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
percentileList = sys.argv[3].split(",")
aliasFileName = sys.argv[4]
geneFileName = sys.argv[5]
ctcfFileName = sys.argv[6]
groListFileName = sys.argv[7]
signalFileName = sys.argv[8]
outputFileName = sys.argv[9]

# Initialization
outLoc = "/".join(outputFileName.split("/")[:-1])+"/"
command = "mkdir -p "+outLoc
os.system(command)
bamExt = 5

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

###################################################################################################
# Expression stratification
###################################################################################################

# Fetching GRO for all genes
if(groListFileName.split(".")[-1] == "tsv"): # Sergi's file
  groList = []
  groDict = dict()
  groListFile = open(groListFileName, "rU")
  groListFile.readline()
  for line in groListFile:
    ll = line.strip().split("\t")
    groGene = ll[0].split("\"")[1]
    try: gene = aliasDict[groGene]
    except Exception: continue
    grovalue = float(ll[1])
    groDict[gene] = grovalue
    groList.append(grovalue)
  groListFile.close()
#######################################################
else: # Natasa's file
  groList = []
  groDict = dict()
  groListFile = open(groListFileName, "rU")
  groListFile.readline()
  for line in groListFile:
    ll = line.strip().split("\t")
    try:
      gene = aliasDict[ll[0]]
      gr = geneDict[gene]
    except Exception: continue
    region = [gr[0], int(gr[1]), int(gr[2])]
    groValue = float(ll[1]) / (float(gr[2]) - float(gr[1]))
    groDict[gene] = groValue
    groList.append(groValue)
  groListFile.close()

# Percentile values dictionary
percValueDict = dict()
groListNp = np.array(groList)
for percentile in percentileList: percValueDict[percentile] = np.percentile(groListNp, int(percentile))

# Percentile dictionary
percentileDict = dict()
groDictKeys = sorted(groDict.keys())
for gene in groDictKeys:
  for percentile in percentileList:
    percValue = percValueDict[percentile]
    groValue = groDict[gene]
    if(groValue >= percValue):
      percentileDict[gene] = percentile
      break

###################################################################################################
# Creating table
###################################################################################################

# Initialization
signalFile = Samfile(signalFileName, "rb")

# Fetching the bam signal in all categories
# GENE, GENE_CHR, GENE_P1, GENE_P2, GENE_STR, CTCF_CHR, CTCF_P1, CTCF_P2, CTCF_STR, GRO_VALUE, GRO_PERC, [SIGNAL...]
ctcfFile = open(ctcfFileName, "rU")
outputFile = open(outputFileName,"w")
for line in ctcfFile:

  # Initialization
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  ctcfChrom = ll[0]; ctcfP1 = int(ll[1]); ctcfP2 = int(ll[2]); ctcfGeneName = ll[3]; ctcfScore = ll[4]; ctcfStrand = ll[5]
  mid = (ctcfP1+ctcfP2)/2

  if(higherSide == "left"): region = [ctcfChrom, mid-ctcfExt, mid]
  elif(higherSide == "right"): region = [ctcfChrom, mid, mid+ctcfExt]
  totalSignal = fetchTotalSignalBam(signalFile, [ctcfChrom, mid-ctcfExt, mid+ctcfExt])

  # Fetching location, gene and gro
  try:
    geneName = aliasDict[ctcfGeneName]
    perc = percentileDict[geneName]
    gg = geneDict[geneName]
    vector = [geneName, gg[0], gg[1], gg[2], gg[5], ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, (float(groDict[geneName])+1) * totalSignal, int(perc) + totalSignal]
  except Exception: vector = ["NA", "NA", "NA", "NA", "NA", ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, 0, 0]

  # Fetching signal
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
  elif(higherSide == "restR"):
    r11p = 0.4
    r12p = 0.4
    r13p = 0.3
    r21p = 0.3
    r22p = 0.4
    r23p = 0.4
  elif(higherSide == "restx"):
    r11p = 0.45
    r12p = 0.32
    r13p = 0.32
    r21p = 0.37
    r22p = 0.41
    r23p = 0.45
  elif(higherSide == "resty"):
    r11p = 0.22
    r12p = 0.2
    r13p = 0.25
    r21p = 0.4
    r22p = 0.13
    r23p = 0.16
  elif(higherSide == "bothfr"):
    r11p = 0.6
    r12p = 0.6
    r13p = 0.6
    r21p = 0.6
    r22p = 0.6
    r23p = 0.6
  elif(higherSide == "bothf"):
    r11p = 0.63
    r12p = 0.63
    r13p = 0.63
    r21p = 0.63
    r22p = 0.63
    r23p = 0.63
  elif(higherSide == "bothr"):
    r11p = 0.57
    r12p = 0.57
    r13p = 0.57
    r21p = 0.57
    r22p = 0.57
    r23p = 0.57
  else:
    r11p = 0.6
    r12p = 0.6
    r13p = 0.6
    r21p = 0.6
    r22p = 0.6
    r23p = 0.6
  signal1 = [e*r11p for e in fetchSignalBam(signalFile, region11, bamExt)] + [e*r12p for e in fetchSignalBam(signalFile, region12, bamExt)] + [e*r13p for e in fetchSignalBam(signalFile, region13, bamExt)]
  signal2 = [e*r21p for e in fetchSignalBam(signalFile, region21, bamExt)] + [e*r22p for e in fetchSignalBam(signalFile, region22, bamExt)] + [e*r23p for e in fetchSignalBam(signalFile, region23, bamExt)]
  signal = signal1 + signal2

  # Updating vector
  vector = vector + signal

  # Writing vector
  outputFile.write("\t".join([str(e) for e in vector])+"\n")

# Closing all files
ctcfFile.close()
outputFile.close()
signalFile.close()



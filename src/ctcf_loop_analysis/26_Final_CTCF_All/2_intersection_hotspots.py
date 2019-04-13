
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
ctcfExt = int(sys.argv[1])
numberOfCtcf = int(sys.argv[2])
ctcfFileNameList = sys.argv[3].split(",")
hotspotFileName = sys.argv[4]
signalFileName = sys.argv[5]
tempLocation = sys.argv[6]
outputFileBarName = sys.argv[7]
outputFileVioName = sys.argv[8]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def fetchTotalSignalBam(bamFile, region):
  totalSignal = 0.0
  for read in bamFile.fetch(region[0], region[1], region[2]): totalSignal += 1.0
  return totalSignal

def fileLen(fileName):
  if(os.stat(fileName).st_size <= 0): return 0
  i = 0
  with open(fileName) as f:
    for i, l in enumerate(f): pass
  return i + 1

###################################################################################################
# Creating intersection files
###################################################################################################

# Initialization
intersectionLabelList = ["_".join(e.split("/")[-1].split(".")[0].split("_")[1:]) for e in ctcfFileNameList]
intersectionFileListH = [tempLocation+e+"_intH.bed" for e in intersectionLabelList]
intersectionFileListC = [tempLocation+e+"_intC.bed" for e in intersectionLabelList]
newCtcfFileNameList = [tempLocation+e+"_newCTCF.bed" for e in intersectionLabelList]

# Sorting hotspot file
newHotspotFileName = tempLocation+"newHotspotFileName.bed"
command = "sort -k1,1 -k2,2n "+hotspotFileName+" > "+newHotspotFileName
os.system(command)
numberOfHotspots = fileLen(newHotspotFileName)

# Iterating on ctcfFiles
for i in range(0,len(ctcfFileNameList)):

  # Files
  ctcfFileName = ctcfFileNameList[i]
  intersectionFileH = intersectionFileListH[i]
  intersectionFileC = intersectionFileListC[i]
  newCtcfFileName = newCtcfFileNameList[i]

  # Extending CTCF file
  extCtcfFileName = tempLocation+"extCtcfFileName.bed"
  extCtcfFile = open(tempLocation+"extCtcfFileName.bed", "w")
  ctcfFile = open(ctcfFileName, "r")
  for line in ctcfFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    mid = (p1+p2)/2
    extCtcfFile.write("\t".join([chrom, str(max(0,mid-ctcfExt)), str(mid+ctcfExt)])+"\n")
  ctcfFile.close()
  extCtcfFile.close()

  # Sorting ctcfFile
  command = "sort -k1,1 -k2,2n "+extCtcfFileName+" > "+newCtcfFileName
  os.system(command)

  # Performing intersection H
  command = "intersectBed -wa -u -a "+newHotspotFileName+" -b "+newCtcfFileName+" > "+intersectionFileH
  os.system(command)

  # Performing intersection C
  command = "intersectBed -wa -u -a "+newCtcfFileName+" -b "+newHotspotFileName+" > "+intersectionFileC
  os.system(command)

###################################################################################################
# Creating barplot output
###################################################################################################

# Output file
outputFileBar = open(outputFileBarName, "w")
header = ["X", "Y", "Z", "C"]
outputFileBar.write("\t".join(header)+"\n")

# Iterating on CTCF files
for i in range(0,len(newCtcfFileNameList)):

  # Files
  ctcfFileName = newCtcfFileNameList[i]
  intersectionLabel = intersectionLabelList[i]
  intersectionFileH = intersectionFileListH[i]
  intersectionFileC = intersectionFileListC[i]

  # Fetching numbers
  numberOfIntH = fileLen(intersectionFileH)
  numberOfIntC = fileLen(intersectionFileC)
  perc1 = round(100*float(numberOfIntH)/numberOfHotspots, 2)
  perc2 = round(100*float(numberOfIntC)/numberOfCtcf, 2)

  # Writing
  vec1 = [intersectionLabel, perc1, "HOTSPOT_WITH_CTCF", numberOfIntH]
  vec2 = [intersectionLabel, perc2, "CTCF_WITH_HOTSPOT", numberOfIntC]
  outputFileBar.write("\t".join([str(e) for e in vec1])+"\n")
  outputFileBar.write("\t".join([str(e) for e in vec2])+"\n")
  
# Closing files
outputFileBar.close()

###################################################################################################
# Creating violinplot output
###################################################################################################
"""
# Initialization
signalFile = Samfile(signalFileName, "rb")
header = ["X", "Y", "Z"]
outputFileVio = open(outputFileVioName, "w")
outputFileVio.write("\t".join(header)+"\n")

# Iterating on intersection hotspots
for i in range(0,len(intersectionFileListH)):

  # Files
  intersectionLabel = intersectionLabelList[i]
  intersectionFileName = intersectionFileListH[i]
  
  # Iterating on intersection file
  intersectionFile = open(intersectionFileName, "rU")
  for line in intersectionFile:
    ll = line.strip().split("\t")
    region = [ll[0], int(ll[1]), int(ll[2])]
    signal = fetchTotalSignalBam(signalFile, region)
    outputFileVio.write("\t".join([intersectionLabel, str(signal), "HOTSPOT_WITH_CTCF"])+"\n")
  intersectionFile.close()

# Interating on intersection CTCFs
for i in range(0,len(intersectionFileListC)):

  # Files
  intersectionLabel = intersectionLabelList[i]
  intersectionFileName = intersectionFileListC[i]
  
  # Iterating on intersection file
  intersectionFile = open(intersectionFileName, "rU")
  for line in intersectionFile:
    ll = line.strip().split("\t")
    region = [ll[0], int(ll[1]), int(ll[2])]
    signal = fetchTotalSignalBam(signalFile, region)
    outputFileVio.write("\t".join([intersectionLabel, str(signal), "CTCF_WITH_HOTSPOT"])+"\n")
  intersectionFile.close()

# Closing files
signalFile.close()
outputFileVio.close()
"""

###################################################################################################
# Termination
###################################################################################################

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)



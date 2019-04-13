
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input 
loopFileName = sys.argv[1]
geneFileName = sys.argv[2]
tempLoc = sys.argv[3]
outputFileName1 = sys.argv[4]
outputFileName2 = sys.argv[5]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

def intervalDistance(r1, r2):
  x, y = sorted((r1, r2))
  if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)): return y[0] - x[1]
  return 0

def getClosestGene(signalFile, region, ext = 10000):
  geneList = []
  geneFetch = signalFile.fetch(region[0], region[1] - ext, region[2] + ext)
  for read in geneFetch: geneList.append([read.qname, read.pos, read.aend])
  distance = 999999999
  gene = None
  for g in geneList:
    dist = intervalDistance([g[1], g[2]], [region[1], region[2]])
    if(dist < distance):
      distance = dist
      gene = g[0]
  return gene

def getBestMotifFileName(inFileName, outFileName):
  bestDict = dict()
  inFile = open(inFileName, "rU")
  for line in inFile:
    ll = line.strip().split("\t")
    gene = ll[3]
    score = int(ll[4])
    try:
      gg = bestDict[gene].strip().split("\t")
      if(score > int(gg[4])): bestDict[gene] = line
    except Exception: bestDict[gene] = line
  inFile.close()
  outFile = open(outFileName, "w")
  for k in bestDict.keys():
    outFile.write(bestDict[k])
  outFile.close()


# Execution
geneFile = Samfile(geneFileName, "rb")
loopFile = open(loopFileName, "rU")
loopFileName1 = tempLoc+"loopFile1.bed"
loopFileName2 = tempLoc+"loopFile2.bed"
loopFile1 = open(loopFileName1, "w")
loopFile2 = open(loopFileName2, "w")
loopFile.readline()
for line in loopFile:

  # Initialization
  ll = line.strip().split("\t")
  chrom = "chr"+ll[0]
  if(chrom not in chrList): continue
  ctcfx1 = ll[20]; ctcfx2 = ll[21]; ctcfxo = ll[23]
  ctcfy1 = ll[25]; ctcfy2 = ll[26]; ctcfyo = ll[28]
  loopScore = ll[7]

  # Closest gene 1
  if(ctcfxo != "NA"):
    if(ctcfxo == "p"): strand = "+"
    else: strand = "+"
    closestGene = getClosestGene(geneFile, [chrom, int(ctcfx1), int(ctcfx2)])
    if(closestGene): loopFile1.write("\t".join([chrom, ctcfx1, ctcfx2, closestGene.upper(), loopScore, strand])+"\n")

  # Closest gene 2
  if(ctcfyo != "NA"):
    if(ctcfyo == "p"): strand = "+"
    else: strand = "+"
    closestGene = getClosestGene(geneFile, [chrom, int(ctcfy1), int(ctcfy2)])
    if(closestGene): loopFile2.write("\t".join([chrom, ctcfy1, ctcfy2, closestGene.upper(), loopScore, strand])+"\n")

geneFile.close()
loopFile.close()
loopFile1.close()
loopFile2.close()

# Fetching best motifs
bestFileName1 = tempLoc+"bestFileName1.bed"
bestFileName2 = tempLoc+"bestFileName2.bed"
getBestMotifFileName(loopFileName1, bestFileName1)
getBestMotifFileName(loopFileName2, bestFileName2)

# Sorting bed file
command = "sort -k1,1 -k2,2n "+bestFileName1+" > "+outputFileName1
os.system(command)
command = "sort -k1,1 -k2,2n "+bestFileName2+" > "+outputFileName2
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)



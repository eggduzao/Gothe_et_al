
# Import
import os
import sys

# Input 
inFileName1 = sys.argv[1]
inFileName2 = sys.argv[2]
tempLoc = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
command = "mkdir -p "+tempLoc
os.system(command)

# Cat files
catFileName = tempLoc+"catFileName.bed"
command = "cat "+inFileName1+" "+inFileName2+" > "+catFileName
os.system(command)

# Grep chromosomes
grepFileName = tempLoc+"grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+catFileName+" > "+grepFileName
os.system(command)

# Get only best motifs
bestDict = dict() # GENE -> line
grepFile = open(grepFileName, "rU")
for line in grepFile:
  ll = line.strip().split("\t")
  gene = ll[3]
  score = int(ll[4])
  try:
    gg = bestDict[gene].strip().split("\t")
    if(score > int(gg[4])): bestDict[gene] = line
  except Exception: bestDict[gene] = line
grepFile.close()
bestMotifFileName = tempLoc+"bestMotifFileName.bed"
bestMotifFile = open(bestMotifFileName, "w")
for k in bestDict.keys():
  bestMotifFile.write(bestDict[k])
bestMotifFile.close()

# Get unique entries
uniqFileName = tempLoc+"uniqFileName.bed"
command = "sort "+bestMotifFileName+" | uniq > "+uniqFileName
os.system(command)

# Sort bed
command = "sort -k1,1 -k2,2n "+uniqFileName+" > "+outputFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLoc
os.system(command)



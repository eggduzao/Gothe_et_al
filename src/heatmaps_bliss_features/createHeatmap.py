
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
ext = int(sys.argv[1])
yMax = int(sys.argv[2])
featureSummitFileName = sys.argv[3]
bwFileName = sys.argv[4]
bwLabel = sys.argv[5]
tempLocation = sys.argv[6]
outputLocation = sys.argv[7]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Heatmaps
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Creating bed from peak files
featureSummitFile = open(featureSummitFileName,"r")
tempBedFileName = tempLocation+"bedfile.bed"
bedFile = open(tempBedFileName,"w")
for line in featureSummitFile:
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  region = [ll[0], int(ll[1])-ext, int(ll[2])+ext]
  if(int(region[1]) < 0): continue
  bedFile.write("\t".join([str(e) for e in region])+"\n")
featureSummitFile.close()
bedFile.close()

# Creating heatmaps
tempMatFileName = tempLocation+"matrix.mat.gz"

# Creating matrix (heatmaps_eachSignalOrderByDBS)
command = "computeMatrix reference-point -S \""+bwFileName+"\" -R \""+tempBedFileName+"\" -a \""+str(ext)+"\" -b \""+str(ext)+"\" --referencePoint \"center\" --binSize \"10\" --sortRegions \"keep\" --missingDataAsZero --smartLabels --numberOfProcessors \"max/2\" -o \""+tempMatFileName+"\""
os.system(command)

# Creating heatmap with yMax
#command = "plotHeatmap -m \""+tempMatFileName+"\" -out \""+outputLocation+bwLabel+".pdf\" --dpi \"90\" --missingDataColor \"white\" --refPointLabel \"Summit\" --yAxisLabel \""+bwLabel+" Signal\" --yMax "+str(yMax)+" --samplesLabel \""+bwLabel+"\" --legendLocation \"upper-right\" --plotFileFormat \"pdf\""
#os.system(command)

# Creating heatmap without yMax
command = "plotHeatmap -m \""+tempMatFileName+"\" -out \""+outputLocation+bwLabel+".pdf\" --dpi \"90\" --missingDataColor \"white\" --refPointLabel \"Summit\" --yAxisLabel \""+bwLabel+" Signal\" --samplesLabel \""+bwLabel+"\" --legendLocation \"upper-right\" --plotFileFormat \"pdf\""
os.system(command)
# --sortRegions \"keep\"

# Termination
command = "rm "+" ".join([tempBedFileName, tempMatFileName])
os.system(command)



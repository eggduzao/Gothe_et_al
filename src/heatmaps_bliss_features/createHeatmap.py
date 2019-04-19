
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

###################################################################################################
# Create Heatmap
###################################################################################################

def create_heatmap(half_ext, feature_summit_file_name, signal_file_name, signal_label, temp_location, output_file_name):

  # Initialization
  command = "mkdir -p "+temp_location
  os.system(command)
  outLoc = "/".join(output_file_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Allowed chromosomes
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Creating bed from peak files
  featureSummitFile = open(feature_summit_file_name,"r")
  tempBedFileName = temp_location + "bedfile.bed"
  bedFile = open(tempBedFileName,"w")
  for line in featureSummitFile:
    ll = line.strip().split("\t")
    if(ll[0] not in chrList): continue
    region = [ll[0], int(ll[1])-half_ext, int(ll[2])+half_ext]
    if(int(region[1]) < 0): continue
    bedFile.write("\t".join([str(e) for e in region])+"\n")
  featureSummitFile.close()
  bedFile.close()

  # Creating heatmaps
  tempMatFileName = temp_location + "matrix.mat.gz"

  # Creating matrix (heatmaps_eachSignalOrderByDBS)
  command = "computeMatrix reference-point -S \""+signal_file_name+"\" -R \""+tempBedFileName+"\" -a \""+str(half_ext)+"\" -b \""+str(half_ext)+"\" --referencePoint \"center\" --binSize \"10\" --sortRegions \"keep\" --missingDataAsZero --numberOfProcessors \"max/2\" -o \""+tempMatFileName+"\""
  os.system(command)

  # Creating heatmap without yMax
  command = "plotHeatmap -m \""+tempMatFileName+"\" -out \""+output_file_name+"\" --dpi \"90\" --missingDataColor \"white\" --refPointLabel \"Summit\" --yAxisLabel \""+signal_label+" Signal\" --samplesLabel \""+signal_label+"\" --legendLocation \"upper-right\" --plotFileFormat \"pdf\""
  os.system(command)


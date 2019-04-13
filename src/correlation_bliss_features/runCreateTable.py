
# Import
import os
import sys

# Input
inputNameListFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/input/name_list.txt"
inputCountListFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/input/count_list.txt"
inputBamListFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/input/bam_list.txt"
inputBwListFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/input/bw_list.txt"

# Creating nameDict
nameDict = dict()
name = None
inputNameListFile = open(inputNameListFileName,"r")
for line in inputNameListFile:
  if(line[0] == "#"):
    name = line.strip()[2:]
    nameDict[name] = []
    continue
  ll = line.strip().split("\t")
  nameDict[name].append(ll[2])
inputNameListFile.close()

# Creating countDict
countDict = dict()
name = None
inputCountListFile = open(inputCountListFileName,"r")
for line in inputCountListFile:
  if(line[0] == "#"):
    name = line.strip()[2:]
    countDict[name] = []
    continue
  ll = line.strip().split("\t")
  countDict[name].append(ll[1])
inputCountListFile.close()

# Creating bamDict
bamDict = dict()
name = None
inputBamListFile = open(inputBamListFileName,"r")
for line in inputBamListFile:
  if(line[0] == "#"):
    name = line.strip()[2:]
    bamDict[name] = []
    continue
  ll = line.strip()
  if(ll == "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/mnase/wgEncodeSydhNsomeK562AlnRep1.bam"):
    bamDict[name].append("/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/mnase/new_bw/K562_MNAse.bw")
  else:
    bamDict[name].append(ll)
inputBamListFile.close()

# Iterating on the categories
inputList = sorted(bamDict.keys())
for k in inputList:

  if(k != "SelectedFactors" and k != "SelectedHistones"): continue
  
  # Parameters
  ext = "500"
  featureSummitFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/results/macs2/K562_Etop20uM_rep2.3.deep_summits.bed"
  bamNames = ",".join(nameDict[k])
  bamCounts = ",".join(countDict[k])
  bamList = ",".join(bamDict[k])
  outputFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/correlation_selected/"+k+".txt"

  # Execution
  command = "python createTable.py "+" ".join([ext, featureSummitFileName, bamNames, bamCounts, bamList, outputFileName])
  os.system(command)



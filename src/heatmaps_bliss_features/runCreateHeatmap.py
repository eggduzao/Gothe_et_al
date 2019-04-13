
# Import
import os
import sys

# Input
inputNameListFileName = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/input/name_list.txt"
inputBwListFileName = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/input/bw_list.txt"

# Y Limit
yMaxDict = {"Broad": "180", "Haib": "18", "Histones": "70", "Occ": "80", "Sydh": "130", "Uc": "30", "Bliss": "20", "SelectedFactors": "130", "SelectedHistones": "70"}

# Creating nameList
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

# Creating bwList
bwDict = dict()
name = None
inputBwListFile = open(inputBwListFileName,"r")
for line in inputBwListFile:
  if(line[0] == "#"):
    name = line.strip()[2:]
    bwDict[name] = []
    continue
  ll = line.strip()
  bwDict[name].append(ll)
inputBwListFile.close()

# Iterating on the categories
inputList = sorted(bwDict.keys())
for k in inputList:

  if(k != "SelectedFactors" and k != "SelectedHistones"): continue

  bwList = bwDict[k]
  nameList = nameDict[k]
  yMax = yMaxDict[k]

  for i in range(0,len(bwList)):
  
    # Parameters
    ext = "1500"
    featureSummitFileName = "/media/egg/sbc/AG_Papantonis/eduardo/Roukos_BLISS/data/BLISS2_corrected/results/macs2/K562_Etop20uM_rep2.3.deep_summits.bed"
    bwFileName = bwList[i]
    bwLabel = nameList[i]
    tempLocation = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/heatmap_selected/TEMP/"+k+"/"
    outputLocation = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/heatmap_selected/"+k+"/"

    # Execution
    command = "python createHeatmap.py "+" ".join([ext, yMax, featureSummitFileName, bwFileName, bwLabel, tempLocation, outputLocation])
    os.system(command)




# Import
import os
import sys

# Input
#il = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/correlation/"
#inputList = [il+"Bliss.txt", il+"Broad.txt", il+"Haib.txt", il+"Histones.txt", il+"Occ.txt", il+"Sydh.txt", il+"Uc.txt"]
#widList = ["4", "20", "25", "15", "8", "45", "8"]
#mList = ["9.5", "9.5", "9.5", "9.5", "9.5", "16", "9.5"]

# Selected input
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/correlation_selected/"
inputList = [il+"SelectedFactors.txt"]
widList = ["10"]
mList = ["9.5"]

# Selected input
#il = "/home/egg/Projects/Roukos_Bliss/Results/3_heatmap_and_spearman/correlation_selected/"
#inputList = [il+"SelectedHistones.txt"]
#widList = ["13"]
#mList = ["9.5"]

for i in range(0,len(inputList)):

  # Parameters
  graphWidth = widList[i]
  marginX = mList[i]
  blissTableFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/3_heatmap_and_spearman/correlation/Bliss.txt"
  inputTableFileName = inputList[i]
  outputFileName = inputTableFileName[:-3]+"pdf"
  outputLocation = il+"correlation_plots/"

  # Execution
  command = "R CMD BATCH '--args '"+graphWidth+"' '"+marginX+"' '"+blissTableFileName+"' '"+inputTableFileName+"' '"+outputFileName+"' '"+outputLocation+" correlation.R correlation.Rout"
  os.system(command)



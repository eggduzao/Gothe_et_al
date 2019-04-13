
# Import
import os
import sys

# Cell List
ll = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/4_Extended_Anchors/"
cellList = ["K562", "TK6", "CD34"]

# Cell Loop
for cell in cellList:

  # Execution anchors_with_ctcf
  aliasFileName = "/Users/egg/rgtdata/hg19/alias_human_booster.txt"
  allGenesFileName = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/all_genes_and_mll.bed"
  anchorAllFilePrefix = ll + cell + "/anchors_with_ctcf"
  outputFileName = ll + cell + "_anchors_with_ctcf.txt"
  print cell, "anchors_with_ctcf"
  command = "python 4_createDistanceTable.py "+" ".join([aliasFileName, allGenesFileName, anchorAllFilePrefix, outputFileName])
  os.system(command)

  # Execution anchors_with_and_wo_ctcf
  aliasFileName = "/Users/egg/rgtdata/hg19/alias_human_booster.txt"
  allGenesFileName = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/all_genes_and_mll.bed"
  anchorAllFilePrefix = ll + cell + "/anchors_with_and_wo_ctcf"
  outputFileName = ll + cell + "_anchors_with_and_wo_ctcf.txt"
  print cell, "anchors_with_and_wo_ctcf"
  command = "python 4_createDistanceTable.py "+" ".join([aliasFileName, allGenesFileName, anchorAllFilePrefix, outputFileName])
  os.system(command)


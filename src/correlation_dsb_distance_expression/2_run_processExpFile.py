
# Import
import os
import sys

# Exp List
aliasFile = "/home/egg/rgtdata/hg19/alias_human_booster.txt"
chromSizesFile = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
geneLocationFile = "/home/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/all_genes.bed"
il = "/home/egg/Projects/Roukos_Bliss/Data/CD34/cd34/exp/"
ol = "/home/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/2_Exp/"
nameList = ["CD34_expression_readcounts"]
labelList = ["CD34_expression_filtered"]

# Exp Loop
for i in range(0,len(nameList)):

  # Input
  aliasFileName = aliasFile
  chromSizesFileName = chromSizesFile
  geneLocationFileName = geneLocationFile
  expFileName = il + nameList[i] + ".txt"
  outputFileName = ol + labelList[i] + ".txt"

  # Execution
  command = "python 2_processExpFile.py "+" ".join([aliasFileName, chromSizesFileName, geneLocationFileName, expFileName, outputFileName])
  os.system(command)



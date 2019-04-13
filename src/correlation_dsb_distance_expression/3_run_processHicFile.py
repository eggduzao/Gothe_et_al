
# Import
import os
import sys

# Hic List
chromSizesFile = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
ctcfPeaksFile = "/home/egg/Projects/Roukos_Bliss/Data/CD34/ctcf/CTCF_peaks.bam"
ctcfMotifsFile = "/home/egg/Projects/Roukos_Bliss/Data/CTCF_motifs.bam"
il = "/home/egg/Projects/Roukos_Bliss/Data/CD34/cd34/hic/"
ol = "/home/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/3_Hic/"
nameList = ["CD34_150Kb_10Mb"]
labelList = ["CD34_loops"]

# Hic Loop
for i in range(0,len(nameList)):

  # Input
  chromSizesFileName = chromSizesFile
  ctcfBamPeaksFileName = ctcfPeaksFile
  ctcfBamMotifsFileName = ctcfMotifsFile
  loopsFileName = il + nameList[i] + ".txt"
  loopsHiccupsOutFileName = ol + labelList[i] + ".txt"

  # Execution
  command = "python 3_processHicFile.py "+" ".join([chromSizesFileName, ctcfBamPeaksFileName, ctcfBamMotifsFileName, loopsFileName, loopsHiccupsOutFileName])
  os.system(command)



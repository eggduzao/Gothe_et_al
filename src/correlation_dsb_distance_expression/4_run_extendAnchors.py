
# Import
import os
import sys

# Cell List
chromFile = "/Users/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
il = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/3_Hic/"
ol = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/4_Extended_Anchors/"
cellList = ["K562", "TK6", "CD34"]

# Cell Loop
for cell in cellList:

  # Ext List
  extList = [str(e*1000) for e in range(0, 1001)]
  extList = ["249000"]

  # Ext Loop
  for ext in extList:

    # Input
    largestLengthHalf = ext
    loopFileName = il + cell + "_loops.txt"
    chromSizesFileName = chromFile
    tempLoc = "./TEMP/"
    outputFileWithAndWoCtcfName = ol + cell + "/anchors_with_and_wo_ctcf_" + ext
    outputFileWithCtcfName = ol + cell + "/anchors_with_ctcf_" + ext
    outputFileWoCtcfName = ol + cell + "/anchors_wo_ctcf_" + ext

    # Execution
    command = "python 4_extendAnchors.py "+" ".join([largestLengthHalf, loopFileName, chromSizesFileName, tempLoc, outputFileWithAndWoCtcfName, outputFileWithCtcfName, outputFileWoCtcfName])
    os.system(command)

/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/4_Extended_Anchors/K562/anchors_with_and_wo_ctcf_249000.bam
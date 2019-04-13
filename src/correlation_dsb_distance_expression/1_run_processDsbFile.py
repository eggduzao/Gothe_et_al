
# Import
import os
import sys

# Rep List
chromSizesFile = "/home/egg/rgtdata/hg19/chrom.sizes.hg19.filter"
il = "/home/egg/Projects/Roukos_Bliss/Data/CD34/cd34/dsb/"
tl = "./TEMP/"
ol = "/home/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/1_Dsb/"
repList = [["BB187_CD34pDMSOmixedrp1_GCTTGGAG_chr-loc-countDifferentUMI", "BB188_CD34pDMSOmixedrp2_GCTTGGAG_chr-loc-countDifferentUMI"], ["BB187_CD34pETOmixedrp1_TAGGTGAG_chr-loc-countDifferentUMI", "BB188_CD34pETOmixedrp2_TAGGTGAG_chr-loc-countDifferentUMI"]]
labelList = ["CD34_DMSO", "CD34_ETO"]

# Cell Loop
for i in range(0,len(repList)):

  # Input
  chromSizesFileName = chromSizesFile
  dsbBedFileNameList = ",".join([il + e + ".bed" for e in repList[i]])
  tempLocation = tl + str(i) + "/"
  dsbBamFileName = ol + labelList[i] + ".bam"

  # Execution
  command = "python 1_processDsbFile.py "+" ".join([chromSizesFileName, dsbBedFileNameList, tempLocation, dsbBamFileName])
  os.system(command)



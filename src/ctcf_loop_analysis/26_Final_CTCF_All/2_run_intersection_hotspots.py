
# Import
import os
import sys
from glob import glob
from pysam import Samfile

# Folder List
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/0_input/intersection/"
olb = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/2_hotspotBarplot/"
olv = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/3_hotspotViolinplot/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
hl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/results/macs2/"
folderList = ["all", "a", "b"]

# Folder Loop
for folder in folderList:

  # Subfolder List
  subfolderList = ["active_genes", "all_genes", "inactive_genes", "intergenic", "active_promoters", "all_promoters", "inactive_promoters"]

  # Subfolder Loop
  for subfolder in subfolderList:

    # Cell List
    cellList = ["K562", "TK6"]

    # Cell Loop
    for cell in cellList:
 
      # Number of CTCF
      if(cell == "K562"): nCtcf = "32913"
      elif(cell == "TK6"): nCtcf = "45439"

      # Files
      dsbFile = dl+"TK6_wt_ETO.bam"
      hotspotFile = hl+"TK6_WT_ETO_filtered.narrowPeak"
      if(cell == "K562"):
        dsbFile = dl+"K562_Etop20uM_rep2.3.deep.bam"
        hotspotFile = hl+"K562_Etop20uM_rep2.3.deep_peaks.narrowPeak"
      inLoc = il+folder+"/"+subfolder+"/"
      outLocBar = olb+folder+"/"
      outLocVio = olv+folder+"/"
      os.system("mkdir -p "+outLocBar+" "+outLocVio)

      # Input
      ctcfExt = "10000"
      numberOfCtcf = nCtcf
      ctcfFileNameList = ",".join(glob(inLoc+cell+"_*.bed"))
      hotspotFileName = hotspotFile
      signalFileName = dsbFile
      tempLocation = inLoc+cell+"/TEMP/"
      outputFileBarName = outLocBar+cell+"_"+subfolder+".txt"
      outputFileVioName = outLocVio+cell+"_"+subfolder+".txt"

      # Execution
      print folder, subfolder, cell
      command = "python 2_intersection_hotspots.py "+" ".join([ctcfExt, numberOfCtcf, ctcfFileNameList, hotspotFileName, signalFileName, tempLocation, outputFileBarName, outputFileVioName])
      os.system(command)



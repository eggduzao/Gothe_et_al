
# Import
import os
import sys

# Cell List
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/0_input/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/3_ctcfMotifsInsideLoops/"
cellList = ["K562", "GM12878"]

# Cell Loop
for cell in cellList:

  # Feature List
  fl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/1_ctcfMotifsInsideTranscription/"
  featureList = ["ctcf_genes", "ctcf_promoters"]

  # Feature Loop
  for feature in featureList:

    # Parameters
    if(feature == "ctcf_genes"): name = cell+"_ctcf_loop_genes"
    elif(feature == "ctcf_promoters"): name = cell+"_ctcf_loop_promoters"

    # Input
    ctcfBamFileName = fl+cell+"_"+feature+"_MFR_TFR.bam"
    ctcfBedFileName = fl+cell+"_"+feature+"_MFR_TFR.bed"
    loopFileName = il+cell+"_loops_full.txt"
    tempLoc = ol+"TEMP/"
    outputPrefix = ol+name

    # Execution
    print name
    command = "python 3_ctcfMotifsInsideLoops.py "+" ".join([ctcfBamFileName, ctcfBedFileName, loopFileName, tempLoc, outputPrefix])
    os.system(command)



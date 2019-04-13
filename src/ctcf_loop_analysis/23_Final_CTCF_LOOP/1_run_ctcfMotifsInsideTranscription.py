
# Import
import os
import sys

# Cell List
cellList = ["K562", "GM12878"]

# Cell Loop
for cell in cellList:

  # Feature List
  fl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/0_input/"
  ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/1_ctcfMotifsInsideTranscription/"
  featureList = ["all_genes", "all_promoters"]

  # Feature Loop
  for feature in featureList:

    # Parameters
    if(feature == "all_genes"): name = cell+"_ctcf_genes"
    elif(feature == "all_promoters"): name = cell+"_ctcf_promoters"

    # Input
    ctcfBamFileName = fl+cell+"_CTCF_active_motifs.bam"
    ctcfBedFileName = fl+cell+"_CTCF_active_motifs.bed"
    featureFileName = fl+feature+".bed"
    tempLoc = ol+"TEMP/"
    outputPrefix = ol+name

    # Execution
    print name
    command = "python 1_ctcfMotifsInsideTranscription.py "+" ".join([ctcfBamFileName, ctcfBedFileName, featureFileName, tempLoc, outputPrefix])
    os.system(command)



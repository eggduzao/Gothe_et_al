
# Import
import os
import sys

# Cell List
cellList = ["K562", "GM12878"]

# Cell Loop
for cell in cellList:

  # Feature List
  fl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/0_input/"
  il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/"
  ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/9_ctcfMotifsInsideInactive/"
  featureList = ["inactive_genes", "inactive_promoters"]

  # Feature Loop
  for feature in featureList:

    # Parameters
    if(feature == "inactive_genes"): name = cell+"_ctcf_inactive_genes"
    elif(feature == "inactive_promoters"): name = cell+"_ctcf_inactive_promoters"

    # Input
    ctcfBamFileName = fl+cell+"_CTCF_active_motifs.bam"
    ctcfBedFileName = fl+cell+"_CTCF_active_motifs.bed"
    featureFileName = il+cell+"/"+feature+".bed"
    tempLoc = ol+"TEMP/"
    outputPrefix = ol+name

    # Execution
    print name
    command = "python 9_ctcfMotifsInsideInactive.py "+" ".join([ctcfBamFileName, ctcfBedFileName, featureFileName, tempLoc, outputPrefix])
    os.system(command)



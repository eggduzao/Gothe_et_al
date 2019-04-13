
# Import
import os
import sys

# Cell List
fl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/0_input/"
gl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/5_CTCFatAllFeatures/"
cellList = ["K562", "GM12878"]

# Cell Loop
for cell in cellList:

  # Input 
  ctcfBamFileName = fl+cell+"_CTCF_active_motifs.bam"
  ctcfBedFileName = fl+cell+"_CTCF_active_motifs.bed"
  loopBedFileName = fl+cell+"_loops_full.txt"
  activeGeneBedFileName = gl+cell+"/active_genes.bed"
  inactiveGeneBedFileName = gl+cell+"/inactive_genes.bed"
  activePromBedFileName = gl+cell+"/active_promoters.bed"
  inactivePromBedFileName = gl+cell+"/inactive_promoters.bed"
  intergenicBedFileName = gl+cell+"/intergenic.bed"
  tempLoc = ol+"TEMP/"
  outputPrefix = ol+cell+"_"

  # Execution
  print cell
  command = "python 5_CTCFatAllFeatures.py "+" ".join([ctcfBamFileName, ctcfBedFileName, loopBedFileName, activeGeneBedFileName, inactiveGeneBedFileName, activePromBedFileName, inactivePromBedFileName, intergenicBedFileName, tempLoc, outputPrefix])
  os.system(command)



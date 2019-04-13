
# Import
import os
import sys

# Cell List
ll = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/loops/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/3_ctcf_motifs_loops/"
cellList = ["TK6", "K562"]

# Cell Loop
for cell in cellList:

  # Input 
  loopFileName = ll+cell+"_CTCF_loops.txt"
  geneFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/all_genes.bam"
  tempLoc = ol+"TEMP/"
  outputFileName1 = ol+cell+"_CTCF_loops_L1.txt"
  outputFileName2 = ol+cell+"_CTCF_loops_L2.txt"

  # Execution
  command = "python 3_assign_loop_gene.py "+" ".join([loopFileName, geneFileName, tempLoc, outputFileName1, outputFileName2])
  os.system(command)




# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/10_inactivePlot/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_ETO.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  # Region List
  cell = dsbFile.split("/")[-1].split("_")[0].upper()
  if(cell == "TK6"): cellD = "GM12878"
  else: cellD = "K562"
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/9_ctcfMotifsInsideInactive/"
  regionList = [cellD+"_ctcf_inactive_genes_MFR_TF", cellD+"_ctcf_inactive_genes_MFR_TFR", cellD+"_ctcf_inactive_genes_MFR_TR", cellD+"_ctcf_inactive_genes_MFTF_MRTR", cellD+"_ctcf_inactive_genes_MF_TFR", cellD+"_ctcf_inactive_genes_MFTR_MRTF", cellD+"_ctcf_inactive_genes_MR_TFR", cellD+"_ctcf_inactive_promoters_MFR_TF", cellD+"_ctcf_inactive_promoters_MFR_TFR", cellD+"_ctcf_inactive_promoters_MFR_TR", cellD+"_ctcf_inactive_promoters_MFTF_MRTR", cellD+"_ctcf_inactive_promoters_MF_TFR", cellD+"_ctcf_inactive_promoters_MFTR_MRTF", cellD+"_ctcf_inactive_promoters_MR_TFR"]

  # Region Loop
  for region in regionList:

    if("_TFR" in region): regionO = "both"
    elif("_TF" in region): regionO = "left"
    elif("_TR" in region): regionO = "right"
    else: regionO = "both"

    # Parameters
    if("DMSO" in dsbFile): cond = "DMSO"
    else: cond = "ETO"
    if(cell == "TK6"): groFile = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/external/results/counts/GM12878_GROseq.tsv"
    else: groFile = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/external/results/counts/K562_GROseq.tsv"
    name = "_".join([cell, cond]+region.split("_")[1:])

    # Input
    higherSide = regionO
    ctcfExt = "500"
    percentileList = ",".join([str(e) for e in range(100,-1,-5)])
    aliasFileName = "/home/egg/rgtdata/hg19/alias_human.txt"
    geneFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/"+cellD+"/all_genes.bed"
    ctcfFileName = rl+region+".bed"
    groListFileName = groFile
    signalFileName = dsbFile
    outputFileName = ol+"/"+name+".txt"

    # Execution
    print name
    command = "python 10_inactivePlot.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)



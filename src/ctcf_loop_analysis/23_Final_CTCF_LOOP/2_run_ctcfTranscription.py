
# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/2_CtcfTranscription/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_DMSO.bam", dl+"TK6_wt_ETO.bam", dl+"K562_DMSO_rep2.3.deep.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  # Region List
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/1_ctcfMotifsInsideTranscription/"
  regionList = ["GM12878_ctcf_genes_MFR_TF", "GM12878_ctcf_genes_MFR_TFR", "GM12878_ctcf_genes_MFR_TR", "GM12878_ctcf_genes_MFTF_MRTR", "GM12878_ctcf_genes_MF_TFR", "GM12878_ctcf_genes_MFTR_MRTF", "GM12878_ctcf_genes_MR_TFR", "GM12878_ctcf_genes_restMF", "GM12878_ctcf_genes_restMR", "GM12878_ctcf_promoters_MFR_TF", "GM12878_ctcf_promoters_MFR_TFR", "GM12878_ctcf_promoters_MFR_TR", "GM12878_ctcf_promoters_MFTF_MRTR", "GM12878_ctcf_promoters_MF_TFR", "GM12878_ctcf_promoters_MFTR_MRTF", "GM12878_ctcf_promoters_MR_TFR", "GM12878_ctcf_promoters_restMF", "GM12878_ctcf_promoters_restMR", "K562_ctcf_genes_MFR_TF", "K562_ctcf_genes_MFR_TFR", "K562_ctcf_genes_MFR_TR", "K562_ctcf_genes_MFTF_MRTR", "K562_ctcf_genes_MF_TFR", "K562_ctcf_genes_MFTR_MRTF", "K562_ctcf_genes_MR_TFR", "K562_ctcf_genes_restMF", "K562_ctcf_genes_restMR", "K562_ctcf_promoters_MFR_TF", "K562_ctcf_promoters_MFR_TFR", "K562_ctcf_promoters_MFR_TR", "K562_ctcf_promoters_MFTF_MRTR", "K562_ctcf_promoters_MF_TFR", "K562_ctcf_promoters_MFTR_MRTF", "K562_ctcf_promoters_MR_TFR", "K562_ctcf_promoters_restMF", "K562_ctcf_promoters_restMR"]

  # Region Loop
  for region in regionList:

    if("_TFR" in region): regionO = "both"
    elif("rest" in region): regionO = "rest"
    elif("_TF" in region): regionO = "left"
    elif("_TR" in region): regionO = "right"
    else: regionO = "both"

    # Parameters
    cell = dsbFile.split("/")[-1].split("_")[0].upper()
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
    geneFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/all_genes.bed"
    ctcfFileName = rl+region+".bed"
    groListFileName = groFile
    signalFileName = dsbFile
    outputFileName = ol+"/"+name+".txt"

    # Execution
    print name
    command = "python 2_ctcfTranscription.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)



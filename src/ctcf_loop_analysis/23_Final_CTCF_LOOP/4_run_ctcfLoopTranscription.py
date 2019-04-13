
# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/4_CtcfLoopTranscription/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_DMSO.bam", dl+"TK6_wt_ETO.bam", dl+"K562_DMSO_rep2.3.deep.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  # Region List
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/3_ctcfMotifsInsideLoops/"
  regionList = ["GM12878_ctcf_loop_genes_rest", "GM12878_ctcf_loop_genes_restIN", "GM12878_ctcf_loop_genes_restOUT", "GM12878_ctcf_loop_genes_TF", "GM12878_ctcf_loop_genes_TFR", "GM12878_ctcf_loop_genes_TI", "GM12878_ctcf_loop_genes_TIO", "GM12878_ctcf_loop_genes_TO", "GM12878_ctcf_loop_genes_TR", "GM12878_ctcf_loop_promoters_rest", "GM12878_ctcf_loop_promoters_restIN", "GM12878_ctcf_loop_promoters_restOUT", "GM12878_ctcf_loop_promoters_TF", "GM12878_ctcf_loop_promoters_TFR", "GM12878_ctcf_loop_promoters_TI", "GM12878_ctcf_loop_promoters_TIO", "GM12878_ctcf_loop_promoters_TO", "GM12878_ctcf_loop_promoters_TR", "K562_ctcf_loop_genes_rest", "K562_ctcf_loop_genes_restIN", "K562_ctcf_loop_genes_restOUT", "K562_ctcf_loop_genes_TF", "K562_ctcf_loop_genes_TFR", "K562_ctcf_loop_genes_TI", "K562_ctcf_loop_genes_TIO", "K562_ctcf_loop_genes_TO", "K562_ctcf_loop_genes_TR", "K562_ctcf_loop_promoters_rest", "K562_ctcf_loop_promoters_restIN", "K562_ctcf_loop_promoters_restOUT", "K562_ctcf_loop_promoters_TF", "K562_ctcf_loop_promoters_TFR", "K562_ctcf_loop_promoters_TI", "K562_ctcf_loop_promoters_TIO", "K562_ctcf_loop_promoters_TO", "K562_ctcf_loop_promoters_TR"]

  # Region Loop
  for region in regionList:

    if("_restIN" in region): regionO = "resty"
    elif("_restOUT" in region): regionO = "restx"
    elif("_TIO" in region): regionO = "both"
    elif("_rest" in region): regionO = "restR"
    elif("_TO" in region): regionO = "left"
    elif("_TI" in region): regionO = "right"
    elif("_TFR" in region): regionO = "bothfr"
    elif("_TF" in region): regionO = "bothf"
    elif("_TR" in region): regionO = "bothr"
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
    geneFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/K562/all_genes.bed"
    ctcfFileName = rl+region+".bed"
    groListFileName = groFile
    signalFileName = dsbFile
    outputFileName = ol+"/"+name+".txt"

    # Execution
    print name
    command = "python 2_ctcfTranscription.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)




# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/1_ctcf_centric_signal/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_DMSO.bam", dl+"TK6_wt_ETO.bam", dl+"K562_DMSO_rep2.3.deep.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  # Region List
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/0_ctcf_motifs/"
  regionList = ["K562_Genes_negativeMotifs", "K562_Genes_plusMotifs", "K562_Promoters_negativeMotifs", "K562_Promoters_plusMotifs", 
"TK6_Genes_negativeMotifs", "TK6_Genes_plusMotifs", "TK6_Promoters_negativeMotifs", "TK6_Promoters_plusMotifs"]

  # Region Loop
  for region in regionList:

    if("_negativeMotifs" in region): regionO = "right"
    else: regionO = "left"

    # Parameters
    cell = dsbFile.split("/")[-1].split("_")[0].upper()
    if("DMSO" in dsbFile): cond = "DMSO"
    else: cond = "ETO"
    if(cell == "TK6"): groFile = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/external/results/counts/GM12878_GROseq.tsv"
    else: groFile = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/external/results/counts/K562_GROseq.tsv"
    name = "_".join([cell, cond, region.split("_")[1], region.split("_")[2]])

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
    command = "python 1_ctcf_centric_signal.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)


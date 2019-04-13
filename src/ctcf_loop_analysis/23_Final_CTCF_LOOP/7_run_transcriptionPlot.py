
# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/7_transcriptionPlot/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_DMSO.bam", dl+"TK6_wt_ETO.bam", dl+"K562_DMSO_rep2.3.deep.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  cell = dsbFile.split("/")[-1].split("_")[0].upper()

  # Region List
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/5_CTCFatAllFeatures/"
  if(cell == "TK6"):
    regionList = ["GM12878_allCtcfInLoopMinusGene", "GM12878_allCtcfInLoopPlusGene"]
    cellD = "GM12878"
  elif(cell == "K562"):
    regionList = ["K562_allCtcfInLoopMinusGene", "K562_allCtcfInLoopPlusGene"]
    cellD = "GM12878"

  # Region Loop
  for region in regionList:

    # regionO
    if("allCtcfInLoopMinusGene" in region): regionO = "rest"
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
    command = "python 7_transcriptionPlot.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)



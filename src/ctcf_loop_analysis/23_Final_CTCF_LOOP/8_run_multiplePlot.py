
# Import
import os
import sys

# DSB List
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/8_multiplePlot/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
dsbFileList = [dl+"TK6_wt_ETO.bam", dl+"K562_Etop20uM_rep2.3.deep.bam"]

# DSB Loop
for dsbFile in dsbFileList:

  # Region List
  cell = dsbFile.split("/")[-1].split("_")[0].upper()
  if(cell == "TK6"): cellD = "GM12878"
  else: cellD = "K562"
  rl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/5_CTCFatAllFeatures/"
  regionList = [cellD+"_allCtcf", cellD+"_allCtcfInLoop", cellD+"_allCtcfOutLoop", cellD+"_allCtcfAllGene", cellD+"_allCtcfPlusGene", cellD+"_allCtcfMinusGene", cellD+"_allCtcfIntergenic"]

  # Region Loop
  for region in regionList:

    if(region == cellD+"_allCtcf"): regionO = "both"
    elif(region == cellD+"_allCtcfInLoop"): regionO = "both"
    elif(region == cellD+"_allCtcfOutLoop"): regionO = "both"
    elif(region == cellD+"_allCtcfAllGene"): regionO = "both"
    elif(region == cellD+"_allCtcfPlusGene"): regionO = "both"
    elif(region == cellD+"_allCtcfMinusGene"): regionO = "rest"
    elif(region == cellD+"_allCtcfIntergenic"): regionO = "both"

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
    command = "python 8_multiplePlot.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
    os.system(command)



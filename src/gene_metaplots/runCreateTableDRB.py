
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Peak list
pn = "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Results/7_facrna_corr/input/CountTablesNewDRB/"
peakList = [pn+"count_gene_drb.txt", pn+"count_gene_wt.txt"]
outLoc = "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Results/7_facrna_corr/tables/DRB_new/count_gene/"

# Peak loop
for featurePeakFileName in peakList:

  # Bam list
  bn = "/media/egg/sbc/AG_Papantonis/eduardo/Roukos_BLISS/data/BLISS2_corrected/mapped/"
  #bamList = [bn+"TK6_DRB_DMSO.bam", bn+"TK6_DRB_ETO.bam", bn+"TK6_wt_DMSO.bam", bn+"TK6_wt_ETO.bam"]
  #bamCountList = [18483127, 18454888, 6715854, 1765410]
  bamList = [bn+"TK6_drb200_DMSO_reampli_LIB1_org7.deep.bam", bn+"TK6_drb200_ETO20_reampli_LIB1_org8.deep.bam", bn+"TK6_untr_DMSO_reampli_LIB1_org3.deep.bam", bn+"TK6_untr_ETO20_reampli_LIB1_org4.deep.bam"]
  bamCountList = [1000000, 1000000, 1000000, 1000000]

  # Bam loop
  for i in range(0,len(bamList)):

    # Folder
    bamFileName = bamList[i]
    bamCount = bamCountList[i]
    outFName = "_".join(featurePeakFileName.split("/")[-1].split(".")[0].split("_")[:-1])

    # Name
    if(featurePeakFileName.split("/")[-1].split(".")[0].split("_")[-1] == "wt"): peakName = "FACTORY_WT"
    else: peakName = "FACTORY_DRB"
    if("_drb200_" in bamFileName): btype = "DRB"
    else:  btype = "WT"
    if("_ETO20_" in bamFileName): bcond = "ETO"
    else:  bcond = "DMSO"
    #btype = bamFileName.split("/")[-1].split(".")[0].split("_")[1].upper()
    #bcond = bamFileName.split("/")[-1].split(".")[0].split("_")[2].upper()
    bamName = "DSB_"+"_".join([btype,bcond])
    name = "__".join([peakName, bamName])
    print name

    # Parameters
    nBins = "15"
    tssExt = "3000"
    bamCount = str(bamCount)
    percentileList = ",".join([str(e) for e in range(100,-1,-5)])
    aliasFileName = "/home/egg/rgtdata/hg19/alias_human.txt"
    genesFileName = "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Results/4_genome_distribution/input/refseq_genes_hg19_filtered.txt"
    featurePeakFileName = featurePeakFileName
    bamFileName = bamFileName
    pol2BamFileName = 
# "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Data/tf/GM12878/sydh/GM12878_POL2.bam"
# "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Data/groseq/GSM12878_GROseq.bw"
# ""
# "/media/egg/sbc/AG_Papantonis/Eduardo/Roukos_Bliss/Data/external/mapped/GM12878_RNAseq.R1.bam"
    tempLocation = outLoc+outFName+"/"+"TEMP/"
    outputFileName = outLoc+outFName+"/"+name+".txt"
    outputFilePlot2Name = outLoc+outFName+"/"+name+"_POL2.txt"



    # Execution
    command = "python createTable.py "+" ".join([nBins, tssExt, bamCount, percentileList, aliasFileName, genesFileName, featurePeakFileName, bamFileName, pol2BamFileName, tempLocation, outputFileName, outputFilePlot2Name])
    os.system(command)



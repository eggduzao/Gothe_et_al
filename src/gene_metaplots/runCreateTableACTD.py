
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Peak list
pn = "/home/egg/Projects/Roukos_Bliss/Results/7_facrna_corr/input/CountTables/"
peakList = [pn+"count_gene_actd.txt", pn+"count_gene_wt.txt", pn+"count_intron_actd.txt", pn+"count_intron_wt.txt", pn+"counts_intron_single_actd.txt", pn+"counts_intron_single_wt.txt"]
outLoc = "/home/egg/Projects/Roukos_Bliss/Results/7_facrna_corr/tables/ACTD/"

# Peak loop
for featurePeakFileName in peakList:

  # Bam list
  bn = "/media/egg/sbc/AG_Papantonis/eduardo/Roukos_BLISS/data/BLISS2_corrected/mapped/"
  bamList = [bn+"TK6_ACTD_DMSO.bam", bn+"TK6_ACTD_ETO.bam", bn+"TK6_wt_DMSO.bam", bn+"TK6_wt_ETO.bam"]
  bamCountList = [12344748, 13353074, 6715854, 1765410]

  # Bam loop
  for i in range(0,len(bamList)):

    # Folder
    bamFileName = bamList[i]
    bamCount = bamCountList[i]
    outFName = "_".join(featurePeakFileName.split("/")[-1].split(".")[0].split("_")[:-1])

    # Name
    if(featurePeakFileName.split("/")[-1].split(".")[0].split("_")[-1] == "wt"): peakName = "FACTORY_WT"
    else: peakName = "FACTORY_ACTD"
    if(len(bamFileName.split("/")[-1].split(".")[0].split("_")) == 4): rep = bamFileName.split("/")[-1].split(".")[0].split("_")[-1].upper()
    else: rep = "MERGED"
    btype = bamFileName.split("/")[-1].split(".")[0].split("_")[1].upper()
    bcond = bamFileName.split("/")[-1].split(".")[0].split("_")[2].upper()
    bamName = "DSB_"+"_".join([rep,btype,bcond])
    name = "::".join([peakName, bamName])
    print name

    # Parameters
    nBins = "15"
    tssExt = "3000"
    bamCount = str(bamCount)
    percentileList = ",".join([str(e) for e in range(95,-1,-5)])
    aliasFileName = "/home/egg/Projects/Roukos_Bliss/Results/4_genome_distribution/input/alias_human.txt"
    genesFileName = "/home/egg/Projects/Roukos_Bliss/Results/4_genome_distribution/input/refseq_genes_hg19_filtered.txt"
    featurePeakFileName = featurePeakFileName
    bamFileName = bamFileName
    tempLocation = outLoc+outFName+"/"+"TEMP/"
    outputFileName = outLoc+outFName+"/"+name+".txt"

    # Execution
    command = "python createTable.py "+" ".join([nBins, tssExt, bamCount, percentileList, aliasFileName, genesFileName, featurePeakFileName, bamFileName, tempLocation, outputFileName])
    os.system(command)



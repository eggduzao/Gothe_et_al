
# Import
import os
import sys

# Folder List
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/0_input/signal/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/1_ctcfSignalPlots/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
gl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/external/results/counts/"
aliasFile = "/home/egg/rgtdata/hg19/alias_human.txt"
geneFile = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/0_genomic_regions/K562/all_genes.bed"
folderList = ["all", "a", "b"]

# Folder Loop
for folder in folderList:

  # Subfolder List
  subfolderList = ["active_genes", "all_genes", "inactive_genes", "intergenic", "active_promoters", "all_promoters", "inactive_promoters"]

  # Subfolder Loop
  for subfolder in subfolderList:

    # Cell List
    cellList = ["K562", "TK6"]

    # Cell Loop
    for cell in cellList:

      # Files
      dsbFile = dl+"TK6_wt_ETO.bam"
      groFile = gl+"GM12878_GROseq.tsv"
      if(cell == "K562"):
        dsbFile = dl+"K562_Etop20uM_rep2.3.deep.bam"
        groFile = gl+"K562_GROseq.tsv"

      # Condition List
      if(folder == "b" and subfolder == "intergenic"):
        condList = ["intergenic_all", "intergenic_inside", "intergenic_outside"]
        roList = ["restR", "resty", "restx"]
      elif(folder == "b"):
        condList = ["TIO", "TO", "TI", "TF", "TR", "TFR"]
        if(subfolder == "inactive_genes" or subfolder == "inactive_promoters"): roList = ["rest", "leftI", "rightI", "bothfI", "bothrI", "bothfrI"]
        else: roList = ["both", "left", "right", "bothf", "bothr", "bothfr"]
      elif(subfolder == "intergenic"):
        condList = ["intergenic_all", "intergenic_forward", "intergenic_reverse"]
        roList = ["rest", "rest", "rest"]
      else:
        condList = ["MFR_TFR", "MFR_TF", "MFR_TR", "MF_TFR", "MR_TFR", "MFTF_MRTR", "MFTR_MRTF"]
        if(subfolder == "inactive_genes" or subfolder == "inactive_promoters"): roList = ["rest", "leftI", "rightI", "rest", "rest", "rest", "rest"]
        else: roList = ["both", "left", "right", "both", "both", "both", "both"]

      # Condition Loop
      for i in range(0,len(condList)):

        cond = condList[i]
        ro = roList[i]
        inLoc = il+folder+"/"+subfolder+"/"
        outLoc = ol+folder+"/"+subfolder+"/"
        os.system("mkdir -p "+outLoc)

        # Input
        higherSide = ro
        ctcfExt = "500"
        percentileList = ",".join([str(e) for e in range(100,-1,-5)])
        aliasFileName = aliasFile
        geneFileName = geneFile
        ctcfFileName = inLoc+cell+"_"+cond+".bed"
        groListFileName = groFile
        signalFileName = dsbFile
        outputFileName = outLoc+cell+"_"+cond+".txt"

        # Execution
        print folder, subfolder, cell, cond
        command = "time python 1_ctcfSignal.py "+" ".join([higherSide, ctcfExt, percentileList, aliasFileName, geneFileName, ctcfFileName, groListFileName, signalFileName, outputFileName])
        os.system(command)




# Import
import os
import sys

# Cell List
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/15_split_ctcf_categories/1_ctcf_gene_strand_categories/HINT/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/0_ctcf_motifs/"
cellList = ["GM12878", "K562"]

# Cell Loop
for cell in cellList:

  if(cell == "GM12878"): cellN = "TK6"
  else: cellN = cell

  # Dataset List
  datasetList = [["negativeMotifsInsideNegativeGenes", "negativeMotifsInsidePlusGenes"],
                 ["plusMotifsInsideNegativeGenes", "plusMotifsInsidePlusGenes"],
                 ["negativeMotifsInsidePlusPromoters", "negativeMotifsInsideNegativePromoters"],
                 ["plusMotifsInsidePlusPromoters", "plusMotifsInsideNegativePromoters"]]

  # Dataset Loop
  for dataset in datasetList:

    # Name
    motif = dataset[0].split("Inside")[0]
    feature = dataset[0].split("Negative")[-1].split("Plus")[-1]
    name = "_".join([cellN, feature, motif])

    # Input 
    inFileName1 = il+cell+"/"+dataset[0]+".bed"
    inFileName2 = il+cell+"/"+dataset[1]+".bed"
    tempLoc = ol+"TEMP/"
    outputFileName = ol+name+".bed"

    # Execution
    print name
    command = "python 0_ctcf_motifs.py "+" ".join([inFileName1, inFileName2, tempLoc, outputFileName])
    os.system(command)



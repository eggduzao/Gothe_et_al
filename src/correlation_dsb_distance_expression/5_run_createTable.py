
# Import
import os
import sys

mode = "regular" # "regular" "raw"

# Cell List
aliasFile = "/Users/egg/rgtdata/hg19/alias_human_booster.txt"
mllFile = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/list_of_genes.txt"
allGenesFile = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/all_genes.bed"
ml = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/0_Input/mll_helpers/"
dl = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/1_Dsb/"
gl = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/2_Exp/"
ll = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/4_Extended_Anchors/"
ol = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/5_Tables/"
cellList = ["K562", "TK6", "CD34"]

# Cell Loop
for cell in cellList:

  targetTableList = []
  outputHelperList = []
  outputTableList = []
  outputTableSingle = []

  # Count Type List
  ctList = ["promoter_gene", "promoter", "gene"]

  # Count Type Loop
  for ct in ctList:

    # Type of Anchor
    tlList = ["anchors_with_ctcf", "anchors_with_and_wo_ctcf"]

    # Type of Anchor Loop
    for tl in tlList:

      # Maximum Distance List (from the loop anchors)
      maxDistList = ["201", "501", "1001"]

      # Maximum Distance Loop
      for md in maxDistList:

        # Parameters
        outLoc = ol + ct + "/" + tl + "/" + md + "/"

        # Input
        maxDist = md
        countType = ct
        aliasFileName = aliasFile
        mllGenesFileName = mllFile
        allGenesFileName = allGenesFile
        expressionFileName = gl + cell + "_expression_filtered.txt"
        dsbFileName = dl + cell + "_ETO.bam"
        distFileName = ll + cell + "_" + tl + ".txt"
        outputFileName = outLoc + cell + "_table.txt"
  
        if(mode == "regular"):

          # Execution
          print mode, cell, ct, tl, md
          command = "python 5_createTable.py "+" ".join([maxDist, countType, aliasFileName, mllGenesFileName, allGenesFileName, expressionFileName, dsbFileName, distFileName, outputFileName])
          os.system(command)

          # Input
          inputTableFileName = outputFileName
          outputTableFileName = outLoc + cell + "_table_filt.txt"

          command = "python 5_filterGenes.py "+" ".join([inputTableFileName, outputTableFileName])
          os.system(command)
          os.system("mv "+outputTableFileName+" "+inputTableFileName)

          targetTableList.append(inputTableFileName)
          outputHelperList.append(outLoc + cell + "_helper.txt")
          outputTableList.append(outLoc + cell + "_table_filt.txt")
          outputTableSingle.append(outLoc + cell + "_table_single.txt")

        if(mode == "raw"):
         
          # Execution
          print mode, cell, ct, tl, md
          outputFileName = outLoc + cell + "_table_raw.txt"
          command = "python 5_createTableRaw.py "+" ".join([maxDist, countType, aliasFileName, mllGenesFileName, allGenesFileName, expressionFileName, dsbFileName, anchorAllFilePrefix, outputFileName])
          os.system(command)

  # Input
  aliasFileName = aliasFile
  originalMllListFileName = mllFile
  inputMllFileName = ml + cell + ".txt"
  targetTableFileNameList = ",".join(targetTableList)
  outputHelperFileNameList = ",".join(outputHelperList)
  outputTableFileNameList = ",".join(outputTableList)
  outputTableSingleFileNameList = ",".join(outputTableSingle)

  # Execution
  print "HELPER: ", cell
  command = "python 5_putHelperMll.py "+" ".join([aliasFileName, originalMllListFileName, inputMllFileName, targetTableFileNameList, outputHelperFileNameList, outputTableFileNameList, outputTableSingleFileNameList])
  os.system(command)


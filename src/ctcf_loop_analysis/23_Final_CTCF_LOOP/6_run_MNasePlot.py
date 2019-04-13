
# Import
import os
import sys

cml = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/5_CTCFatAllFeatures/"
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
cl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/tf/"
ml = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/mnase/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/"
cellList = ["GM12878", "K562"]

# Cell Loop
for cell in cellList:

  if(cell == "GM12878"):
    etoName = "TK6_wt_ETO.bam"
    mnaseName = "wgEncodeSydhNsomeGm12878Sig.bigWig"
  elif(cell == "K562"):
    etoName = "K562_Etop20uM_rep2.3.deep.bam"
    mnaseName = "wgEncodeSydhNsomeK562Sig.bigWig"

  # Input
  higherSide = "both"
  ctcfExt = "500"
  ctcfFileName = cml+cell+"_allCtcfInLoop.bed"
  etoBamFileName = dl+etoName
  ctcfBamFileName = cl+cell+"/sydh/"+cell+"_CTCF.bam"
  mnaseBwFileName = ml+mnaseName
  outputEtoFileName = ol+cell+"_ETO.txt"
  outputCtcfFileName = ol+cell+"_CTCF.txt"
  outputMNaseFileName = ol+cell+"_MNASE.txt"

  # Execution
  print cell
  command = "python 6_MNasePlot.py "+" ".join([higherSide, ctcfExt, ctcfFileName, etoBamFileName, ctcfBamFileName, mnaseBwFileName, outputEtoFileName, outputCtcfFileName, outputMNaseFileName])
  os.system(command)



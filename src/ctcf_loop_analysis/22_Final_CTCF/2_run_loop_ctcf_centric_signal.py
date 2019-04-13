
# Import
import os
import sys

# Signal List
dl = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/BLISS2_corrected/mapped/"
ll = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Data/loops/"
ol = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/2_ctcf_at_loops/"
signalList = ["TK6_wt_DMSO.bam", "TK6_wt_ETO.bam", "K562_DMSO_rep2.3.deep.bam", "K562_Etop20uM_rep2.3.deep.bam"]

# Signal Loop
for signal in signalList:

  # Parameters
  cell = signal.split("_")[0]
  cond = "ETO"
  if("DMSO" in signal): cond = "DMSO"
  name = "_".join([cell, cond])
  print name

  # Input
  ctcfExt = "500"
  loopFileName = ll+cell+"_CTCF_loops.txt"
  signalFileName = dl+signal
  outputFileName = ol+name+".txt"

  # Execution
  command = "python 2_loop_ctcf_centric_signal.py "+" ".join([ctcfExt, loopFileName, signalFileName, outputFileName])
  os.system(command)



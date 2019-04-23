
###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import os
import sys
from optparse import SUPPRESS_HELP
import warnings
warnings.filterwarnings("ignore")

# Internal
from src import __version__
from ..Util import PassThroughOptionParser
from createTable import create_table

"""
gene_metaplots

This program creates a meta-plot of a gene based on the signals and regions given.

Dependencies:
- Numpy >= 1.13.1
- Pysam >= 0.11.2.2
- pyBigWig

Authors: Eduardo G. Gusmao.
"""

###################################################################################################
# Functions
###################################################################################################

# Uncompressing files
def uncompressing_files(compressed_file_name, uncompressed_file_name):

  # Taking suffix
  cc = compressed_file_name.split(".")

  # Uncompressing
  if(cc[-1] == "gz" or cc[-1] == "tar" or cc[-2] == "tar" or cc[-1] == "zip"):
    if(cc[-1] == "tar"):
      command = "tar -O "+compressed_file_name+" > "+uncompressed_file_name
      os.system(command)
    if(cc[-2] == "tar"):
      if(cc[-1] == "gz"):
        command = "tar -xO "+compressed_file_name+" > "+uncompressed_file_name
        os.system(command)
      else: print("ERROR: Unrecognized tarball.")
    elif(cc[-1] == "gz"):
      command = "gzip -cd "+dsb_file_name+" > "+uncompressed_file_name
      os.system(command)
    elif(cc[-1] == "zip"):
      command = "unzip -p "+dsb_file_name+" > "+uncompressed_file_name
      os.system(command)
    else: print("ERROR: We only support tar.gz, .gz and .zip compressions.")
  else: uncompressed_file_name = compressed_file_name

###################################################################################################
# Main
###################################################################################################

def main():
  """
  Main function that creates a gene meta-plot and plot it using R scripts.

  Keyword arguments: None

  Return: None
  """

  ###################################################################################################
  # Processing Input Arguments
  ###################################################################################################

  # Parameters
  usage_message = ("\n--------------------------------------------------\n"
                   "This program creates a meta-plot of a gene based on \n"
                   "the signals and regions given.\n\n"

                   "The program should be called as:\n"
                   "%prog <args>\n\n"

                   "For more information on the arguments please type:\n"
                   "%prog --help\n\n"

                   "For more information, please refer to the original paper:\n"
                   "Placeholder.\n\n"

                   "For further questions or comments please contact:\n"
                   "eduardogade@gmail.com\n"
                   "--------------------------------------------------")
  version_message = "Gothe et al. - Spatial chromosome folding and active transcription drive DNA fragility and formation of oncogenic MLL translocations. Version: " + str(__version__)

  # Initializing Option Parser
  parser = PassThroughOptionParser(usage=usage_message, version=version_message)

  """ Sample input options:

  # string (file, flags, etc):
  parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19", help=("Placeholder."))

  # int
  parser.add_option("--dnase-sg-window-size", dest="dnase_sg_window_size", type="int", metavar="INT", default=9, help=("Placeholder."))

  # float (hidden from help menu)
  parser.add_option("--dnase-norm-per", dest="dnase_norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)

  # boolean
  parser.add_option("--estimate-bias-correction", dest="estimate_bias_correction", action="store_true", default=False, help=("Placeholder."))
  """

  # Input Options
  parser.add_option("--nBins", dest="nBins", type="int", metavar="INT", default=15, help=("Number of bins in which the meta-plot region (6Kbp) is going to be divided to average the signal. The original was set to 15"))
  parser.add_option("--tssExt", dest="tssExt", type="int", metavar="INT", default=3000, help=("The size, in bp, of promoter regions. The original was set to 2Kbp."))
  parser.add_option("--bamCount", dest="bamCount", type="int", metavar="INT", default=1000000, help=("The total number of reads in the BAM file containing the signal to plot."))
  parser.add_option("--aliasFileName", dest="aliasFileName", type="string", metavar="FILE", default=None, help=("File containing gene aliases."))
  parser.add_option("--genesFileName", dest="genesFileName", type="string", metavar="FILE", default=None, help=("A file containing the location of genes. In this particular case the format has to be ENCODE's refseq table."))
  parser.add_option("--expressionList", dest="expListFileName", type="string", metavar="FILE", default=None, help=("A plain text (tab-separated) file containing the genes in the first column and their expression in the second column."))
  parser.add_option("--bamFileName", dest="bamFileName", type="string", metavar="FILE", default=None, help=("A BAM file containing the signal in which the meta-plot will be calculated."))
  parser.add_option("--temp", dest="tempLocation", type="string", metavar="PATH", default=None, help=("Temporary location to aid in the execution."))
  parser.add_option("--outputFileName", dest="outputFileName", type="string", metavar="FILE", default=None, help=("Output file name."))

  # Processing Options
  options, arguments = parser.parse_args() 

  # General options
  nBins = options.nBins
  tssExt = options.tssExt
  bamCount = options.bamCount
  percentileList = range(100,-1,-1)
  aliasFileName = options.aliasFileName
  genesFileName = options.genesFileName
  expListFileName = options.expListFileName
  bamFileName = options.bamFileName
  tempLocation = options.tempLocation
  outputFileName = options.outputFileName

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not nBins): print(argument_error_message)
  if(not tssExt): print(argument_error_message)
  if(not bamCount): print(argument_error_message)
  if(not percentileList): print(argument_error_message)
  if(not aliasFileName): print(argument_error_message)
  if(not genesFileName): print(argument_error_message)
  if(not expListFileName): print(argument_error_message)
  if(not bamFileName): print(argument_error_message)
  if(not tempLocation): print(argument_error_message)
  if(not outputFileName): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Uncompress aliasFileName
  aliasFileNameUnc = tempLocation + "aliasFileNameUnc.bed"
  uncompressing_files(aliasFileName, aliasFileNameUnc)

  # Uncompress genesFileName
  genesFileNameUnc = tempLocation + "genesFileNameUnc.bed"
  uncompressing_files(genesFileName, genesFileNameUnc)

  # Uncompress expListFileName
  expListFileNameUnc = tempLocation + "expListFileNameUnc.txt"
  uncompressing_files(expListFileName, expListFileNameUnc)

  # Uncompress bamFileName
  bamFileNameUnc = tempLocation + "bamFileNameUnc.bam"
  uncompressing_files(bamFileName, bamFileNameUnc)

  # Creating table
  create_table(nBins, tssExt, bamCount, percentileList, aliasFileNameUnc, genesFileNameUnc, expListFileNameUnc, bamFileNameUnc, tempLocation, outputFileName)

  # Script path
  script_path = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"

  # Creating plot
  plotFileName = ".".join(outputFileName.split(".")[:-1]) + ".pdf"
  command = "Rscript "+script_path+"correlation.R "+" ".join([str(nBins), outputFileName, plotFileName])
  os.system(command)


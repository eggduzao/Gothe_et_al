
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
from optparse import SUPPRESS_HELP
import warnings
warnings.filterwarnings("ignore")

# Internal
from src import __version__
from ..Util import PassThroughOptionParser
from createTable import create_table

# External
# TODO

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
  if(cc[-2] == "tar"):
    if(cc[-1] == "gz"):
      command = "tar -xO "+compressed_file_name+" > "+uncompressed_file_name
      os.system(command)
    else:
      command = "tar -O "+compressed_file_name+" > "+uncompressed_file_name
      os.system(command)
  elif(cc[-1] == "gz"):
    command = "gzip -cd "+dsb_file_name+" > "+uncompressed_file_name
    os.system(command)
  elif(cc[-1] == "zip"):
    command = "unzip -p "+dsb_file_name+" > "+uncompressed_file_name
    os.system(command)

  else: print("ERROR: We only support tar.gz, .gz and .zip compressions.")


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
  parser.add_option("--nBins", dest="nBins", type="int", metavar="INT", default=15, help=("Placeholder."))
  parser.add_option("--tssExt", dest="tssExt", type="int", metavar="INT", default=3000, help=("Placeholder."))
  parser.add_option("--bamCount", dest="bamCount", type="int", metavar="INT", default=1000000, help=("Placeholder."))
  parser.add_option("--percentileList", dest="percentileList", type="string", metavar="LIST", default=",".join([str(e) for e in range(100,-1,-20)]), help=("Placeholder."))
  parser.add_option("--aliasFileName", dest="aliasFileName", type="string", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--genesFileName", dest="genesFileName", type="string", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--featurePeakFileName", dest="featurePeakFileName", type="FILE", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--bamFileName", dest="bamFileName", type="string", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--tempLocation", dest="tempLocation", type="string", metavar="PATH", default=None, help=("Placeholder."))
  parser.add_option("--outputFileName", dest="outputFileName", type="string", metavar="FILE", default=None, help=("Placeholder."))

  # Processing Options
  options, arguments = parser.parse_args() 

  # General options
  nBins = options.nBins
  tssExt = options.tssExt
  bamCount = options.bamCount
  percentileList = options.percentileList
  aliasFileName = options.aliasFileName
  genesFileName = options.genesFileName
  featurePeakFileName = options.featurePeakFileName
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
  if(not output_file_name): print(argument_error_message)
  if(not bamFileName): print(argument_error_message)
  if(not tempLocation): print(argument_error_message)
  if(not outputFileName): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Uncompress aliasFileName
  aliasFileNameUnc = aliasFileName
  uncompressing_files(aliasFileName, aliasFileNameUnc)

  # Uncompress genesFileName
  genesFileNameUnc = genesFileName
  uncompressing_files(genesFileName, genesFileNameUnc)

  # Uncompress featurePeakFileName
  featurePeakFileNameUnc = featurePeakFileName
  uncompressing_files(featurePeakFileName, featurePeakFileNameUnc)

  # Uncompress bamFileName
  bamFileNameUnc = bamFileName
  uncompressing_files(bamFileName, bamFileNameUnc)

  # Creating table
  create_table(nBins, tssExt, bamCount, percentileList, aliasFileNameUnc, genesFileNameUnc, featurePeakFileNameUnc, bamFileNameUnc, tempLocation, outputFileName)

  # Creating plot
  plotFileName = ".".join(outputFileName.split(".")[:-1]) + ".pdf"


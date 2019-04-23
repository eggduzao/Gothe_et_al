
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
corr-dsb-feat

This program creates scatterplots and lineplots with the correlation of user-defined
features and DSBs.

Dependencies:
- numpy
- pysam
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
  Main function that creates scatterplots and lineplots with the correlation of user-defined
  features and DSBs.

  Keyword arguments: None

  Return: None
  """

  ###################################################################################################
  # Processing Input Arguments
  ###################################################################################################

  # Parameters
  usage_message = ("\n--------------------------------------------------\n"
                   "This program creates scatterplots and lineplots with the\n"
                   "correlation of user-defined features and DSBs.\n\n"

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
  parser.add_option("--half-ext", dest="half_ext", type="int", metavar="INT", default=500, help=("Half the distance (in bp) from the middle of the feature to calculate the correlation."))
  parser.add_option("--regions", dest="feature_summit_file_name", type="string", metavar="FILE", default=None, help=("A file containing the regions in which the features will be centered at to calculate the correlation."))
  parser.add_option("--signal-label-list", dest="bam_names", type="string", metavar="NAME_1[,NAME_2,...,NAME_N]", default=None, help=("A comma-separated list of labels for each signal (features) to be plot."))
  parser.add_option("--signal-count-list", dest="bam_counts", type="string", metavar="INT_1[,INT_2,...,INT_N]", default=None, help=("A comma-separated list containing the total read count of each signal's (features's) BAM file."))
  parser.add_option("--signal-file-list", dest="bam_list", type="string", metavar="FILE_1[,FILE_2,...,FILE_N]", default=None, help=("A comma-separated list of BAM files for each signal (features) to be plot."))
  parser.add_option("--temp", dest="temp_location", type="string", metavar="PATH", default=None, help=("Temporary location to aid in the execution."))
  parser.add_option("--output-file", dest="output_file_name", type="string", metavar="FILE", default=None, help=("Output file name."))

  # Processing Options
  options, arguments = parser.parse_args()

  # General options
  half_ext = options.half_ext
  feature_summit_file_name = options.feature_summit_file_name
  bam_names = options.bam_names
  bam_counts = options.bam_counts
  bam_list = options.bam_list
  temp_location = options.temp_location
  output_file_name = options.output_file_name

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not half_ext): print(argument_error_message)
  if(not feature_summit_file_name): print(argument_error_message)
  if(not bam_names): print(argument_error_message)
  if(not bam_counts): print(argument_error_message)
  if(not bam_list): print(argument_error_message)
  if(not temp_location): print(argument_error_message)
  if(not output_file_name): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Uncompress feature_summit_file_name
  feature_summit_file_name_unc = temp_location + "feature_summit_file_name_unc.bed"
  uncompressing_files(feature_summit_file_name, feature_summit_file_name_unc)

  # Uncompress bam_list
  bam_list_unc = []
  counter = 1
  for bam_file_name in bam_list:
    bam_file_name_unc = temp_location + "bam_file_name_unc" + str(counter) + ".bam"
    uncompressing_files(bam_file_name, bam_file_name_unc)
    bam_list_unc.append(bam_file_name_unc)
    counter += 1

  # Create table
  create_table(half_ext, feature_summit_file_name_unc, bam_names, bam_counts, bam_list_unc, output_file_name)

  # Script path
  script_path = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"

  # Creating plot
  graphWidth = "10"
  marginX = "9.5"
  inputTableFileName = output_file_name
  outputFileName = ".".join(output_file_name.split(".")[:-1]) + ".pdf"
  outputLocation = "/".join(output_file_name.split("/")[:-1]) + "/scatterplots/"
  command = "mkdir -p "+outputLocation
  os.system(command)
  command = "Rscript "+script_path+"correlation.R "+" ".join([graphWidth, marginX, inputTableFileName, outputFileName, outputLocation])
  os.system(command)


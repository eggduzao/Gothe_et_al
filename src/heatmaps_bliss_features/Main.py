
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
from createHeatmap import create_heatmap

"""
heatmap-dsb

This program creates a centered around the features given for the signal also given.
Furthermore, it sorts the heatmap in decreasing order by intensity of a given list (e.g. DSBs).

Dependencies:
- deeptools

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


###################################################################################################
# Main
###################################################################################################

def main():
  """
  Main function that creates a centered around the features given for the signal also given.
  Furthermore, it sorts the heatmap in decreasing order by intensity of a given list (e.g. DSBs)

  Keyword arguments: None

  Return: None
  """

  ###################################################################################################
  # Processing Input Arguments
  ###################################################################################################

  # Parameters
  usage_message = ("\n--------------------------------------------------\n"
                   "This program calculates:\n"
                   "This program creates a centered around the features given for\n"
                   "the signal also given. Furthermore, it sorts the heatmap in\n"
                   "decreasing order by intensity of a given list (e.g. DSBs).\n\n"

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
  parser.add_option("--half-ext", dest="half_ext", type="int", metavar="INT", default=200, help=("Placeholder."))
  parser.add_option("--regions", dest="feature_summit_file_name", type="string", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--signal-file", dest="signal_file_name", type="string", metavar="FILE", default=None, help=("Placeholder."))
  parser.add_option("--signal-label", dest="signal_label", type="string", metavar="STRING", default=None, help=("Placeholder."))
  parser.add_option("--temp", dest="temp_location", type="string", metavar="PATH", default=None, help=("Placeholder."))
  parser.add_option("--output-file", dest="output_file_name", type="string", metavar="FILE", default=None, help=("Placeholder."))

  # Processing Options
  options, arguments = parser.parse_args()

  # General options
  half_ext = options.half_ext
  feature_summit_file_name = options.feature_summit_file_name
  signal_file_name = options.signal_file_name
  signal_label = options.signal_label
  temp_location = options.temp_location
  output_file_name = options.output_file_name

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not half_ext): print(argument_error_message)
  if(not feature_summit_file_name): print(argument_error_message)
  if(not signal_file_name): print(argument_error_message)
  if(not signal_label): print(argument_error_message)
  if(not temp_location): print(argument_error_message)
  if(not output_file_name): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Uncompress feature_summit_file_name
  feature_summit_file_name_unc = feature_summit_file_name
  uncompressing_files(feature_summit_file_name, feature_summit_file_name_unc)

  # Uncompress signal_file_name
  signal_file_name_unc = signal_file_name
  uncompressing_files(signal_file_name, signal_file_name_unc)

  # Create heatmap
  create_heatmap(half_ext, feature_summit_file_name_unc, signal_file_name_unc, signal_label, temp_location, output_file_name)


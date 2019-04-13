
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
# TODO
from optparse import SUPPRESS_HELP
import warnings
warnings.filterwarnings("ignore")

# Internal
from src import __version__
from ..Util import PassThroughOptionParser
# TODO

# External
# TODO

"""
corr-dsb-dist-exp

This program calculates:
- Distances from a list of genes to their closest anchors."
- The expression of these genes."
- The proportion of double-strand breaks (DSBs) in these genes."

Dependencies:
- 

Authors: Eduardo G. Gusmao.
"""

###################################################################################################
# Functions
###################################################################################################

# TODO

###################################################################################################
# Main
###################################################################################################

def main():
  """
  Main function that creates a triple correlation table and plots them using R scripts.

  Keyword arguments: None

  Return: None
  """

  ###################################################################################################
  # Processing Input Arguments
  ###################################################################################################

  # Initializing ErrorHandler
  error_handler = ErrorHandler()

  # Parameters
  usage_message = ("\n--------------------------------------------------\n"
                   "This program calculates:\n"
                   "- Distances from a list of genes to their closest anchors.\n"
                   "- The expression of these genes.\n\n"
                   "- The proportion of double-strand breaks (DSBs) in these genes.\n\n"

                   "The program should be called as:\n"
                   "%prog <args>\n\n"

                   "For more information on the arguments please type:\n"
                   "%prog --help\n\n"

                   "For more information, please refer to the original paper in :\n"
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
  parser.add_option("--max-dist", dest="max_dist", type="int", metavar="INT", default=200, help=("Placeholder.")) # 1
  parser.add_option("--alias-file", dest="alias_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 2
  parser.add_option("--genes-file", dest="genes_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 3
  parser.add_option("--expression-file", dest="exp_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 4
  parser.add_option("--dsb-file", dest="dsb_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 5
  parser.add_option("--distance-file", dest="dist_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 6
  parser.add_option("--output-file", dest="output_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 7

  # Processing Options
  options, arguments = parser.parse_args()
  if len(arguments) < 7:
    print(usage_message)
    exit(1)

  # General options
  max_dist = options.max_dist
  alias_file_name = options.alias_file_name
  genes_file_name = options.genes_file_name
  exp_file_name = options.exp_file_name
  dsb_file_name = options.dsb_file_name
  dist_file_name = options.dist_file_name
  output_file_name = options.output_file_name

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not max_dist): print(argument_error_message)
  if(not alias_file_name): print(argument_error_message)
  if(not genes_file_name): print(argument_error_message)
  if(not exp_file_name): print(argument_error_message)
  if(not dsb_file_name): print(argument_error_message)
  if(not dist_file_name): print(argument_error_message)
  if(not output_file_name): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # TODO




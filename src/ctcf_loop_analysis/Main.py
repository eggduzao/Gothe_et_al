
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
from ctcfSignal import ctcf_signal
from ..correlation_dsb_distance_expression.processDsbFile import create_bam_file
from ..correlation_dsb_distance_expression.processExpFile import create_exp_file
from ..correlation_dsb_distance_expression.processHicFile import create_hic_file

"""
ctcf-loop

This program replicates the CTCF-DSB plots present in the original paper.

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
  Main function that replicates the CTCF-DSB plots present in the original paper.

  Keyword arguments: None

  Return: None
  """

  ###################################################################################################
  # Processing Input Arguments
  ###################################################################################################

  # Parameters
  usage_message = ("\n--------------------------------------------------\n"
                   "This program replicates the CTCF-DSB plots present\n"
                   " in the original paper.\n\n"

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
  parser.add_option("--region-type", dest="region_type", type="string", metavar="STRING", default=0, help=("The type of CTCF file generated. For all the options check our 'full documentation' page. This option is important to ensure a correct and fast analysis."))
  parser.add_option("--alias-file", dest="alias_file_name", type="string", metavar="FILE", default=None, help=("File containing gene aliases."))
  parser.add_option("--gene-file", dest="gene_file_name", type="string", metavar="FILE", default=None, help=("A simple BED file containing the location of genes."))
  parser.add_option("--ctcf-file", dest="ctcf_file_name", type="string", metavar="FILE", default=None, help=("A file containing the particular genes (or other elements) that overlapped a CTCF factor."))
  parser.add_option("--expression-file", dest="expression_file_name", type="string", metavar="FILE", default=None, help=("A plain text (tab-separated) file containing the genes in the first column and their expression in the second column."))
  parser.add_option("--dsb-file", dest="dsb_file_name", type="string", metavar="FILE", default=None, help=("A BAM file containing all the DSBs."))
  parser.add_option("--temp", dest="temp_location", type="string", metavar="PATH", default=None, help=("Temporary location to aid in the execution."))
  parser.add_option("--output-file", dest="output_file_name", type="string", metavar="FILE", default=None, help=("Output file name."))

  # Processing Options
  options, arguments = parser.parse_args()

  # General options
  region_type = options.region_type
  ctcf_resolution = 500
  percentile_list = range(100,-1,-1)
  alias_file_name = options.alias_file_name
  gene_file_name = options.gene_file_name
  ctcf_file_name = options.ctcf_file_name
  expression_file_name = options.expression_file_name
  dsb_file_name = options.dsb_file_name
  temp_location = options.temp_location
  output_file_name = options.output_file_name

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not region_type): print(argument_error_message)
  if(not ctcf_resolution): print(argument_error_message)
  if(not percentile_list): print(argument_error_message)
  if(not alias_file_name): print(argument_error_message)
  if(not gene_file_name): print(argument_error_message)
  if(not ctcf_file_name): print(argument_error_message)
  if(not expression_file_name): print(argument_error_message)
  if(not dsb_file_name): print(argument_error_message)
  if(not temp_location): print(argument_error_message)
  if(not output_file_name): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Uncompress alias_file_name
  alias_file_name_unc = temp_location + "alias_file_name_unc.bed"
  uncompressing_files(alias_file_name, alias_file_name_unc)

  # Uncompress gene_file_name
  gene_file_name_unc = temp_location + "gene_file_name_unc.bed"
  uncompressing_files(gene_file_name, gene_file_name_unc)

  # Uncompress ctcf_file_name
  ctcf_file_name_unc = temp_location + "ctcf_file_name_unc.bed"
  uncompressing_files(ctcf_file_name, ctcf_file_name_unc)

  # Uncompress expression_file_name
  expression_file_name_unc = temp_location + "expression_file_name_unc.txt"
  uncompressing_files(expression_file_name, expression_file_name_unc)

  # Uncompress dsb_file_name
  dsb_file_name_unc = temp_location + "dsb_file_name_unc.bam"
  uncompressing_files(dsb_file_name, dsb_file_name_unc)

  # Create ctcf table
  ctcf_signal(region_type, ctcf_resolution, percentile_list, alias_file_name_unc, gene_file_name_unc, ctcf_file_name_unc, expression_file_name_unc, dsb_file_name_unc, output_file_name)

  # Script path
  script_path = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"

  # Creating plot
  inputTableFileName = output_file_name
  outputFileNameAggr = ".".join(output_file_name.split(".")[:-1]) + "_aggr.pdf"
  outputFileNameHeat = ".".join(output_file_name.split(".")[:-1]) + "_heat.pdf"
  outputFileNameCorr = ".".join(output_file_name.split(".")[:-1]) + "_corr.pdf"
  command = "Rscript "+script_path+"lineplot_correlation_heatmap.R "+" ".join([region_type, str(ctcf_resolution), inputTableFileName, outputFileNameAggr, outputFileNameHeat, outputFileNameCorr])
  os.system(command)


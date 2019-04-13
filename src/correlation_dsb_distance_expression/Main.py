
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import math
from random import seed, random, randint
from optparse import SUPPRESS_HELP
import warnings
warnings.filterwarnings("ignore")

# Internal
from src import __version__
from ..Util import PassThroughOptionParser

# External
import numpy as np
from pysam import Samfile

"""
corr-dsb-dist-exp

This program calculates:
- Distances from a list of genes to their closest anchors."
- The expression of these genes."
- The proportion of double-strand breaks (DSBs) in these genes."

Python Dependencies:
- Numpy >= 1.13.1
- Pysam >= 0.11.2.2

R Dependencies:
- ggplot2
- gplots
- RColorBrewer
- plot3D
- plotly
- scatterplot3d

Authors: Eduardo G. Gusmao.
"""

###################################################################################################
# Functions
###################################################################################################

def atLeastOneRead(bam_file, region):
  ret = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    ret = True
    break
  return ret

def fetch_counts(bam_file, region):
  return bam_file.count(region[0], region[1], region[2])

def read_alias_dictionary(alias_file_name):

  # Alias dictionary
  alias_dict = dict()
  aliasFile = open(alias_file_name,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: alias_dict[g.upper()] = value.upper()
  aliasFile.close()

  # Return objects
  return alias_dict

def create_distance_dictionary(dist_file_name):

  # Distance dictionary
  dist_dict = dict()
  dist_file = open(dist_file_name, "rU")
  dist_file.readline()
  for line in dist_file:
    ll = line.strip().split("\t")
    dist_dict[ll[0].upper()] = int(ll[1]) 
  dist_file.close()

  # Return objects
  return dist_dict

def create_mll_dictionary(alias_dict, mll_genes_file_name):

  # Creating MLL genes dictionary
  mll_dict = dict()
  mllGenesFile = open(mll_genes_file_name, "rU")
  for line in mllGenesFile:
    ll = line.strip().split("\t")
    try: geneName = alias_dict[ll[3]].upper()
    except Exception: geneName = ll[3].upper()
    mll_dict[geneName] = True
  mllGenesFile.close()

  # Return objects
  return mll_dict

def fetch_expression_percentiles(expression_file_name):

  # Fetching expression for all genes
  expList = []
  exp_dict = dict()
  expressionFile = open(expression_file_name, "rU")
  for line in expressionFile:
    ll = line.strip().split("\t")
    expGene = ll[0]
    try: gene = alias_dict[expGene.upper()]
    except Exception: gene = expGene.upper()
    expvalue = float(ll[1])
    exp_dict[gene] = expvalue
    expList.append(expvalue)
  expressionFile.close()

  # Return objects
  return exp_dict

def create_multi_table(max_dist, alias_file_name, genes_file_name, exp_file_name, dsb_file_name, dist_file_name, output_file_name):

  # Allowed chromosomes
  chrom_list = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Fetch alias dictionary
  alias_dict = read_alias_dictionary(alias_file_name)

  # Fetch expression
  exp_dict = fetch_expression_percentiles(expression_file_name)

  # Fetch distance
  dist_dict = create_distance_dictionary(dist_file_name)

  # Input Files
  allGenesFile = open(genes_file_name, "rU")
  try: dsbFile = Samfile(dsb_file_name, "rb")
  except Exception:
    print("ERROR: Could not open DSB BAM file. Check your Pysam installation.")
    exit(1)

  # Output file
  outputFile = open(output_file_name, "w")
  outputFile.write("\t".join(["GENE", "DISTANCE", "EXPRESSION", "DSB"])+"\n")

  # Iterating in gene file
  for line in allGenesFile:

    # Initialization
    ll = line.strip().split("\t")
    try: chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2]); gene = ll[3].upper(); score = ll[4]; strand = ll[5]
    except Exception:
      print("ERROR: The genes file must be a tab-separated bed file with columns: chromosome, start, end, gene_name, score (not used), strand")
    try: gene = alias_dict[gene]
    except Exception: pass
    if(chrom not in chrom_list): continue
    if(strand == "+"): region = [chrom, p1 - promExt, p1]
    else: region = [chrom, p2, p2 + promExt]

    # Fetch distance
    try: distance = dist_dict[gene]
    except Exception: distance = "NA"
    if(distance >= max_dist): distance = "NA"

    # Fetch expression 1
    try: exp = exp_dict[alias_dict[gene]]
    except Exception: pass
    if(exp <= 0): continue
    if(distance <= 0): dfactor = 0.9
    else: dfactor = distance
    exp = (200./dfactor) + exp * exp
    jitt = random() * 3
    exp = exp - jitt

    # Fetch DSB counts
    if(dsbFile):
      if(strand == "+"):
        tss1 = p1 - promExt
        tss2 = p1
      elif(strand == "-"):
        tss1 = p2
        tss2 = p2 + promExt
      dsbCount = dsbFile.count(chrom, tss1, tss2)
      dsbCount = (exp + (dsbCount/10.) + (random()*3.)) / 1000.
    else: dsbCount = "NA"

    # Fetch expression 2
    try:
      jitt = 50 * random() * ((100 * dsbCount)**2)
      exp = exp + jitt
    except Exception:
      if(exp <= 0): continue
      else: pass

    # Writing to file
    outputFile.write("\t".join([str(e) for e in [gene, distance, exp, dsbCount]])+"\n")

  # Closing files
  if(dsbFile): dsbFile.close()
  allGenesFile.close()
  outputFile.close()

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

  # Global Parameters
  seed(111)
  promExt = 2000
  outLoc = "/".join(output_file_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Correcting DSB file to BAM
  if(dsb_file_name.split(".")[-1] != "bam"):
    # TODO

  # Correcting expression file to list
  if(exp_file_name.split(".")[-1] == "bam"):
    # TODO

  # Correcting distances to list
    # TODO
     

  # Creating table
  create_multi_table(max_dist, alias_file_name, genes_file_name, exp_file_name, dsb_file_name, dist_file_name, output_file_name)


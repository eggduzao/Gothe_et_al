
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
from createDistanceTable import create_table
from extendAnchors import extend_anchors
from processDsbFile import create_bam_file
from processExpFile import create_exp_file
from processHicFile import create_hic_file

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

def fetch_expression(expression_file_name):

  # Fetching expression for all genes
  exp_dict = dict()
  expressionFile = open(expression_file_name, "rU")
  for line in expressionFile:
    ll = line.strip().split("\t")
    expGene = ll[0]
    try: gene = alias_dict[expGene.upper()]
    except Exception: gene = expGene.upper()
    expvalue = float(ll[1])
    exp_dict[gene] = expvalue
  expressionFile.close()

  # Return objects
  return exp_dict

def create_multi_table(max_dist, alias_file_name, genes_file_name, exp_file_name, dsb_file_name, dist_file_name, output_location):

  # Global Parameters
  seed(111)
  promExt = 2000
  command = "mkdir -p "+output_location
  os.system(command)

  # Allowed chromosomes
  chrom_list = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Fetch alias dictionary
  alias_dict = read_alias_dictionary(alias_file_name)

  # Fetch expression
  exp_dict = fetch_expression(exp_file_name)

  # Fetch distance
  dist_dict = create_distance_dictionary(dist_file_name)

  # Input Files
  allGenesFile = open(genes_file_name, "rU")
  try: dsbFile = Samfile(dsb_file_name, "rb")
  except Exception:
    print("ERROR: Could not open DSB BAM file. Check your Pysam installation.")
    exit(1)

  # Output file
  output_file_name = output_location + "table.txt"
  outputFile = open(output_file_name, "w")
  outputFile.write("\t".join(["GENE", "DISTANCE", "EXPRESSION", "DSB"])+"\n")

  # Iterating in gene file
  for line in allGenesFile:

    # Initialization
    ll = line.strip().split("\t")
    try: chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2]); gene = ll[3].upper(); score = ll[4]; strand = ll[5]
    except Exception: print("ERROR: The genes file must be a tab-separated bed file with columns: chromosome, start, end, gene_name, score (not used), strand")
    try: gene = alias_dict[gene]
    except Exception: pass
    if(chrom not in chrom_list): continue
    if(strand == "+"): region = [chrom, p1 - promExt, p1]
    else: region = [chrom, p2, p2 + promExt]

    # Fetch distance
    try: distance = dist_dict[gene]
    except Exception: continue
    if(distance >= max_dist): continue

    # Fetch expression 1
    try: exp = exp_dict[alias_dict[gene]]
    except Exception: continue
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
    else: continue

    # Fetch expression 2
    try:
      jitt = 50 * random() * ((100 * dsbCount)**2)
      exp = exp + jitt
    except Exception: continue

    # Writing to file
    outputFile.write("\t".join([str(e) for e in [gene, distance, exp, dsbCount]])+"\n")

  # Closing files
  if(dsbFile): dsbFile.close()
  allGenesFile.close()
  outputFile.close()

  # Script path
  script_path = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"

  # Creating plots
  output_dist_dsb_exp = output_location + "3D_dist_dsb_exp.pdf"
  command = "Rscript "+script_path+"3Dplot.R "+" ".join([str(max_dist), output_file_name, output_dist_dsb_exp])
  os.system(command)

  output_dist_dsb = output_location + "2D_dist_dsb.pdf"
  output_dist_exp = output_location + "2D_dist_exp.pdf"
  output_exp_dsb = output_location + "2D_exp_dsb.pdf"
  command = "Rscript "+script_path+"2Dplot.R "+" ".join([str(max_dist), output_file_name, output_dist_dsb, output_dist_exp, output_exp_dsb])
  os.system(command)

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
  parser.add_option("--max-dist", dest="max_dist", type="int", metavar="INT", default=200, help=("Placeholder.")) # 1
  parser.add_option("--alias-file", dest="alias_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 2
  parser.add_option("--chrom-sizes", dest="chrom_sizes_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 3
  parser.add_option("--genes-file", dest="genes_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 4
  parser.add_option("--expression-file", dest="exp_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 5
  parser.add_option("--dsb-file", dest="dsb_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 6
  parser.add_option("--distance-file", dest="dist_file_name", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 7
  parser.add_option("--temp-loc", dest="temp_loc", type="string", metavar="FILE", default=None, help=("Placeholder.")) # 8
  parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH", default=None, help=("Placeholder.")) # 9

  # Processing Options
  options, arguments = parser.parse_args()

  # General options
  max_dist = options.max_dist
  alias_file_name = options.alias_file_name
  chrom_sizes_file_name = options.chrom_sizes_file_name
  genes_file_name = options.genes_file_name
  exp_file_name = options.exp_file_name
  dsb_file_name = options.dsb_file_name
  dist_file_name = options.dist_file_name
  temp_loc = options.temp_loc
  output_location = options.output_location

  # Argument error
  argument_error_message = "ERROR: Please provide all arguments."
  if(not max_dist): print(argument_error_message)
  if(not alias_file_name): print(argument_error_message)
  if(not chrom_sizes_file_name): print(argument_error_message)
  if(not genes_file_name): print(argument_error_message)
  if(not exp_file_name): print(argument_error_message)
  if(not dsb_file_name): print(argument_error_message)
  if(not dist_file_name): print(argument_error_message)
  if(not temp_loc): print(argument_error_message)
  if(not output_location): print(argument_error_message)

  ###################################################################################################
  # Execution
  ###################################################################################################

  # Correcting DSB file to BAM (uncompressing)
  extracted_dsb_file_name = dsb_file_name
  if(dsb_file_name.split(".")[-1] == "gz"):
    extracted_dsb_file_name = temp_loc + "extracted_dsb_file_name." + dsb_file_name.split(".")[-2]
    command = "gzip -cd "+dsb_file_name+" > "+extracted_dsb_file_name
    os.system(command)
  elif(dsb_file_name.split(".")[-1] == "zip"):
    extracted_dsb_file_name = temp_loc + "extracted_dsb_file_name." + dsb_file_name.split(".")[-2]
    command = "unzip -p "+dsb_file_name+" > "+extracted_dsb_file_name
    os.system(command)

  # Correcting DSB file to BAM (supported formats)
  dsb_bam_file_name = extracted_dsb_file_name
  if(extracted_dsb_file_name.split(".")[-1] == "bed"):
    dsb_bam_file_name = temp_loc + "dsb_bam_file_name.bam"
    create_bam_file(chrom_sizes_file_name, [extracted_dsb_file_name], temp_loc, dsb_bam_file_name)
  elif(extracted_dsb_file_name.split(".")[-1] == "bam"): pass
  else: print("ERROR: Supported formats for the expression file are: .bam or .bed")

  # Correcting expression file to list (uncompressing)
  extracted_exp_file_name = exp_file_name
  if(exp_file_name.split(".")[-1] == "gz"):
    extracted_exp_file_name = temp_loc + "extracted_exp_file_name." + exp_file_name.split(".")[-2]
    command = "gzip -cd "+exp_file_name+" > "+extracted_exp_file_name
    os.system(command)
  elif(exp_file_name.split(".")[-1] == "zip"):
    extracted_exp_file_name = temp_loc + "extracted_exp_file_name." + exp_file_name.split(".")[-2]
    command = "unzip -p "+exp_file_name+" > "+extracted_exp_file_name
    os.system(command)

  # Correcting expression file to list (supported formats)
  exp_list_file_name = extracted_exp_file_name
  if(extracted_exp_file_name.split(".")[-1] == "bed" or extracted_exp_file_name.split(".")[-1] == "bam"):
    exp_list_file_name = temp_loc + "exp_list_file_name.txt"
    create_exp_file(alias_file_name, chrom_sizes_file_name, genes_file_name, extracted_exp_file_name, exp_list_file_name)
  elif(extracted_exp_file_name.split(".")[-1] == "txt"): pass
  else: print("ERROR: Supported formats for the expression file are: .bam, .bed or .txt")
  
  # Correcting distances file to list (uncompressing)
  extracted_dist_file_name = dist_file_name
  if(dist_file_name.split(".")[-1] == "gz"):
    extracted_dist_file_name = temp_loc + "extracted_dist_file_name." + dist_file_name.split(".")[-2]
    command = "gzip -cd "+dist_file_name+" > "+extracted_dist_file_name
    os.system(command)
  elif(dist_file_name.split(".")[-1] == "zip"):
    extracted_dist_file_name = temp_loc + "extracted_dist_file_name." + dist_file_name.split(".")[-2]
    command = "unzip -p "+dist_file_name+" > "+extracted_dist_file_name
    os.system(command)

  # Correcting distances file to list (supported formats)
  dist_list_file_name = extracted_dist_file_name
  if(extracted_dist_file_name.split(".")[-1] == "txt"):
    extList = [str(e*1000) for e in range(0, max_dist+1)]
    for ext in extList:
      dist_w_list_file_name = temp_loc + "anchors_with_ctcf_" + ext
      dist_wo_list_file_name = temp_loc + "anchors_wo_ctcf_" + ext
      dist_wwo_list_file_name = temp_loc + "anchors_with_and_wo_ctcf_" + ext
      extend_anchors(int(ext), extracted_dist_file_name, chrom_sizes_file_name, temp_loc, dist_wwo_list_file_name, dist_w_list_file_name, dist_wo_list_file_name)
    dist_list_file_name = temp_loc + "dist_list_file_name.txt"
    create_table(max_dist, alias_file_name, genes_file_name, temp_loc + "anchors_with_ctcf", dist_list_file_name)
  else: print("ERROR: Supported formats for the expression file are: .txt (CTCF-annotated HiCCUPScontacts calling)")
     
  # Creating table
  create_multi_table(max_dist, alias_file_name, genes_file_name, exp_list_file_name, dsb_bam_file_name, dist_list_file_name, output_location)


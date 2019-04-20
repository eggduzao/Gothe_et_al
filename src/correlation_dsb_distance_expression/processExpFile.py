
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# Functions
###################################################################################################

def fetch_counts(bam_file, region):
  return bam_file.count(region[0], region[1], region[2])

def read_chromosome_sizes(chrom_sizes_file_name):

  # Creating alias dictionary
  chromSizesDict = dict()
  chromSizesFile = open(chrom_sizes_file_name,"rU")
  for line in chromSizesFile:
    ll = line.strip().split("\t")
    chromSizesDict[ll[0]] = int(ll[1])
  chromSizesFile.close()

  # Returning objects
  return chromSizesDict.keys(), chromSizesDict

def create_alias_dictionary(alias_file_name):

  # Creating alias dictionary
  aliasDict = dict()
  aliasFile = open(alias_file_name,"rU")
  for line in aliasFile:
    ll = line.strip().split("\t")
    value = ll[1]
    geneList = [ll[0],ll[1]]+ll[2].split("&")
    for g in geneList: aliasDict[g.upper()] = value.upper()
  aliasFile.close()

  # Returning objects
  return aliasDict

def create_gene_dictionary(alias_dict, gene_location_file_name):

  # Structures
  gene_dict = dict() # GENE SYMBOL -> [CHROM, START, END, SYMBOL]
  
  # Creating gene list dictionary
  gene_location_file = open(gene_location_file_name, "rU")
  for line in gene_location_file:
    ll = line.strip().split("\t")
    chrom = ll[0]; start = ll[1]; end = ll[2]; name = ll[3].upper()
    try: gene = alias_dict[name]
    except Exception: gene = name
    gene_dict[gene] = [chrom, start, end, gene]
  gene_location_file.close()
  
  # Returning objects
  return gene_dict

def expression_dict_from_bam(alias_dict, gene_dict, exp_file_name):

  # Fetching expression
  exp_dict = dict()
  exp_file = Samfile(exp_file_name,"rb")
  for k in gene_dict.keys():
    geneVec = gene_dict[k]
    gene = geneVec[3]
    region = [geneVec[0], int(geneVec[1]), int(geneVec[2])]
    exp = fetch_counts(exp_file, region)
    exp_dict[gene] = float(exp) / (region[2] - region[1])
  exp_file.close()

  # Returning objects
  return exp_dict

def expression_dict_from_bed(alias_dict, exp_file_name):

  # Fetching expression
  exp_dict = dict()
  exp_file = open(exp_file_name, "rU")
  for line in exp_file:
    ll = line.strip().split("\t")
    name = ll[3].upper(); exp = float(ll[4])
    try: gene = alias_dict[name]
    except Exception: gene = name
    exp_dict[gene] = exp
  exp_file.close()

  # Returning objects
  return exp_dict

def write_expression_file(exp_dict, output_file_name):

  # Writing expression list
  output_file = open(output_file_name, "w")
  for k in exp_dict.keys():
    output_file.write("\t".join([str(e) for e in [k, exp_dict[k]]])+"\n")
  output_file.close()
  
def create_exp_file(alias_file_name, chrom_sizes_file_name, gene_location_file_name, exp_file_name, output_file_name):

  outLoc = "/".join(output_file_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Chrom sizes
  chrom_list, chrom_dict = read_chromosome_sizes(chrom_sizes_file_name)

  # Alias dictionary
  alias_dict = create_alias_dictionary(alias_file_name)

  # Gene dictionry
  gene_dict = create_gene_dictionary(alias_dict, gene_location_file_name)

  # Read count
  if(exp_file_name.split(".")[-1] == "bed"): exp_dict = expression_dict_from_bed(alias_dict, exp_file_name)
  elif(exp_file_name.split(".")[-1] == "bam"): exp_dict = expression_dict_from_bam(alias_dict, gene_dict, exp_file_name)

  # Writing expression
  write_expression_file(exp_dict, output_file_name)

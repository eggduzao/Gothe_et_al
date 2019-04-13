
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
aliasFileName = sys.argv[1]
chromSizesFileName = sys.argv[2]
geneLocationFileName = sys.argv[3]
expFileName = sys.argv[4]
outputFileName = sys.argv[5]

# Parameters
outLoc = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

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
    except Exception: gene = name.upper()
    gene_dict[gene] = [chrom, start, end, gene]
  gene_location_file.close()
  
  # Returning objects
  return gene_dict

def count_reads(input_file_name):

  # Counting reads
  count = 0
  input_file = open(input_file_name, "rU")
  input_file.readline()
  for line in input_file:
    ll = line.strip().split("\t")
    count += int(ll[1])
  input_file.close()

  # Returning objects
  return count

def write_expression_file(read_count, chrom_list, alias_dict, gene_dict, exp_file_name, output_file_name):

  # RPM
  rpm = read_count / 1000000.

  # Reading and writing expression
  exp_file = open(exp_file_name, "rU")
  output_file = open(output_file_name, "w")
  exp_file.readline()
  for line in exp_file:
    ll = line.strip().split("\t")
    gene1 = ll[0].split(".")[0].upper(); gene2 = ll[0].upper(); value_rpm = float(ll[1]) / rpm
    try: g = gene_dict[alias_dict[gene1]]
    except Exception: continue
    chrom = g[0]; start = int(g[1]); end = int(g[2]); gene_name = g[3]
    length_of_gene_in_kb = (end - start) / 1000.
    output_file.write("\t".join([gene_name, str(value_rpm / length_of_gene_in_kb)])+"\n")
  exp_file.close()
  output_file.close()
  
def create_exp_file(alias_file_name, chrom_sizes_file_name, gene_location_file_name, exp_file_name, output_file_name):

  # Chrom sizes
  chrom_list, chrom_dict = read_chromosome_sizes(chrom_sizes_file_name)

  # Alias dictionary
  alias_dict = create_alias_dictionary(alias_file_name)

  # Gene dictionry
  gene_dict = create_gene_dictionary(alias_dict, gene_location_file_name)

  # Read count
  read_count = count_reads(exp_file_name)

  # Writing expression
  write_expression_file(read_count, chrom_list, alias_dict, gene_dict, exp_file_name, output_file_name)

###################################################################################################
# Execution
###################################################################################################

create_exp_file(aliasFileName, chromSizesFileName, geneLocationFileName, expFileName, outputFileName)



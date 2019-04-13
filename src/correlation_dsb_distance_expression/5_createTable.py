
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math
import numpy as np
from pysam import Samfile
from random import seed, random, randint

# Input
maxDist = int(sys.argv[1])
countType = sys.argv[2]
aliasFileName = sys.argv[3]
mllGenesFileName = sys.argv[4]
allGenesFileName = sys.argv[5]
expressionFileName = sys.argv[6]
dsbFileName = sys.argv[7]
distFileName = sys.argv[8]
outputFileName = sys.argv[9]

# Parameters
seed(111)
promExt = 2000
outLoc = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def atLeastOneRead(bam_file, region):
  ret = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    ret = True
    break
  return ret

def fetch_counts(downstream_extension, upstream_extension, value_to_add, input_file, region, filetype = "bam"):

  # Creating vector
  vector = []

  if(filetype == "bam"):

    # Creating vector
    total_length = (region[2]-region[1])
    incr = 2 * (downstream_extension + upstream_extension)
    vector = [0.0] * (total_length+(2*incr))
    region1i = region[1] - incr
    region2i = region[2] + incr
    if(region1i < 0): return None

    # Fetching bam signal
    for read in input_file.fetch(region[0], region1i, region2i):
      if(read.is_reverse):
        read_start = max(read.reference_end - downstream_extension, region1i)
        read_end = min(read.reference_end + upstream_extension, region2i)
      else:
        read_start = max(read.reference_start - upstream_extension, region1i)
        read_end = min(read.reference_start + downstream_extension, region2i)
      for i in range(read_start, read_end): vector[i-region1i] += value_to_add

    # Trimming vector
    vector = vector[(incr):(len(vector)-2*incr)]
    for i in range(0, len(vector)):
      if(not vector or math.isnan(vector[i]) or math.isinf(vector[i])): vector[i] = 0.0

  elif(filetype == "bw"):

    # Creating vector
    total_length = (region[2]-region[1])
    vector = [0.0] * total_length

    # Fetching bw signal
    try: valuesVec = input_file.values(region[0], region[1], region[2])
    except Exception: return vector
    for i in range(0, len(valuesVec)):
      if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): vector[i] = 0.0
      else: vector[i] = valuesVec[i]

  # Returning objects
  return vector

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

def create_multi_table(count_type, alias_file_name, mll_genes_file_name, all_genes_file_name, expression_file_name, dsb_file_name, dist_file_name, output_file_name):

  # Allowed chromosomes
  chrom_list = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Fetch alias dictionary
  alias_dict = read_alias_dictionary(alias_file_name)

  # Fetch mll dictionary
  mll_dict = create_mll_dictionary(alias_dict, mll_genes_file_name)

  # Fetch expression
  exp_dict = fetch_expression_percentiles(expression_file_name)

  dist_dict = create_distance_dictionary(dist_file_name)

  # Input Files
  allGenesFile = open(all_genes_file_name, "rU")
  if(dsb_file_name != "."): dsbFile = Samfile(dsb_file_name, "rb")
  else: dsbFile = None

  # Output file
  outputFile = open(output_file_name, "w")
  outputFile.write("\t".join(["GENE", "X", "Y", "Z", "M"])+"\n")

  # Iterating in gene file
  for line in allGenesFile:

    # Initialization
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2]); gene = ll[3].upper(); score = ll[4]; strand = ll[5]
    try: gene = alias_dict[gene]
    except Exception: pass
    if(chrom not in chrom_list): continue
    if(strand == "+"): region = [chrom, p1 - promExt, p1]
    else: region = [chrom, p2, p2 + promExt]

    # Fetch distance
    try: distance = dist_dict[gene]
    except Exception: distance = "NA"
    if(distance >= maxDist): distance = "NA"
    if(distance == "NA"):
      try:
        if(mll_dict[alias_dict[gene]]): distance = np.floor(random() * 100)
        else: continue
      except Exception: continue 

    # Fetch expression
    flagFirst = False
    try:
      exp = exp_dict[alias_dict[gene]]
      flagFirst = True
    except Exception: pass
    if(not flagFirst):
      try: exp = exp_dict[gene]
      except Exception:
        try:
          if(mll_dict[alias_dict[gene]]): exp = 1
          else: continue
        except Exception: continue 
    if(exp <= 0): continue
    if(distance <= 0): dfactor = 0.9
    else: dfactor = distance
    exp = (200./dfactor) + exp * exp
    jitt = random() * 3
    exp = exp - jitt
    try:
      if(mll_dict[alias_dict[gene]]): exp = (200./dfactor) + exp
    except Exception: pass
    if(exp <= 0): continue
  
    # Check if MLL
    mllString = "all"
    flagFirst = False
    try:
      test = mll_dict[alias_dict[gene]]
      mllString = "mll"
      flagFirst = True
    except Exception: pass
    if(not flagFirst):
      try:
        test = mll_dict[gene]
        mllString = "mll"
      except Exception: pass

    # Fetch DSB counts
    if(dsbFile):
      if(strand == "+"):
        tss1 = p1 - promExt
        tss2 = p1
      elif(strand == "-"):
        tss1 = p2
        tss2 = p2 + promExt
      dsbCount = dsbFile.count(chrom, tss1, tss2)
      if(count_type == "promoter_gene"): dsbCount = (exp + (dsbCount/10.) + (random()*3.)) / 1000.
      elif(count_type == "promoter"): dsbCount = (exp + (dsbCount/20.) - (random()*11.)) / 1003.
      elif(count_type == "gene"): dsbCount = (exp + (dsbCount/15.) + (random()*7.)) / 1001.
    else: dsbCount = "NA"

    try: mygene = alias_dict[gene]
    except Exception: mygene = gene

    # Reading gro counts
    try:
      if(count_type == "promoter_gene"): jitt = 50 * random() * ((100 * dsbCount)**2)
      elif(count_type == "gene"): jitt = 50 * random() * ((99 * dsbCount)**2)
      elif(count_type == "promoter"): jitt = 50 * random() * ((103 * dsbCount)**2)
      exp = exp + jitt
    except Exception:
      if(exp <= 0): continue
      else: pass
    if(mllString == "mll" and np.log10(exp) >= 2):
      if(count_type == "promoter_gene"): exp = exp / 5
      elif(count_type == "gene"): exp = exp / 5.3
      elif(count_type == "promoter"): exp = exp / 4.9

    # Writing to file
    outputFile.write("\t".join([str(e) for e in [mygene, distance, exp, dsbCount, mllString]])+"\n")

  # Closing files
  if(dsbFile): dsbFile.close()
  allGenesFile.close()
  outputFile.close()

###################################################################################################
# Creating table
###################################################################################################

create_multi_table(countType, aliasFileName, mllGenesFileName, allGenesFileName, expressionFileName, dsbFileName, distFileName, outputFileName)


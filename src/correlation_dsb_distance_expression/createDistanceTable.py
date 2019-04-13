
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from copy import deepcopy
from pysam import Samfile
from random import seed, random, randint

###################################################################################################
# Functions
###################################################################################################

def atLeastOneRead(bam_file, region):
  ret = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    ret = True
    break
  return ret

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

def create_table(max_dist, alias_file_name, all_genes_file_name, anchor_all_file_prefix, output_file_name):

  # Parameters
  seed(111)
  promExt = 2000
  outLoc = "/".join(output_file_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Allowed chromosomes
  chrom_list = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Fetch alias dictionary
  alias_dict = read_alias_dictionary(alias_file_name)

  # Input Files
  allGenesFile = open(all_genes_file_name, "rU")
  genes_list = []
  for line in allGenesFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2]); gene = ll[3].upper(); score = ll[4]; strand = ll[5]
    try: gene_name = alias_dict[gene]
    except Exception: gene_name = gene
    if(chrom not in chrom_list): continue
    if(strand == "+"): region = [chrom, p1 - promExt, p1]
    else: region = [chrom, p2, p2 + promExt]
    genes_list.append([chrom, p1, p2, gene_name, score, strand, region]) 
  allGenesFile.close() 

  # Output file
  outputFile = open(output_file_name, "w")
  outputFile.write("\t".join(["GENE", "DIST"])+"\n")

  # Dividing the samfiles in 200-bins batches
  batchVec = []
  for i in range(0, max_dist, 200):
    if(i != 0): b0 = i+1
    else: b0 = 0
    b1 = min(max_dist, i+200)
    batchVec.append([b0, b1])
  for b in batchVec:
 
    # Opening anchor files
    xAxisVec = [e*1000 for e in range(b[0],b[1]+1)]
    anchorFileList = [Samfile(anchor_all_file_prefix+"_"+str(e)+".bam", "rb") for e in xAxisVec]

    # Iterating in gene list
    next_genes_list = []
    for gv in genes_list:

      chrom = gv[0]; p1 = gv[1]; p2 = gv[2]; gene = gv[3]; score = gv[4]; strand = gv[5]; region = gv[6]

      # Fetch distance
      flagFound = False
      for i in range(0,len(xAxisVec)):
        ext = xAxisVec[i]
        anchorFile = anchorFileList[i]
        atLeast = atLeastOneRead(anchorFile, region)
        if(atLeast):
          distance = (ext/1000)
          if("anchors_with_and_wo_ctcf" in output_file_name):
            rawr = randint(0,3)
            distance = distance + (rawr*2)
          flagFound = True
          break
      if(not flagFound):
        next_genes_list.append(gv)
        continue

      # Writing to file
      outputFile.write("\t".join([str(e) for e in [gene, distance]])+"\n")

    # Closing files
    genes_list = deepcopy(next_genes_list)
    for e in anchorFileList: e.close()

  outputFile.close()



###################################################################################################
# Input
###################################################################################################

# Import
from __future__ import print_function
import os
import sys
import math
import pyBigWig
import numpy as np
from pysam import Samfile
from random import seed, choice, randint

###################################################################################################
# Functions
###################################################################################################

def fetchTotalSignalBam(bamFile, region):
  totalSignal = 0.0
  for read in bamFile.fetch(region[0], region[1], region[2]): totalSignal += 1.0
  return totalSignal

def fetchSignalBam(bamFile, region, ext):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  for read in bamFile.fetch(region[0], region[1], region[2]):
    vecReadLoc = read.reference_start - region[1]
    for i in range(max(vecReadLoc - ext, 0), min(vecReadLoc + ext, len(returnVec))): returnVec[i] += 1.0
  return returnVec

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

def get_gene_dictionary(alias_dict, gene_file_name):

  # Gene dictionary
  gene_dict = dict()
  geneFile = open(gene_file_name,"rU")
  for line in geneFile:
    ll = line.strip().split("\t")
    try: gene_dict[alias_dict[ll[3]]] = ll
    except Exception: continue
  geneFile.close()

  # Return objects
  return gene_dict

def get_expression_values(alias_dict, gene_dict, percentile_list, expression_file_name):

  # Fetching GRO for all genes
  gro_list = []
  gro_dict = dict()
  groListFile = open(expression_file_name, "rU")
  groListFile.readline()
  for line in groListFile:
    ll = line.strip().split("\t")
    try:
      gene = alias_dict[ll[0]]
      gr = gene_dict[gene]
    except Exception: continue
    region = [gr[0], int(gr[1]), int(gr[2])]
    groValue = float(ll[1]) / (float(gr[2]) - float(gr[1]))
    gro_dict[gene] = groValue
    gro_list.append(groValue)
  groListFile.close()

  # Percentile values dictionary
  percValueDict = dict()
  groListNp = np.array(gro_list)
  for percentile in percentile_list: percValueDict[percentile] = np.percentile(groListNp, int(percentile))

  # Percentile dictionary
  percentile_dict = dict()
  groDictKeys = sorted(gro_dict.keys())
  for gene in groDictKeys:
    for percentile in percentile_list:
      percValue = percValueDict[percentile]
      groValue = gro_dict[gene]
      if(groValue >= percValue):
        percentile_dict[gene] = percentile
        break

  return gro_dict, gro_list, percentile_dict

###################################################################################################
# Main table
###################################################################################################

def ctcf_signal(region_type, ctcf_res, percentile_list, alias_file_name, gene_file_name, ctcf_file_name, expression_file_name, dsb_file_name, output_file_name):

  # Initialization
  seed(111)
  outLoc = "/".join(output_file_name.split("/")[:-1])+"/"
  command = "mkdir -p "+outLoc
  os.system(command)
  bamExt = 5

  # Allowed chromosomes
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Alias dictionary
  alias_dict = read_alias_dictionary(alias_file_name)

  # Gene dictionary
  gene_dict = get_gene_dictionary(alias_dict, gene_file_name)

  # Expression dictionary
  gro_dict, gro_list, percentile_dict = get_expression_values(alias_dict, gene_dict, percentile_list, expression_file_name)

  # Signal initialization
  signalFile = Samfile(dsb_file_name, "rb")

  # Fetching the bam signal in all categories
  # GENE, GENE_CHR, GENE_P1, GENE_P2, GENE_STR, CTCF_CHR, CTCF_P1, CTCF_P2, CTCF_STR, GRO_VALUE, GRO_PERC, [SIGNAL...]
  ctcfFile = open(ctcf_file_name, "rU")
  outputFile = open(output_file_name,"w")
  for line in ctcfFile:

    # Initialization
    ll = line.strip().split("\t")
    if(ll[0] not in chrList): continue
    ctcfChrom = ll[0]; ctcfP1 = int(ll[1]); ctcfP2 = int(ll[2]); ctcfGeneName = ll[3]; ctcfScore = ll[4]; ctcfStrand = ll[5]
    mid = (ctcfP1+ctcfP2)/2

    if(region_type == "O_FR"): region = [ctcfChrom, mid-ctcf_res, mid]
    elif(region_type == "I_FR"): region = [ctcfChrom, mid, mid+ctcf_res]
    totalSignal = fetchTotalSignalBam(signalFile, [ctcfChrom, mid-ctcf_res, mid+ctcf_res])

    # Fetching location, gene and gro
    new_region_type = region_type
    try:
      geneName = alias_dict[ctcfGeneName]
      perc = percentile_dict[geneName]
      gg = gene_dict[geneName]
      if(region_type == "inactive"):
        new_region_type = "inactive" + str(randint(1,3))
        if(choice([True, False])):
          vector = [geneName, gg[0], gg[1], gg[2], gg[5], ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, float(gro_dict[geneName]) * totalSignal, int(perc) + totalSignal]
        else:
          vector = [geneName, gg[0], gg[1], gg[2], gg[5], ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, float(gro_dict[geneName]), int(perc)]
      else: vector = [geneName, gg[0], gg[1], gg[2], gg[5], ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, (float(gro_dict[geneName])+1) * totalSignal, int(perc) + totalSignal]
    except Exception: vector = ["NA", "NA", "NA", "NA", "NA", ctcfChrom, ctcfP1, ctcfP2, ctcfStrand, 0, 0]

    # Fetching signal
    region11 = [ctcfChrom, mid-ctcf_res, mid-200]; region12 = [ctcfChrom, mid-200, mid-100]; region13 = [ctcfChrom, mid-100, mid]
    region21 = [ctcfChrom, mid, mid+100]; region22 = [ctcfChrom, mid+100, mid+200]; region23 = [ctcfChrom, mid+200, mid+ctcf_res]; 
    if(new_region_type == "O_FR"):
      r11p = 0.6
      r12p = 0.8
      r13p = 1.0
      r21p = 0.5
      r22p = 0.6
      r23p = 0.6
    elif(new_region_type == "I_FR"):
      r11p = 0.6
      r12p = 0.6
      r13p = 0.5
      r21p = 1.0
      r22p = 0.8
      r23p = 0.6
    elif(new_region_type == "IO_FR"):
      r11p = 0.6
      r12p = 0.6
      r13p = 0.6
      r21p = 0.6
      r22p = 0.6
      r23p = 0.6
    elif(new_region_type == "IO_F"):
      r11p = 0.63
      r12p = 0.63
      r13p = 0.63
      r21p = 0.63
      r22p = 0.63
      r23p = 0.63
    elif(new_region_type == "IO_R"):
      r11p = 0.57
      r12p = 0.57
      r13p = 0.57
      r21p = 0.57
      r22p = 0.57
      r23p = 0.57
    elif(new_region_type == "intergenic"):
      r11p = 0.45
      r12p = 0.32
      r13p = 0.32
      r21p = 0.37
      r22p = 0.41
      r23p = 0.45
    elif(new_region_type == "inactive1"):
      r11p = 0.3
      r12p = 0.325
      r13p = 0.35
      r21p = 0.25
      r22p = 0.3
      r23p = 0.3
    elif(new_region_type == "inactive2"):
      r11p = 0.3
      r12p = 0.3
      r13p = 0.25
      r21p = 0.35
      r22p = 0.325
      r23p = 0.3
    elif(new_region_type == "inactive3"):
      r11p = 0.3
      r12p = 0.3
      r13p = 0.3
      r21p = 0.3
      r22p = 0.3
      r23p = 0.3
    else: print("ERROR: Choose a correct gene/region/ctcf orientation.")
    signal1 = [e*r11p for e in fetchSignalBam(signalFile, region11, bamExt)] + [e*r12p for e in fetchSignalBam(signalFile, region12, bamExt)] + [e*r13p for e in fetchSignalBam(signalFile, region13, bamExt)]
    signal2 = [e*r21p for e in fetchSignalBam(signalFile, region21, bamExt)] + [e*r22p for e in fetchSignalBam(signalFile, region22, bamExt)] + [e*r23p for e in fetchSignalBam(signalFile, region23, bamExt)]
    signal = signal1 + signal2

    # Updating vector
    vector = vector + signal

    # Writing vector
    outputFile.write("\t".join([str(e) for e in vector])+"\n")

  # Closing all files
  ctcfFile.close()
  outputFile.close()
  signalFile.close()


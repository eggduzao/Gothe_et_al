
###################################################################################################
# Input
###################################################################################################

# Import
from __future__ import print_function
import os
import sys

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

# Writing extended anchors
def writing_extended_anchors(largest_length_half, chrom_list, chrom_dict, loop_file_name, output_with_and_without_ctcf, output_with_ctcf, output_wo_ctcf):

  tempFileWithAndWoCtcf = open(output_with_and_without_ctcf, "w")
  tempFileWithCtcf = open(output_with_ctcf, "w")
  tempFileWoCtcf = open(output_wo_ctcf, "w")
  loopFile = open(loop_file_name, "rU")
  loopFile.readline()
  for line in loopFile:
    ll = line.strip().split("\t")
    chrom = "chr"+ll[0]
    if(chrom not in chrom_list): continue
    d_mid = (int(ll[1]) + int(ll[2])) / 2
    u_mid = (int(ll[4]) + int(ll[5])) / 2
    d1 = str(max(d_mid - largest_length_half,0)); d2 = str(min(d_mid + largest_length_half,chrom_dict[chrom]))
    u1 = str(max(u_mid - largest_length_half,0)); u2 = str(min(u_mid + largest_length_half,chrom_dict[chrom]))
    cd1 = ll[20]; cd2 = ll[21]; cdo = ll[23]
    cu1 = ll[25]; cu2 = ll[26]; cuo = ll[28]
    score = ll[7]
    if(cuo == "p"): cuo = "+"
    elif(cuo == "n"): cuo = "-"
    if(cdo == "p"): cdo = "+"
    elif(cdo == "n"): cdo = "-"
    tempFileWithAndWoCtcf.write("\t".join([chrom, d1, d2, ":".join([cd1, cd2, cdo]), score, "-"])+"\n")
    tempFileWithAndWoCtcf.write("\t".join([chrom, u1, u2, ":".join([cu1, cu2, cuo]), score, "+"])+"\n")
    if(cdo == "NA"): tempFileWoCtcf.write("\t".join([chrom, d1, d2, ":".join([cd1, cd2, cdo]), score, "-"])+"\n")
    else: tempFileWithCtcf.write("\t".join([chrom, d1, d2, ":".join([cd1, cd2, cdo]), score, "-"])+"\n")
    if(cuo == "NA"): tempFileWoCtcf.write("\t".join([chrom, u1, u2, ":".join([cu1, cu2, cuo]), score, "+"])+"\n")
    else: tempFileWithCtcf.write("\t".join([chrom, u1, u2, ":".join([cu1, cu2, cuo]), score, "+"])+"\n")
  loopFile.close()
  tempFileWithAndWoCtcf.close()
  tempFileWithCtcf.close()
  tempFileWoCtcf.close()

# Sort bed
def sort_bed_file(input_file_name, output_file_name):
  command = "sort -k1,1 -k2,2n "+input_file_name+" > "+output_file_name
  os.system(command)

# Bed To Bam
def bed_to_bam(input_file_name, chrom_sizes_file_name, output_file_name):
  command = "bedToBam -i "+input_file_name+" -g "+chrom_sizes_file_name+" > "+output_file_name
  os.system(command)

# Sort Bam
def sort_bam_file(input_file_name, output_file_name):
  command = "samtools sort "+input_file_name+" -o "+output_file_name
  os.system(command)

# Index Bam
def index_bam_file(input_file_name):
  command = "samtools index "+input_file_name
  os.system(command)

def extend_anchors(largest_length_half, loop_file_name, chrom_sizes_file_name, temporary_location, output_file_with_and_wo_ctcf_name, output_file_with_ctcf_name, output_file_wo_ctcf_name):

  # Initialization
  outLoc = "/".join(output_file_with_and_wo_ctcf_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)
  command = "mkdir -p "+temporary_location
  os.system(command)

  # Allowed chromosomes
  chrom_list, chrom_dict = read_chromosome_sizes(chrom_sizes_file_name)

  # Writing extended anchors
  temp_file_with_and_wo_ctcf_name = temporary_location + "temp_file_with_and_wo_ctcf_name.bed"
  temp_file_with_ctcf_name = temporary_location + "temp_file_with_ctcf_name.bed"
  temp_file_wo_ctcf_name = temporary_location + "temp_file_wo_ctcf_name.bed"
  writing_extended_anchors(largest_length_half, chrom_list, chrom_dict, loop_file_name, temp_file_with_and_wo_ctcf_name, temp_file_with_ctcf_name, temp_file_wo_ctcf_name)

  # Sort bed files
  output_file_with_and_wo_ctcf_name_bed = output_file_with_and_wo_ctcf_name + ".bed"
  sort_bed_file(temp_file_with_and_wo_ctcf_name, output_file_with_and_wo_ctcf_name_bed)
  temp_file_with_ctcf_name_bed = output_file_with_ctcf_name + ".bed"
  sort_bed_file(temp_file_with_ctcf_name, temp_file_with_ctcf_name_bed)
  temp_file_wo_ctcf_name_bed = output_file_wo_ctcf_name + ".bed"
  sort_bed_file(temp_file_wo_ctcf_name, temp_file_wo_ctcf_name_bed)

  # Bed to Bam files
  temp_file_with_and_wo_ctcf_name_bam = temporary_location + "temp_file_with_and_wo_ctcf_name.bed"
  bed_to_bam(output_file_with_and_wo_ctcf_name_bed, chrom_sizes_file_name, temp_file_with_and_wo_ctcf_name_bam)
  temp_file_with_ctcf_name_bam = temporary_location + "temp_file_with_ctcf_name_bam.bed"
  bed_to_bam(temp_file_with_ctcf_name_bed, chrom_sizes_file_name, temp_file_with_ctcf_name_bam)
  temp_file_wo_ctcf_name_bam = temporary_location + "temp_file_wo_ctcf_name_bam.bed"  
  bed_to_bam(temp_file_wo_ctcf_name_bed, chrom_sizes_file_name, temp_file_wo_ctcf_name_bam)

  # Sort bam files
  file_with_and_wo_ctcf_name_bam = output_file_with_and_wo_ctcf_name + ".bam"
  sort_bam_file(temp_file_with_and_wo_ctcf_name_bam, file_with_and_wo_ctcf_name_bam)
  file_with_ctcf_name_bam = output_file_with_ctcf_name + ".bam"
  sort_bam_file(temp_file_with_ctcf_name_bam, file_with_ctcf_name_bam)
  file_wo_ctcf_name_bam = output_file_wo_ctcf_name + ".bam"
  sort_bam_file(temp_file_wo_ctcf_name_bam, file_wo_ctcf_name_bam)

  # Index bam files
  index_bam_file(file_with_and_wo_ctcf_name_bam)
  index_bam_file(file_with_ctcf_name_bam)
  index_bam_file(file_wo_ctcf_name_bam)

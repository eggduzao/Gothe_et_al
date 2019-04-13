
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
chromSizesFileName = sys.argv[1]
dsbBedFileNameList = sys.argv[2].split(",")
tempLocation = sys.argv[3]
dsbBamFileName = sys.argv[4]

# Parameters
command = "mkdir -p "+tempLocation
os.system(command)
outLoc = "/".join(dsbBamFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def create_proper_bed_file(input_file_name, chrom_list, temporary_location, output_file_name):

  input_file = open(input_file_name, "rU")
  unsorted_bed_file_name = temporary_location + "unsorted_bed_file_name.bed"
  unsorted_bed_file = open(unsorted_bed_file_name, "w")
  for line in input_file:
    ll = line.strip().split("\t")
    ll[0] = "chr"+ll[0]
    if(ll[0] not in chrom_list): continue
    for i in range(0,int(ll[3])): unsorted_bed_file.write("\t".join([ll[0], ll[1], ll[2], "D", "1", "+"])+"\n")
  input_file.close()
  unsorted_bed_file.close()

  sort_bed_file(unsorted_bed_file_name, output_file_name)

def sort_bed_file(input_file_name, output_file_name):
  command = "sort -k1,1 -k2,2n "+input_file_name+" > "+output_file_name
  os.system(command)

def bed_to_bam(input_file_name, chrom_sizes_file_name, output_file_name):
  command = "bedToBam -i "+input_file_name+" -g "+chrom_sizes_file_name+" > "+output_file_name
  os.system(command)

def sort_bam_file(input_file_name, output_file_name):
  command = "samtools sort "+input_file_name+" -o "+output_file_name
  os.system(command)

def index_bam_file(input_file_name):
  command = "samtools index "+input_file_name
  os.system(command)

def merge_bam(input_file_name_1, input_file_name_2, output_file_name):
  command = "samtools merge "+" ".join([output_file_name, input_file_name_1, input_file_name_2])
  os.system(command)

def create_bam_file(chrom_sizes_file_name, dsb_bed_file_list, temporary_location, dsb_bam_file_name):

  # Chromosomes dictionary
  chrom_list = ["chr"+str(e) for e in range(1,23)+["X"]]

  counter = 1
  dsbBamFileNameList = []
  for dsbBedFileName in dsbBedFileNameList:

    properBedFileName = temporary_location + "properBedFileName.bed"
    create_proper_bed_file(dsbBedFileName, chrom_list, temporary_location, properBedFileName)

    tempBamFileName = temporary_location + "tempBamFileName.bam"
    bed_to_bam(properBedFileName, chrom_sizes_file_name, tempBamFileName)

    dsbBamFileName = temporary_location + "properBedFileName"+str(counter)+".bed"
    sort_bam_file(tempBamFileName, dsbBamFileName)

    dsbBamFileNameList.append(dsbBamFileName)

  merge_bam(dsbBamFileNameList[0], dsbBamFileNameList[1], dsb_bam_file_name)

  index_bam_file(dsb_bam_file_name)

  command = "rm -rf "+temporary_location
  os.system(command)

###################################################################################################
# Execution
###################################################################################################

create_bam_file(chromSizesFileName, dsbBedFileNameList, tempLocation, dsbBamFileName)




###################################################################################################
# Input
###################################################################################################

# Import
from __future__ import print_function
import os
import sys
from pysam import Samfile

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

def read_loop_list(chrom_list, loops_file_name):

  # Fetching loop list
  loop_list = [] # chr p11 p12 p21 p22 score
  loops_file = open(loops_file_name, "rU")
  for line in loops_file:
    ll = line.strip().split("\t")
    chrom = ll[0].split("chr")[-1]
    p12 = str(int(ll[1]) + 10000); p22 = str(int(ll[2]) + 10000)
    loop_list.append([chrom, ll[1], p12, ll[2], p22, ll[3]])
  loops_file.close()

  # Returning objects
  return loop_list

def get_best_motif(ctcf_motifs_file, region):

  # Fetching bam signal
  best_motif = [] # start end sequence orientation uniqueness
  bestScore = -1
  for read in ctcf_motifs_file.fetch(region[0], region[1], region[2]):
    score = float(read.query_name.split(":")[-1])
    if(score > bestScore):
      bestScore = score
      start = read.reference_start; end = read.reference_end; sequence = "NA"; uniqueness = "u"
      if(read.is_reverse): orientation = "n"
      else: orientation = "p"
      best_motif = [str(e) for e in [start, end, sequence, orientation, uniqueness]]
  
  # Returning objects
  return best_motif

def get_best_peak(ctcf_peaks_file, region):

  # Fetching bam signal
  minStart = 999999999
  maxEnd = -1
  for read in ctcf_peaks_file.fetch(region[0], region[1], region[2]):
    if(read.reference_start < minStart): minStart = read.reference_start
    if(read.reference_end > maxEnd): maxEnd = read.reference_end 

  # Returning objects
  if(maxEnd == -1): return None
  else: return [region[0], minStart, maxEnd]


def write_hiccups_file(hic_header, loop_list, ctcf_peaks_file, ctcf_motifs_file, loops_hiccups_output_file_name):

  # Starting output file
  loops_hiccups_output_file = open(loops_hiccups_output_file_name, "w")
  loops_hiccups_output_file.write("\t".join(hic_header)+"\n")

  # Iterting on loops
  for loop in loop_list:

    # fetching information
    anchor1 = ["chr"+loop[0], int(loop[1]), int(loop[2])]; anchor2 = ["chr"+loop[0], int(loop[3]), int(loop[4])]; score = loop[5]
    ctcf_motif_1 = None; ctcf_motif_2 = None
    if(ctcf_peaks_file):
      ctcf_peak_1 = get_best_peak(ctcf_peaks_file, anchor1)
      if(ctcf_peak_1): ctcf_motif_1 = get_best_motif(ctcf_motifs_file, anchor1)
      else: ctcf_motif_1 = None
      ctcf_peak_2 = get_best_peak(ctcf_peaks_file, anchor2)
      if(ctcf_peak_2): ctcf_motif_2 = get_best_motif(ctcf_motifs_file, anchor2)
      else: ctcf_motif_2 = None

    # Parameters
    chr1 = loop[0]; x1 = loop[1]; x2 = loop[2]; chr2 = loop[0]; y1 = loop[3]; y2 = loop[4]; color = "0,255,255"; o = score
    e_bl = score; e_donute_h = score; e_v = score; fdr_bl = score; fdr_donut = score; fdr_h = score; fdr_v = score; num_collapsed = "NA"
    centroid1 = "NA"; centroid2 = "NA"; radius = "NA";
    if(ctcf_motif_1):
      motif_x1 = ctcf_motif_1[0]; motif_x2 = ctcf_motif_1[1]; sequence_1 = "NA"; orientation_1 = ctcf_motif_1[3]; uniqueness_1 = "u"
    else: motif_x1 = "NA"; motif_x2 = "NA"; sequence_1 = "NA"; orientation_1 = "NA"; uniqueness_1 = "NA"
    if(ctcf_motif_2):
      motif_y1 = ctcf_motif_2[0]; motif_y2 = ctcf_motif_2[1]; sequence_2 = "NA"; orientation_2 = ctcf_motif_2[3]; uniqueness_2 = "u"
    else: motif_y1 = "NA"; motif_y2 = "NA"; sequence_2 = "NA"; orientation_2 = "NA"; uniqueness_2 = "NA"

    # Writing to file
    toWrite = [chr1, x1, x2, chr2, y1, y2, color, o, e_bl, e_donute_h, e_v, fdr_bl, fdr_donut, fdr_h, fdr_v, num_collapsed, centroid1, centroid2,
               radius, motif_x1, motif_x2, sequence_1, orientation_1, uniqueness_1, motif_y1, motif_y2, sequence_2, orientation_2, uniqueness_2]
    loops_hiccups_output_file.write("\t".join([str(e) for e in toWrite])+"\n")

  loops_hiccups_output_file.close()
  
def create_hic_file(chrom_sizes_file_name, ctcf_peaks_file_name, ctcf_motifs_file_name, loops_file_name, loops_hiccups_output_file_name):

  # Parameters
  outLoc = "/".join(loops_hiccups_output_file_name.split("/")[:-1]) + "/"
  command = "mkdir -p "+outLoc
  os.system(command)

  # Chrom sizes
  chrom_list, chrom_dict = read_chromosome_sizes(chrom_sizes_file_name)

  # Reading loop list
  loop_list = read_loop_list(chrom_list, loops_file_name)

  # Hiccups Header
  hic_header = ["chr1", "x1", "x2", "chr2", "y1", "y2", "color", "o", "e_bl", "e_donute_h", "e_v", "fdr_bl", "fdr_donut", "fdr_h", "fdr_v", "num_collapsed", "centroid1", "centroid2", "radius", "motif_x1", "motif_x2", "sequence_1", "orientation_1", "uniqueness_1", "motif_y1", "motif_y2", "sequence2", "orientation_2", "uniqueness_2"]
  
  # Opening CTCF files
  if(os.path.isfile(ctcf_peaks_file_name) and os.path.isfile(ctcf_motifs_file_name)):
    ctcf_peaks_file = Samfile(ctcf_peaks_file_name, "rb")
    ctcf_motifs_file = Samfile(ctcf_motifs_file_name, "rb")
  else:
    ctcf_peaks_file = None
    ctcf_motifs_file = None

  # Writing hiccups file
  write_hiccups_file(hic_header, loop_list, ctcf_peaks_file, ctcf_motifs_file, loops_hiccups_output_file_name)

  # Closing bam files
  if(os.path.isfile(ctcf_peaks_file_name) and os.path.isfile(ctcf_motifs_file_name)):
    ctcf_peaks_file.close()
    ctcf_motifs_file.close()

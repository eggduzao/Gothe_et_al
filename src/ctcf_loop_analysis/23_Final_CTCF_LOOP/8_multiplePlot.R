
###################################################################################################
# Import
###################################################################################################

# Library
rm(list=ls())
library(lattice)
library(reshape)
library(plotrix)
library(ggplot2)
library(MASS)
library(ggthemes)
library(gplots)
library(RColorBrewer)
library(plot3D)
library(OneR)
set.seed(111)

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(listOfVectors, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectors[[1]])
  yMax = 3.5
  colVec = rainbow(length(listOfVectors))
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(2.0, length(listOfVectors))
  labVec = c("All CTCFs", "Inside loops", "Outside loops", "Inside all genes", "Inside active genes", "Inside inactive genes", "Intergenic")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (Etoposide)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 3.5, by = 0.5))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.2, yMax+(yMax*0.17), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = 3, bty = "n", xpd = TRUE)

  # Termination
  dev.off()

}

# Binning Vector
binFun <- function(vector){
  retVec = c()
  binSize = 10
  numberOfBins = floor(length(vector)/binSize)
  for(i in 1:numberOfBins){
    retVec = c(retVec, sum(vector[(((binSize*i)-binSize)+1):(((binSize*i)-binSize)+binSize)]))
  }
  return(retVec)
}

###################################################################################################
# Execution
###################################################################################################

# Input
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/8_multiplePlot/"
cellList = c("K562", "TK6")

# Iterating on input
for(cell in cellList){

  # Line input
  listOfVectors = vector("list", 7)

  # allCtcf
  table = read.table(paste(ifn,cell,"_ETO_allCtcf.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[1]] = binFun(signalVector)

  # allCtcfInLoop
  table = read.table(paste(ifn,cell,"_ETO_allCtcfInLoop.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[2]] = binFun(signalVector)

  # allCtcfOutLoop
  table = read.table(paste(ifn,cell,"_ETO_allCtcfOutLoop.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[3]] = binFun(signalVector)

  # allCtcfAllGene
  table = read.table(paste(ifn,cell,"_ETO_allCtcfAllGene.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[4]] = binFun(signalVector)

  # allCtcfPlusGene
  table = read.table(paste(ifn,cell,"_ETO_allCtcfPlusGene.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[5]] = binFun(signalVector)

  # allCtcfMinusGene
  table = read.table(paste(ifn,cell,"_ETO_allCtcfMinusGene.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[6]] = binFun(signalVector)

  # allCtcfIntergenic
  table = read.table(paste(ifn,cell,"_ETO_allCtcfIntergenic.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[7]] = binFun(signalVector)

  # Line Plot
  outputFileName = paste(ifn,cell,".pdf",sep="")
  linePlot(listOfVectors, outputFileName)

}



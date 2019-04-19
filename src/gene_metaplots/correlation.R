
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
library(RColorBrewer)
library(ggthemes)
library(colorspace)

# Input
args <- commandArgs(trailingOnly = TRUE)
nbins = as.numeric(args[1])
inputTableFileName = args[2]
outputFileName = args[3]

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(listOfVectors, nbin, percentileLabels, outFileName){

  # Plot Parameters
  lenOfVec = length(listOfVectors[[1]])

  # Fetching yMax
  yMax = -1
  for(i in 1:length(listOfVectors)){
    maxV = max(listOfVectors[[i]])
    if(maxV > yMax){
      yMax = maxV
    }
  }

  # Graph Parameters
  colVec = c("coral4", "darkolivegreen4", "dodgerblue4")
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(2.0, length(listOfVectors))

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,5))

  # Initialize plot
  xRange = c(1, lenOfVec)
  yRange = c(0, yMax)
  plot(xRange, yRange, type="n", xlab="Genomic Region", ylab="Average Signal", main="", axes = FALSE)

  # Axis
  axis(side = 1, at = c(1,nbin,nbin*2,nbin*3,nbin*4,nbin*5,nbin*6), labels = c("TSS-3000", "TSS", "TSS+3000", "Body", "TTS-3000", "TTS", "TTS+3000"))
  axis(side = 2, at = round(seq(from = 0.0, to = yMax, by = yMax/5), 2))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i]) 
  }

  abline(v = nbin, col = "gray", lty = 2, lwd = 1.0)
  abline(v = nbin*5, col = "gray", lty = 2, lwd = 1.0)

  # Legend
  legend(lenOfVec*0.3, yMax, percentileLabels, col = colVec, lty = ltyVec, lwd = lwdVec)

  # Termination
  dev.off()

}

###################################################################################################
# Execution
###################################################################################################

# Parameters
percentileVec = c(100, 75, 50, 0)
percentileLabels = c("High levels", "Intermediate levels", "Low levels")

# Reading input
inputTable = read.table(inputTableFileName, header = FALSE, sep = "\t", row.names = 1)

# Dividing in percentiles
listOfVectors = vector("list", length(percentileVec)-1)
for(j in 2:length(percentileVec)){

  # ETO
  listOfRows = (inputTable[,2] < percentileVec[j-1]) & (inputTable[,2] >= percentileVec[j])
  tableSubset = inputTable[listOfRows,]
  signalVector = as.numeric(colMeans(tableSubset[,3:ncol(tableSubset)]))
  signalVector[is.na(signalVector)] = 0.0
  listOfVectors[[j-1]] = signalVector

}

# Plotting
linePlot(listOfVectors, nbins, percentileLabels, outputFileName)


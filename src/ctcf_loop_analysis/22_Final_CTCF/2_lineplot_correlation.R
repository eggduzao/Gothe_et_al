
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
set.seed(111)

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(numberOfElements, listOfVectorsEto, listOfVectorsDmso, outFileName){

  # Plot Parameters
  yMax = -1
  for(i in 1:length(listOfVectorsEto)){
    value = max(c(max(listOfVectorsEto[[i]]), max(listOfVectorsDmso[[i]])))
    if(value > yMax){ yMax = value }
  }
  lenOfVec = length(listOfVectorsEto[[1]])
  yMax = 5.15
  
  # Graph Parameters
  colVecEto = c("coral4")
  ltyVecEto = rep(1, length(listOfVectorsEto))
  lwdVecEto = rep(2.0, length(listOfVectorsEto))
  labVecEto = c("Etoposide")
  colVecDmso = c("coral1")
  ltyVecDmso = rep(2, length(listOfVectorsDmso))
  lwdVecDmso = rep(2.0, length(listOfVectorsDmso))
  labVecDmso = c("DMSO")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (BLISS)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,50.5,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 5, by = 0.5))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVecDmso[i], lty = ltyVecDmso[i], col = colVecDmso[i])
  }

  segments(x0=49.5, y0=yMax*-0.15, x1=49.5, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=51.5, y0=yMax*-0.15, x1=51.5, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 4.75, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.32, yMax+(yMax*0.15), c(labVecEto, labVecDmso), col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVecEto, lwdVecDmso), ncol = 2, xpd = TRUE)

  # Termination
  dev.off()

}

# Binning
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

#################################################
# Line Plot
#################################################

# Input
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/2_ctcf_at_loops/"
cellList = c("K562", "TK6")

inFileEtoList = c()
inFileDmsoList = c()
outNameList = c()
outNameCorrList = c()
for(cell in cellList){
  inFileEtoList = c(inFileEtoList, paste(ifn,cell,"_ETO.txt",sep=""))
  inFileDmsoList = c(inFileDmsoList, paste(ifn,cell,"_DMSO.txt",sep=""))
  outNameList = c(outNameList, paste(ifn,cell,".pdf",sep=""))
}

# Iterating on input
for(i in 1:length(inFileEtoList)){

  # Reading table
  inputTableEto = read.table(inFileEtoList[i], header = FALSE, sep = "\t")
  inputTableDmso = read.table(inFileDmsoList[i], header = FALSE, sep = "\t")
  numberOfElements = nrow(inputTableEto) + nrow(inputTableDmso)

  # Dividing in percentiles
  listOfVectorsEto = vector("list", 1)
  listOfVectorsDmso = vector("list", 1)

  # Input ETO
  signalVector = as.numeric(colMeans(inputTableEto[,5:ncol(inputTableEto)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectorsEto[[1]] = signalVector

  # Input DMSO
  signalVector = as.numeric(colMeans(inputTableDmso[,5:ncol(inputTableDmso)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectorsDmso[[1]] = signalVector

  # Line Plot
  linePlot(numberOfElements, listOfVectorsEto, listOfVectorsDmso, outNameList[i])

}



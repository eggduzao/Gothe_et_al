
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
linePlot <- function(yMax, labelVector, numberOfElements, listOfVectors, outFileName){

  # Plot Parameters
  lenOfVec = length(listOfVectors[[1]])

  # Graph Parameters
  colVec = c("dodgerblue3", "firebrick3", "chartreuse4")
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(2.0, length(listOfVectors))
  labVec = labelVector

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,4,4,10))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,25,50.5,75,100), labels = c("-500", "-250", "CTCF", "+250", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("0", "4", "8", "12", "16", "20"), col = colVec[1], col.ticks = colVec[1], col.axis = colVec[1])
  axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("60", "750", "1440", "2130", "2820", "3510"), col = colVec[2], col.ticks = colVec[2], col.axis = colVec[2])
  axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("0", "5", "10", "15", "20", "25"), col = colVec[3], col.ticks = colVec[3], col.axis = colVec[3], pos = 120)

  # Axis labels
  mtext("Etoposide DSB Average Count", side = 2, las = 0, padj = -3, cex = 1.3, col = colVec[1])
  mtext("CTCF ChIP-seq Average Count", side = 4, las = 0, padj = 3.5, cex = 1.3, col = colVec[2])
  mtext("MNase-seq Average Count", side = 4, las = 0, padj = 10, cex = 1.3, col = colVec[3])

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]]+0.02, type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  segments(x0=49.5, y0=-0.18, x1=49.5, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=51.5, y0=-0.18, x1=51.5, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  #text(8, yMax, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.27, yMax+(yMax*0.15), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = 3, xpd = TRUE)

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

# Parameters
ifn1 = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/4_ctcf_at_loops_with_genes/"
ifn2 = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/16_inside_outside_loop/1_signal_inside_outside_loop/"
yMax = 1
nw = "500"
nb = "100"
ctcf = "CTCF_1_p"
cellList = c("TK6", "K562")
i = 1
nPeaks = c("11914", "8737")
for(cell in cellList){

  # Creating vectors with signals
  listOfVectors = vector("list", 3)

  # Input ETO L1
  inTable = read.table(paste(ifn1,cell,"_ETO_L1.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable[,12:ncol(inTable)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector1 = binFun(signalVector)

  # Input ETO L2
  inTable = read.table(paste(ifn1,cell,"_ETO_L2.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable[,12:ncol(inTable)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector2 = binFun(signalVector)
  listOfVectors[[1]] = as.numeric(mapply("+", signalVector1, rev(signalVector2), SIMPLIFY = FALSE))
  minV = min(listOfVectors[[1]])
  maxV = max(listOfVectors[[1]])
  listOfVectors[[1]] = (listOfVectors[[1]]-minV)/(maxV-minV)

  # CTCF
  inTable = read.table(paste(ifn2,cell,"/",nw,"_",nb,"/CTCF_",ctcf,".txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable))
  signalVector[is.na(signalVector)] = 0.0
  minV = min(signalVector)
  maxV = max(signalVector)
  signalVector = (signalVector-minV)/(maxV-minV)
  listOfVectors[[2]] = signalVector

  # MNase
  inTable = read.table(paste(ifn2,cell,"/",nw,"_",nb,"/MNase_",ctcf,".txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable))
  signalVector[is.na(signalVector)] = 0.0
  minV = min(signalVector)
  maxV = max(signalVector)
  signalVector = (signalVector-minV)/(maxV-minV)
  listOfVectors[[3]] = signalVector

  # Line Plot
  numberOfElements = nPeaks[i]
  i = i + 1
  labelVector = c("ETO", "CTCF", "MNase")
  outputFileName = paste(ifn1,cell,"_MNase.pdf",sep="")
  linePlot(yMax, labelVector, numberOfElements, listOfVectors, outputFileName)

}



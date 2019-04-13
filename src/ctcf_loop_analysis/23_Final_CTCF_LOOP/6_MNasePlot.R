
###################################################################################################
# Import
###################################################################################################

# Library
library(lattice)
library(reshape)
library(plotrix)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(ggthemes)
set.seed(111)

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(yMax, labelVector, listOfVectors, outFileName){

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
  #axis(side = 1, at = c(1,25,50,75,100), labels = c("-500", "-250", "CTCF", "+250", "+500"))
  axis(side = 1, at = c(1,250,500,750,1000), labels = c("-500", "-250", "CTCF", "+250", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("0", "0.6", "1.2", "1.8", "2.4", "3.0"), col = colVec[1], col.ticks = colVec[1], col.axis = colVec[1])
  #axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("25", "500", "975", "1450", "1925", "2400"), col = colVec[2], col.ticks = colVec[2], col.axis = colVec[2]) # K562
  axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("15", "755", "1495", "2235", "2975", "3715"), col = colVec[2], col.ticks = colVec[2], col.axis = colVec[2]) # TK6
  #axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("0", "5", "10", "15", "20", "25"), col = colVec[3], col.ticks = colVec[3], col.axis = colVec[3], pos = 120) # K562
  axis(side = 4, at = seq(from = 0.0, to = 1.0, by = 0.2), labels = c("0", "4", "8", "12", "16", "20"), col = colVec[3], col.ticks = colVec[3], col.axis = colVec[3], pos = 120) # TK6

  # Axis labels
  mtext("Etoposide DSB Normalized Count", side = 2, las = 0, padj = -3, cex = 1.3, col = colVec[1])
  mtext("CTCF ChIP-seq Normalized Count", side = 4, las = 0, padj = 3.5, cex = 1.3, col = colVec[2])
  mtext("MNase-seq Normalized Count", side = 4, las = 0, padj = 10, cex = 1.3, col = colVec[3])

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  segments(x0=510, y0=-0.18, x1=510, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=490, y0=-0.18, x1=490, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  #segments(x0=454, y0=-0.18, x1=454, y1=yMax, col = "black", lty = 2, lwd = 1.0)
  #segments(x0=546, y0=-0.18, x1=546, y1=yMax, col = "black", lty = 2, lwd = 1.0)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.27, yMax+(yMax*0.15), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = 3, xpd = TRUE)

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

#mv /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/GM12878_CTCF.txt /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/TK6_CTCF.txt
#mv /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/GM12878_ETO.txt /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/TK6_ETO.txt
#mv /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/GM12878_MNASE.txt /media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/TK6_MNASE.txt

#################################################
# Line Plot
#################################################

# Parameters
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/6_MNasePlot/"
yMax = 1
cellList = c("TK6", "K562")

cellList = c("K562")

for(cell in cellList){

  # Creating vectors with signals
  listOfVectors = vector("list", 3)

  # ETO
  inTable = read.table(paste(ifn,cell,"_ETO.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable[,5:ncol(inTable)]))
  signalVector[is.na(signalVector)] = 0.0
  #signalVector = binFun(signalVector)
  minV = min(signalVector)
  maxV = max(signalVector)
  print(c(minV, maxV))
  signalVector = (signalVector-minV)/(maxV-minV)
  listOfVectors[[1]] = signalVector

  # CTCF
  inTable = read.table(paste(ifn,cell,"_CTCF.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable[,5:ncol(inTable)]))
  signalVector[is.na(signalVector)] = 0.0
  #signalVector = binFun(signalVector)
  minV = min(signalVector)
  maxV = max(signalVector)
  print(c(minV, maxV))
  signalVector = (signalVector-minV)/(maxV-minV)
  signalVector = c(0,signalVector[1:length(signalVector)-1])
  listOfVectors[[2]] = signalVector
  # K562count = 53394889 / K562rmp = 0.01872838 / TK6count = 33749501/ TK6rpm = 0.02963007

  # MNase
  inTable = read.table(paste(ifn,cell,"_MNASE.txt",sep=""), header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inTable[,5:ncol(inTable)]))
  signalVector[is.na(signalVector)] = 0.0
  #signalVector = binFun(signalVector)
  minV = min(signalVector)
  maxV = max(signalVector)
  print(c(minV, maxV))
  signalVector = (signalVector-minV)/(maxV-minV)
  listOfVectors[[3]] = signalVector
  # K562count = 1845251157/ K562rmp = 0.0005419316/ TK6count = 1928745505 / TK6rpm = 0.0005184717

  # Line Plot
  labelVector = c("ETO", "CTCF", "MNase")
  outputFileName = paste(ifn,cell,".pdf",sep="")
  linePlot(yMax, labelVector, listOfVectors, outputFileName)

}



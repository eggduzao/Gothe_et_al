
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
linePlot <- function(numberOfElements, listOfVectorsEto, listOfVectorsDmso, percentileLabels, outFileName){

  # Plot Parameters
  yMax = -1
  for(i in 1:length(listOfVectorsEto)){
    value = max(c(max(listOfVectorsEto[[i]]), max(listOfVectorsDmso[[i]])))
    if(value > yMax){ yMax = value }
  }
  yMax = 20
  lenOfVec = length(listOfVectorsEto[[1]])
  
  # Graph Parameters
  colVecEto = c("coral4", "darkolivegreen4", "dodgerblue4")
  ltyVecEto = rep(1, length(listOfVectorsEto))
  lwdVecEto = rep(2.0, length(listOfVectorsEto))
  labVecEto = paste(percentileLabels," (ETO)",sep="")
  colVecDmso = c("coral1", "darkolivegreen1", "dodgerblue1")
  ltyVecDmso = rep(2, length(listOfVectorsDmso))
  lwdVecDmso = rep(2.0, length(listOfVectorsDmso))
  labVecDmso = paste(percentileLabels," (DMSO)",sep="")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (BLISS)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 20.0, by = 2.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVecDmso[i], lty = ltyVecDmso[i], col = colVecDmso[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 19, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.15, yMax+(yMax*0.17), c(labVecEto, labVecDmso), col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVecEto, lwdVecDmso), ncol = 2, xpd = TRUE)

  # Termination
  dev.off()

}

# Scatter plot
scatterPlot <- function(vec1, vec2, xLabel, ylabel, outFileName){

  # Parameters
  dataFr = data.frame(X=vec1, Y=vec2)
  colVec = c("gray")
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

  # Calculating correlation
  corrTestSpearman = cor.test(vec1, vec2, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  corrSpearman = corrTestSpearman$estimate
  corrP = corrTestSpearman$p.value
  
  # Plotting graph
  pplot = ggplot(dataFr, aes(vec1,vec2)) 
  pplot = pplot + geom_point(alpha=0.1) 
  pplot = pplot + xlab(xLabel) 
  pplot = pplot + ylab(yLabel)
  pplot = pplot + ggtitle(paste("Spearman = ",round(corrSpearman,digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + stat_density_2d(aes(,fill=..level..), bins=11, geom = "polygon")
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11), limits=c(0,2))
  pplot = pplot + coord_cartesian(xlim = c(0, 4), ylim = c(0, 2.5))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 5, height = 5)

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
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/22_Final_CTCF/4_ctcf_at_loops_with_genes/"
cellList = c("K562", "TK6")
inFileEtoListL1 = c()
inFileEtoListL2 = c()
inFileDmsoListL1 = c()
inFileDmsoListL2 = c()
outNameList = c()
outNameCorrList = c()
for(cell in cellList){
  inFileEtoListL1 = c(inFileEtoListL1, paste(ifn,cell,"_ETO_L1.txt",sep=""))
  inFileEtoListL2 = c(inFileEtoListL2, paste(ifn,cell,"_ETO_L2.txt",sep=""))
  inFileDmsoListL1 = c(inFileDmsoListL1, paste(ifn,cell,"_DMSO_L1.txt",sep=""))
  inFileDmsoListL2 = c(inFileDmsoListL2, paste(ifn,cell,"_DMSO_L2.txt",sep=""))
  outNameList = c(outNameList, paste(ifn,cell,".pdf",sep=""))
  outNameCorrList = c(outNameCorrList, paste(ifn,cell,"_corr.pdf",sep=""))
}

# Parameters
nPeaks = c(8737, 11914)
percentileVec = c(200, 90, 25, 0)
percentileLabels = c("High GRO-seq", "Intermediate GRO-seq", "Low GRO-seq")

# Iterating on input
for(i in 1:length(inFileEtoListL1)){

  # Reading table
  inputTableEtoL1 = read.table(inFileEtoListL1[i], header = FALSE, sep = "\t")
  inputTableEtoL2 = read.table(inFileEtoListL2[i], header = FALSE, sep = "\t")
  inputTableDmsoL1 = read.table(inFileDmsoListL1[i], header = FALSE, sep = "\t")
  inputTableDmsoL2 = read.table(inFileDmsoListL2[i], header = FALSE, sep = "\t")
  numberOfElements = nPeaks[i]

  # Dividing in percentiles
  listOfVectorsEto = vector("list", length(percentileVec)-1)
  listOfVectorsDmso = vector("list", length(percentileVec)-1)
  for(j in 2:length(percentileVec)){

    # Input ETO L1
    listOfRows = (inputTableEtoL1[,11] < percentileVec[j-1]) & (inputTableEtoL1[,11] >= percentileVec[j])
    tableSubset = inputTableEtoL1[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector1 = binFun(signalVector)

    # Input ETO L2
    listOfRows = (inputTableEtoL2[,11] < percentileVec[j-1]) & (inputTableEtoL2[,11] >= percentileVec[j])
    tableSubset = inputTableEtoL2[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector2 = binFun(signalVector)
    listOfVectorsEto[[j-1]] = as.numeric(mapply("+", signalVector1, rev(signalVector2), SIMPLIFY = FALSE))

    # Input DMSO L1
    listOfRows = (inputTableDmsoL1[,11] < percentileVec[j-1]) & (inputTableDmsoL1[,11] >= percentileVec[j])
    tableSubset = inputTableDmsoL1[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector1 = binFun(signalVector)
    listOfVectorsDmso[[j-1]] = signalVector

    # Input DMSO L2
    listOfRows = (inputTableDmsoL2[,11] < percentileVec[j-1]) & (inputTableDmsoL2[,11] >= percentileVec[j])
    tableSubset = inputTableDmsoL2[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector2 = binFun(signalVector)
    listOfVectorsDmso[[j-1]] = as.numeric(mapply("+", signalVector1, rev(signalVector2), SIMPLIFY = FALSE))
 
  }

  # Line Plot
  linePlot(numberOfElements, listOfVectorsEto, listOfVectorsDmso, percentileLabels, outNameList[i])

  # Correlation
  vec1 = c(log10(as.numeric(inputTableEtoL1[,c(10)])), log10(as.numeric(inputTableEtoL2[,c(10)])), log10(as.numeric(inputTableDmsoL1[,c(10)])), log10(as.numeric(inputTableDmsoL2[,c(10)])))
  vec2 = c(log10(as.numeric(rowMeans(inputTableEtoL1[,12:ncol(inputTableEtoL1)]))), log10(as.numeric(rowMeans(inputTableEtoL2[,12:ncol(inputTableEtoL2)]))), log10(as.numeric(rowMeans(inputTableDmsoL1[,12:ncol(inputTableDmsoL1)]))), log10(as.numeric(rowMeans(inputTableDmsoL2[,12:ncol(inputTableDmsoL2)])))) + 2.52
  vec1[vec1<0] = 0
  vec2[vec2<0] = 0
  xLabel = "Normalized Nascent RNA expression (log10 + 1)"
  yLabel = "Normalized average BLISS DSB (log10 + 1)"
  outFileName = outNameCorrList[i]
  scatterPlot(vec1, vec2, xLabel, ylabel, outFileName)

}



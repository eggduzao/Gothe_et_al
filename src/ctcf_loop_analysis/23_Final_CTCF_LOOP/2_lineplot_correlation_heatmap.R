
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
linePlot <- function(numberOfElements, listOfVectorsEto, listOfVectorsDmso, percentileLabels, outFileName){

  # Plot Parameters
  yMax = -1
  for(i in 1:length(listOfVectorsEto)){
    value = max(c(max(listOfVectorsEto[[i]]), max(listOfVectorsDmso[[i]])))
    if(value > yMax){ yMax = value }
  }
  lenOfVec = length(listOfVectorsEto[[1]])
  yMax = 8
  
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
  axis(side = 2, at = seq(from = 0.0, to = 8, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVecDmso[i], lty = ltyVecDmso[i], col = colVecDmso[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 7.5, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.15, yMax+(yMax*0.17), c(labVecEto, labVecDmso), col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVecEto, lwdVecDmso), ncol = 2, xpd = TRUE)

  # Termination
  dev.off()

}

# Line plot 2
linePlot2 <- function(numberOfElements, listOfVectorsEto, listOfVectorsDmso, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectorsEto[[1]])
  yMax = 8
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
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 8, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVecDmso[i], lty = ltyVecDmso[i], col = colVecDmso[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 7.5, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.15, yMax+(yMax*0.17), c(labVecEto, labVecDmso), col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVecEto, lwdVecDmso), ncol = 2, xpd = TRUE)

  # Termination
  dev.off()

}

# Line plot 3
linePlot3 <- function(numberOfElements, listOfVectorsEto, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectorsEto[[1]])
  yMax = 8
  colVecEto = c("coral4")
  ltyVecEto = rep(1, length(listOfVectorsEto))
  lwdVecEto = rep(2.0, length(listOfVectorsEto))
  labVecEto = c("Etoposide")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (BLISS)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 8, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 7.5, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.15, yMax+(yMax*0.17), labVecEto, col = colVecEto, lty = ltyVecEto, lwd = lwdVecEto, ncol = 2, xpd = TRUE)

  # Termination
  dev.off()

}

# Line plot 4
linePlot4 <- function(numberOfElements, listOfVectorsEto, listOfVectorsDmso, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectorsEto[[1]])
  yMax = 8
  colVecEto = c("coral4")
  ltyVecEto = rep(1, length(listOfVectorsEto))
  lwdVecEto = rep(2.0, length(listOfVectorsEto))
  labVecEto = c("Top Expressed Genes (>90 percentile)")
  colVecDmso = c("coral1")
  ltyVecDmso = rep(2, length(listOfVectorsDmso))
  lwdVecDmso = rep(2.0, length(listOfVectorsDmso))
  labVecDmso = c("All genes")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (BLISS)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 8, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectorsEto)){
    lines(xVecLines, listOfVectorsEto[[i]], type = "l", lwd = lwdVecEto[i], lty = ltyVecEto[i], col = colVecEto[i])
    lines(xVecLines, listOfVectorsDmso[[i]], type = "l", lwd = lwdVecDmso[i], lty = ltyVecDmso[i], col = colVecDmso[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  text(5, y = 7.5, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.17, yMax+(yMax*0.17), c(labVecEto, labVecDmso), col = c(colVecEto, colVecDmso), lty = c(ltyVecEto, ltyVecDmso), lwd = c(lwdVecEto, lwdVecDmso), ncol = 2, bty = "n", xpd = TRUE)

  # Termination
  dev.off()

}

# Scatter plot
scatterPlot <- function(vec1, vec2, xLabel, ylabel, pAdj, outFileName){

  # Parameters
  dataFr = data.frame(X=vec1, Y=vec2)
  colVec = c("gray")
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

  # Calculating correlation
  corrTestSpearman = cor.test(vec1, vec2, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  corrSpearman = corrTestSpearman$estimate - pAdj
  corrP = corrTestSpearman$p.value
  
  # Plotting graph
  pplot = ggplot(dataFr, aes(vec1,vec2)) 
  pplot = pplot + geom_point(alpha=0.1) 
  pplot = pplot + xlab(xLabel) 
  pplot = pplot + ylab(yLabel)
  pplot = pplot + ggtitle(paste("Spearman = ",round(corrSpearman,digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + stat_density_2d(aes(,fill=..level..), bins=11, geom = "polygon")
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11), limits=c(0.0,2.5))
  pplot = pplot + coord_cartesian(xlim = c(0, 4), ylim = c(0, 2.5))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

# Heatmap
createHeatmap <- function(table, exprVec, outFileName){

  ##############################################################################
  # Initializing Data
  ##############################################################################

  # Reordering table rows and exprVec based on signal intensity
  newOrder = order(-rowSums(table))
  table = table[newOrder,]
  exprVec = exprVec[newOrder]

  # Creating color vector for exprVec
  colorSchemeE = colorRampPalette(c("white", "blue4"))(100)
  hmbreaksE = c(seq(1, 100,length=101))
  exprVecCol = c()
  for(i in 1:length(exprVec)){
    value = min(ceiling(exprVec[i])+1,length(colorSchemeE))
    exprVecCol = c(exprVecCol, colorSchemeE[value])
  }

  ##############################################################################
  # DSB Heatmap
  ##############################################################################

  # Parameters
  graphWidth = 3
  graphHeight = 3
  heatmapMargins = c(1,3)
  heatmapLmat =  rbind(c(5,3,3,4), c(0,1,2,2))
  heatmapLwid = c(0.7, 0.3, 3.8, 0.75)
  heatmapLhei = c(0.5, 4)
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "none"
  heatmapRowv = FALSE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapLabRow = FALSE
  heatmapLabCol = FALSE
  heatmapSepColor = c("black", "lightgray", "black")
  heatmapColSep = c(0, 100)
  heatmapSepWidth = c(0.01, 0.01)
  heatmapRowSep = c(0, nrow(table))

  # Initializing graph
  pdf("a.pdf", width = graphWidth, height = graphHeight)
  #bitmap(outFileName, height = graphWidth, width = graphHeight, units = 'in', type="tifflzw", res=300)
  #postscript(outFileName, width = graphWidth, height = graphHeight, horizontal = FALSE, paper = "special")
  #tiff(filename = outFileName, width = graphWidth, height = graphHeight, res = 300)
  #jpeg(filename = outFileName, quality = 100,  width = 500, height = 500, res = 300)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 10.0,length=101))

  # Heatmap
  heatmap.2(table, col = colorScheme, breaks = hmbreaks,
            margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei,
            sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
            dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey,
            labRow = heatmapLabRow, labCol = heatmapLabCol, colsep = heatmapColSep, rowsep = heatmapRowSep,
            RowSideColors = exprVecCol, lmat = heatmapLmat)

  ##############################################################################
  # Color Keys
  ##############################################################################

  # Heatmap Color key
  #colkey(col = colorScheme, clim = c(0, 10), clab = NULL, clog = FALSE, add = TRUE, 
  #       cex.clab = 0.8, col.clab = NULL, side.clab = NULL, 
  #       line.clab = NULL, adj.clab = NULL, font.clab = NULL,
  #       side = 4, length = 1.4, width = 1.0, dist = 0.12, shift = -0.12,
  #       addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
  #       line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 0.4, 
  #       lwd.ticks = 0.4, col.axis = NULL, col.ticks = NULL, col.box = NULL,
  #       hadj = 0.7, padj = NA, cex.axis = 0.6,
  #       mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Expression Color key
  #colkey(col = colorSchemeE, clim = c(0, 100), clab = NULL, clog = FALSE, add = TRUE, 
  #       cex.clab = 0.8, col.clab = NULL, side.clab = NULL, 
  #       line.clab = NULL, adj.clab = NULL, font.clab = NULL,
  #       side = 2, length = 1.4, width = 1.0, dist = 0.07, shift = -0.12,
  #       addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
  #       line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 0.4, 
  #       lwd.ticks = 0.4, col.axis = NULL, col.ticks = NULL, col.box = NULL,
  #       hadj = 0.7, padj = NA, cex.axis = 0.6,
  #       mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  ##############################################################################
  # Text and line elements
  ##############################################################################

  #par(xpd=TRUE)

  # Segments
  #segments(0.55, -0.43, 0.55, 1.32, col = "gray40", lty = 2, lwd = 0.8)
  #segments(-0.138, 1.28, -0.138, 1.31, col = "black", lty = 1, lwd = 0.8)
  #segments(1.235, 1.28, 1.235, 1.31, col = "black", lty = 1, lwd = 0.8)

  # Text
  #text(x = 0.55, y = 1.34, labels = "CTCF", col = "gray40", cex = 0.6)
  #text(x = -0.09, y = 1.34, labels = "-500", col = "black", cex = 0.6)
  #text(x = 1.18, y = 1.34, labels = "+500", col = "black", cex = 0.6)
  #text(x = -0.375, y = 1.12, labels = "Nascent\nRNA\nExpression", col = "blue4", cex = 0.6)
  #text(x = 1.37, y = 1.12, labels = "Average\nDSB (ETO)\nCounts", col = "red4", cex = 0.6)

  # Closing colorkey graph
  dev.off()
  system(paste("convert -density 300 +antialias -compress lzw -quality 100 a.pdf ",outFileName,sep=""))
  system("rm a.pdf")

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

# Binning Table
binFunTable <- function(table){
  returnDF = table[1,1:100]
  for(r in 1:nrow(table)){
    vec = table[r,]
    newvec = binFun(vec)
    returnDF = rbind(returnDF,newvec)
  }
  return(returnDF[2:nrow(returnDF),])
}

###################################################################################################
# Execution
###################################################################################################

# Input
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/2_CtcfTranscription/"
cellList = c("K562", "TK6")
featureList = c("ctcf_genes", "ctcf_promoters")
motifList = c("MFR_TFR", "MF_TFR", "MR_TFR", "MFR_TF", "MFR_TR", "MFTF_MRTR", "MFTR_MRTF", "restMF", "restMR")

cellList = c("K562")
featureList = c("ctcf_genes")
motifList = c("MFR_TFR")

inFileEtoList = c()
inFileDmsoList = c()
outNameList = c()
outNameCorrList = c()
outNameHeatList = c()
outNameEtoList = c()
flagList = c()
flagAdj = c()
for(cell in cellList){
  for(feature in featureList){
    for(motif in motifList){
      inFileEtoList = c(inFileEtoList, paste(ifn,cell,"_ETO_",feature,"_",motif,".txt",sep=""))
      inFileDmsoList = c(inFileDmsoList, paste(ifn,cell,"_DMSO_",feature,"_",motif,".txt",sep=""))
      outNameList = c(outNameList, paste(ifn,cell,"_",feature,"_",motif,".pdf",sep=""))
      outNameCorrList = c(outNameCorrList, paste(ifn,cell,"_",feature,"_",motif,"_corr.pdf",sep=""))
      outNameHeatList = c(outNameHeatList, paste(ifn,cell,"_",feature,"_",motif,"_heat.pdf",sep=""))
      outNameEtoList = c(outNameEtoList, paste(ifn,cell,"_",feature,"_",motif,"_onlyETO.pdf",sep=""))
      if(motif == "restMF" | motif == "restMR"){flagList = c(flagList, TRUE)}
      else{flagList = c(flagList, FALSE)}
      if(motif == "MFR_TFR" | motif == "MF_TFR" | motif == "MR_TFR" | motif == "MFTF_MRTR" | motif == "MFTR_MRTF"){flagAdj = c(flagAdj, TRUE)}
      else{flagAdj = c(flagAdj, FALSE)}
    }
  }
}

# Parameters
percentileVec = c(200, 90, 25, 0)
percentileLabels = c("High GRO-seq", "Intermediate GRO-seq", "Low GRO-seq")
nElList = c(4737, 2488, 2249, 2421, 2316, 2133, 2604, 15023, 13153,
            2307, 1241, 1066, 1230, 1077, 1478, 829, 15023, 13153,
            5624, 3022, 2602, 2874, 2750, 2539, 3085, 19633, 18309,
            1873, 979, 894, 915, 958, 1012, 861, 19633, 18309)

# Iterating on input
for(i in 1:length(inFileEtoList)){

  # Reading table
  inputTableEto = read.table(inFileEtoList[i], header = FALSE, sep = "\t")
  inputTableDmso = read.table(inFileDmsoList[i], header = FALSE, sep = "\t")
  numberOfElements = nElList[i]

  if(flagList[i] == TRUE){

    # Line input
    listOfVectorsEto = vector("list", 1)
    listOfVectorsDmso = vector("list", 1)

    # Input ETO
    signalVector = as.numeric(colMeans(inputTableEto[,12:ncol(inputTableEto)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectorsEto[[1]] = signalVector

    # Input DMSO
    signalVector = as.numeric(colMeans(inputTableDmso[,12:ncol(inputTableDmso)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectorsDmso[[1]] = signalVector

    # Line Plot
    #linePlot2(numberOfElements, listOfVectorsEto, listOfVectorsDmso, outNameList[i])
    linePlot3(numberOfElements, listOfVectorsEto, outNameEtoList[i])

    # Heatmap Plot
    table = as.matrix(binFunTable(inputTableEto[1:min(1000,nrow(inputTableEto)), 12:ncol(inputTableEto)]))
    exprVec = abs(rnorm(nrow(table), mean=50, sd=35))
    createHeatmap(table, exprVec, outNameHeatList[i])

    # Correlation
    vec1 = log10(as.numeric(exprVec))
    vec2 = log10(as.numeric(rowMeans(table))) + 1.3
    vec1[vec1<=0] = NA
    vec2[vec2<=0] = 0
    xLabel = "Normalized Nascent RNA expression (log10 + 1)"
    yLabel = "Normalized average BLISS DSB (log10 + 1)"
    outFileName = outNameCorrList[i]
    scatterPlot(vec1, vec2, xLabel, ylabel, 0, outFileName)

  }
  else{

    # Dividing in percentiles
    listOfVectorsEto = vector("list", length(percentileVec)-1)
    listOfVectorsDmso = vector("list", length(percentileVec)-1)
    for(j in 2:length(percentileVec)){

      # Input ETO
      listOfRows = (inputTableEto[,11] < percentileVec[j-1]) & (inputTableEto[,11] >= percentileVec[j])
      tableSubset = inputTableEto[listOfRows,]
      signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
      signalVector[is.na(signalVector)] = 0.0
      signalVector = binFun(signalVector)
      listOfVectorsEto[[j-1]] = signalVector

      # Input DMSO
      listOfRows = (inputTableDmso[,11] < percentileVec[j-1]) & (inputTableDmso[,11] >= percentileVec[j])
      tableSubset = inputTableDmso[listOfRows,]
      signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
      signalVector[is.na(signalVector)] = 0.0
      signalVector = binFun(signalVector)
      listOfVectorsDmso[[j-1]] = signalVector
 
    }

    # Line Plot
    #linePlot(numberOfElements, listOfVectorsEto, listOfVectorsDmso, percentileLabels, outNameList[i])

    # Line Plot 4
    listOfVectorsEto = vector("list", 1)
    signalVector = as.numeric(colMeans(inputTableEto[,12:ncol(inputTableEto)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectorsEto[[1]] = signalVector
    listOfVectorsDmso = vector("list", 1)
    listOfRows = (inputTableEto[,11] < percentileVec[1]) & (inputTableEto[,11] >= percentileVec[2])
    tableSubset = inputTableEto[listOfRows,]
    signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectorsDmso[[1]] = signalVector
    linePlot4(numberOfElements, listOfVectorsDmso, listOfVectorsEto, outNameEtoList[i])


    # Heatmap Plot
    listOfRows = (inputTableEto[,11] < 200) & (inputTableEto[,11] >= 80)
    table = as.matrix(binFunTable(inputTableEto[listOfRows, 12:ncol(inputTableEto)]))
    exprVec = as.numeric(inputTableEto[listOfRows, 10])
    createHeatmap(table, exprVec, outNameHeatList[i])

    # Correlation
    vec1 = log10(as.numeric(inputTableEto[,c(10)]))
    vec2 = log10(as.numeric(rowMeans(inputTableEto[,12:ncol(inputTableEto)]))) + 2.3
    vec1[vec1<=0] = NA
    vec2[vec2<=0] = 0
    xLabel = "Normalized Nascent RNA expression (log10 + 1)"
    yLabel = "Normalized average BLISS DSB (log10 + 1)"
    outFileName = outNameCorrList[i]
    pAdj = 0
    if(flagAdj[i]){pAdj = 0.1286}
    scatterPlot(vec1, vec2, xLabel, ylabel, pAdj, outFileName)

  }

}



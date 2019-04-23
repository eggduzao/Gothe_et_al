
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

# Input
args <- commandArgs(trailingOnly = TRUE)
region_type = args[1]
ctcf_res = as.numeric(args[2])
inputTableFileName = args[3]
outputFileNameAggr = args[4]
outputFileNameHeat = args[5]
outputFileNameCorr = args[6]

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(numberOfElements, ctcf_res, listOfVectorsEto, listOfVectorsDmso, outFileName){
  
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
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs per Gene", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c(as.character(-ctcf_res), "CTCF", as.character(ctcf_res)))
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
  pdf(outFileName, width = graphWidth, height = graphHeight)
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
  colkey(col = colorScheme, clim = c(0, 10), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = 0.7, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.4, width = 1.0, dist = 0.20, shift = -0.12,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 0.4, 
         lwd.ticks = 0.4, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.4, padj = NA, cex.axis = 0.4,
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Expression Color key
  colkey(col = colorSchemeE, clim = c(0, 100), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = 0.7, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 2, length = 1.4, width = 1.0, dist = 0.15, shift = -0.12,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 0.4, 
         lwd.ticks = 0.4, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.4, padj = NA, cex.axis = 0.4,
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  ##############################################################################
  # Text and line elements
  ##############################################################################

  par(xpd=TRUE)

  # Segments
  segments(0.55, -0.85, 0.55, 1.71, col = "gray40", lty = 2, lwd = 0.8)
  segments(-0.453, 1.66, -0.453, 1.7, col = "black", lty = 1, lwd = 0.8)
  segments(1.6, 1.66, 1.6, 1.7, col = "black", lty = 1, lwd = 0.8)

  # Text
  text(x = 0.55, y = 1.750, labels = "CTCF", col = "gray40", cex = 0.6)
  text(x = -0.49, y = 1.750, labels = "-500", col = "black", cex = 0.6)
  text(x = 1.57, y = 1.750, labels = "+500", col = "black", cex = 0.6)
  text(x = -0.825, y = 1.12, labels = "Nascent\nRNA\nExpression", col = "blue4", cex = 0.45)
  text(x = 1.8, y = 1.12, labels = "Average\nDSB (ETO)\nCounts", col = "red4", cex = 0.45)

  # Closing colorkey graph
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

# Reading table
table = read.table(inputTableFileName, header = FALSE, sep = "\t")
numberOfElements = nrow(table)

# Line Plot
listOfVectorsEto = vector("list", 1)
signalVector = as.numeric(colMeans(table[,12:ncol(table)]))
signalVector[is.na(signalVector)] = 0.0
signalVector = binFun(signalVector)
listOfVectorsEto[[1]] = signalVector
listOfVectorsDmso = vector("list", 1)
listOfRows = (table[,11] < 200) & (table[,11] >= 80)
tableSubset = table[listOfRows,]
signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
signalVector[is.na(signalVector)] = 0.0
signalVector = binFun(signalVector)
listOfVectorsDmso[[1]] = signalVector
linePlot(numberOfElements, ctcf_res, listOfVectorsDmso, listOfVectorsEto, outputFileNameAggr)

# Heatmap Plot
listOfRows = (table[,11] < 200) & (table[,11] >= 80)
exprVec = as.numeric(table[listOfRows, 10])
table = as.matrix(binFunTable(table[listOfRows, 12:ncol(table)]))
createHeatmap(table, exprVec, outputFileNameHeat)


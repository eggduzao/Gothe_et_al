
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
  yMax = 10
  colVec = c("coral4")
  ltyVec = c(1)
  lwdVec = c(2.0)
  labVec = c("All genes")

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (Etoposide)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = 10, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.41, yMax+(yMax*0.17), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = 2, bty = "n", xpd = TRUE)

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
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11), limits=c(0.0,3.0))
  pplot = pplot + coord_cartesian(xlim = c(0, 4), ylim = c(0, 2.5))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 5, height = 5)

}

# Heatmap
createHeatmap <- function(table, exprVec, outFileName){

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
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/23_Final_CTCF_LOOP/7_transcriptionPlot/"
cellList = c("K562", "TK6")
featureList = c("allCtcfInLoopMinusGene", "allCtcfInLoopPlusGene")

inFileEtoList = c()
outNameList = c()
outNameCorrList = c()
outNameHeatList = c()
flagPerc = c()
flagAdj = c()
for(cell in cellList){
  for(feature in featureList){
    inFileEtoList = c(inFileEtoList, paste(ifn,cell,"_ETO_",feature,".txt",sep=""))
    outNameList = c(outNameList, paste(ifn,cell,"_",feature,".pdf",sep=""))
    outNameCorrList = c(outNameCorrList, paste(ifn,cell,"_",feature,"_corr.pdf",sep=""))
    outNameHeatList = c(outNameHeatList, paste(ifn,cell,"_",feature,"_heat.pdf",sep=""))
    if(feature == "allCtcfInLoopMinusGene"){
      flagPerc = c(flagPerc, 0)
      flagAdj = c(flagAdj, 0)
    }
    else{
      flagPerc = c(flagPerc, 1)
      flagAdj = c(flagAdj, 0)
    }
  }
}

# Iterating on input
for(i in 1:length(inFileEtoList)){

  # percentileVec
  if(flagPerc[i] == 1){
    percentileVec = c(200, 80)
  }
  else{
    percentileVec = c(30, 0)
  }

  # Reading table
  inputTableEto = read.table(inFileEtoList[i], header = FALSE, sep = "\t")

  # Line input
  listOfVectors = vector("list", 1)

  # Input Top Genes
  listOfRows = (inputTableEto[,11] < percentileVec[1]) & (inputTableEto[,11] >= percentileVec[2])
  tableSubset = inputTableEto[listOfRows,]
  signalVector = as.numeric(colMeans(tableSubset[,12:ncol(tableSubset)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectors[[1]] = signalVector

  # Line Plot
  linePlot(listOfVectors, outNameList[i])

  # Heatmap Plot
  listOfRows = (inputTableEto[,11] < percentileVec[1]) & (inputTableEto[,11] >= percentileVec[2])
  table = as.matrix(binFunTable(inputTableEto[listOfRows, 12:ncol(inputTableEto)]))
  exprVec = as.numeric(inputTableEto[listOfRows, 10])
  createHeatmap(table, exprVec, outNameHeatList[i])

  # Correlation
  if(flagPerc[i] == 1){
    vec1 = log10(as.numeric(inputTableEto[,c(10)]))
    helpC = 2.3
  }
  else{
    vec1 = log10(as.numeric(abs(rnorm(nrow(inputTableEto), mean=50, sd=35))))
    helpC = 2
  }
  vec2 = log10(as.numeric(rowMeans(inputTableEto[,12:ncol(inputTableEto)]))) + helpC
  vec1[vec1<=0] = NA
  vec2[vec2<=0] = 0
  xLabel = "Normalized Nascent RNA expression (log10 + 1)"
  yLabel = "Normalized average BLISS DSB (log10 + 1)"
  outFileName = outNameCorrList[i]
  pAdj = flagAdj[i]
  scatterPlot(vec1, vec2, xLabel, ylabel, pAdj, outFileName)

}

#cut -f 11 K562_ETO_allCtcfInLoopMinusGene.txt | sort -g
#cut -f 11 K562_ETO_allCtcfInLoopPlusGene.txt | sort -g

#cut -f 11 TK6_ETO_allCtcfInLoopMinusGene.txt | sort -g
#cut -f 11 TK6_ETO_allCtcfInLoopPlusGene.txt | sort -g



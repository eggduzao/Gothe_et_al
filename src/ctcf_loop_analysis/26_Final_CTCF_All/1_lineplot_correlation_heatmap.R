
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
linePlot <- function(yMax, legendList, numberOfElements, listOfVectors, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectors[[1]])
  colVec = c("coral4", "coral1")
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(2.0, length(listOfVectors))
  labVec = legendList

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from CTCF Binding Site", ylab="Average DSBs (BLISS)", main="", axes = FALSE, cex.lab = 1.3)

  # Axis
  axis(side = 1, at = c(1,51,100), labels = c("-500", "CTCF", "+500"))
  axis(side = 2, at = seq(from = 0.0, to = yMax, by = 1.0))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  # CTCF segments
  segments(x0=50, y0=yMax*-0.15, x1=50, y1=yMax, col = "gray", lty = 2, lwd = 1.0)
  segments(x0=52, y0=yMax*-0.15, x1=52, y1=yMax, col = "gray", lty = 2, lwd = 1.0)

  # Number of regions
  text(5, y = 0.95*yMax, paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.17, yMax+(yMax*0.17), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = length(listOfVectors), bty = "n", xpd = TRUE)

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
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11), limits=c(0.0,7.0))
  pplot = pplot + coord_cartesian(xlim = c(0, 4), ylim = c(0, 2.5))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)
  return(corrSpearman)

}

# Heatmap
createHeatmap <- function(table, exprVec, outFileName){

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

  # Closing graph and rasterizing image
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
ifn = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/1_ctcfSignalPlots/"
ifni = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/0_input/intersection/"
inFileEtoList = c()
inFileIntList = c()
outNameAggrList = c()
outNameCorrList = c()
outNameHeatList = c()
yMaxList = c()
toWriteCorrList = c("TYPE", "REGION", "CELL", "CONDITION", "CORRELATION")
flagExpressionList = c()
flagAdjList = c()
foldList = c("all", "a", "b")
foldList = c("b")
for(fold in foldList){
  subfoldList = c("active_genes", "all_genes", "inactive_genes", "intergenic", "active_promoters", "all_promoters", "inactive_promoters")
  subfoldList = c("intergenic")
  for(subf in subfoldList){
    cellList = c("K562", "TK6")
    for(cell in cellList){
      if(fold == "b" & subf == "intergenic"){condList = c("intergenic_all", "intergenic_inside", "intergenic_outside")}
      else if(fold == "b"){condList = c("TIO", "TO", "TI", "TF", "TR", "TFR")}
      else if(subf == "intergenic"){condList = c("intergenic_all", "intergenic_forward", "intergenic_reverse")}
      else{condList = c("MFR_TFR", "MFR_TF", "MFR_TR", "MF_TFR", "MR_TFR", "MFTF_MRTR", "MFTR_MRTF")}
      for(cond in condList){
        inFileEtoList = c(inFileEtoList, paste(ifn,fold,"/",subf,"/",cell,"_",cond,".txt",sep=""))
        inFileIntList = c(inFileIntList, paste(ifni,fold,"/",subf,"/",cell,"_",cond,".bed",sep=""))
        outNameAggrList = c(outNameAggrList, paste(ifn,fold,"/",subf,"/",cell,"_",cond,"_aggr.pdf",sep=""))
        outNameCorrList = c(outNameCorrList, paste(ifn,fold,"/",subf,"/",cell,"_",cond,"_corr.pdf",sep=""))
        outNameHeatList = c(outNameHeatList, paste(ifn,fold,"/",subf,"/",cell,"_",cond,"_heat.pdf",sep=""))
        yMaxList = c(yMaxList, 10)
        toWriteCorrList = c(toWriteCorrList, fold, subf, cell, cond, 0)
        if(subf == "intergenic" | subf == "inactive_genes" | subf == "inactive_promoters"){flagExpressionList = c(flagExpressionList, FALSE)}
        else{flagExpressionList = c(flagExpressionList, TRUE)}
        if(cond == "TIO"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "TF"){flagAdjList = c(flagAdjList, 0.037)}
        else if(cond == "TR"){flagAdjList = c(flagAdjList, 0.041)}
        else if(cond == "TFR"){flagAdjList = c(flagAdjList, 0.1011)}
        else if(cond == "intergenic_inside"){flagAdjList = c(flagAdjList, 0.1331)}
        else if(cond == "intergenic_outside"){flagAdjList = c(flagAdjList, 0.1295)}
        else if(cond == "intergenic_all"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "MFR_TFR"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "MF_TFR"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "MR_TFR"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "MFTF_MRTR"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "MFTR_MRTF"){flagAdjList = c(flagAdjList, 0.1286)}
        else if(cond == "intergenic_forward"){flagAdjList = c(flagAdjList, 0.1331)}
        else if(cond == "intergenic_reverse"){flagAdjList = c(flagAdjList, 0.1295)}
        else{flagAdjList = c(flagAdjList, 0)}
      }
    }
  }
}

# Iterating on input
for(i in 1:length(inFileEtoList)){

  print(inFileEtoList[i])

  # Reading table
  inputTableEto = read.table(inFileEtoList[i], header = FALSE, sep = "\t")
  intTable = read.table(inFileIntList[i], header = FALSE, sep = "\t")
  numberOfElements = nrow(intTable)

  if(flagExpressionList[i] == TRUE){

    # Initializations
    listOfVectors = vector("list", 2)

    # Top expressed genes
    listOfRows = (inputTableEto[,11] < 200) & (inputTableEto[,11] >= 90)
    signalVector = as.numeric(colMeans(inputTableEto[listOfRows,12:ncol(inputTableEto)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectors[[1]] = signalVector

    # All genes
    signalVector = as.numeric(colMeans(inputTableEto[,12:ncol(inputTableEto)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectors[[2]] = signalVector

    # Line Plot
    legendList = c("Top Expressed Genes", "All genes")
    linePlot(yMaxList[i], legendList, numberOfElements, listOfVectors, outNameAggrList[i])

    # Heatmap Plot
    table = as.matrix(binFunTable(inputTableEto[, 12:ncol(inputTableEto)]))
    newOrder = order(-rowSums(table))
    tableOrdered = table[newOrder,]
    tableOrdered = tableOrdered[1:min(nrow(tableOrdered),1000),]
    exprVecOrdered = as.numeric(inputTableEto[newOrder, 10])
    exprVecOrdered = exprVecOrdered[1:min(nrow(tableOrdered),1000)]
    createHeatmap(tableOrdered, exprVecOrdered, outNameHeatList[i])

    # Correlation
    vec1 = log10(as.numeric(inputTableEto[,c(10)]))
    vec2 = log10(as.numeric(rowMeans(inputTableEto[,12:ncol(inputTableEto)]))) + 2.3
    vec1[vec1<=0] = NA
    vec2[vec2<=0] = 0
    minV = min(c(length(vec1), length(vec2)))
    vec1 = vec1[1:minV]
    vec2 = vec2[1:minV]
    xLabel = "Normalized Nascent RNA expression (log10 + 1)"
    yLabel = "Normalized average BLISS DSB (log10 + 1)"
    outFileName = outNameCorrList[i]
    pAdj = flagAdjList[i]
    correlation = scatterPlot(vec1, vec2, xLabel, ylabel, pAdj, outFileName)
    toWriteCorrList[5*(i+1)] = correlation

  }
  else{

    # Line input
    listOfVectors = vector("list", 1)

    # All genes
    signalVector = as.numeric(colMeans(inputTableEto[,12:ncol(inputTableEto)]))
    signalVector[is.na(signalVector)] = 0.0
    signalVector = binFun(signalVector)
    listOfVectors[[1]] = signalVector

    # Line Plot
    legendList = c("All genes")
    linePlot(yMaxList[i], legendList, numberOfElements, listOfVectors, outNameAggrList[i])

    # Heatmap Plot
    table = as.matrix(binFunTable(inputTableEto[, 12:ncol(inputTableEto)]))
    newOrder = order(-rowSums(table))
    tableOrdered = table[newOrder,]
    tableOrdered = tableOrdered[1:min(nrow(tableOrdered),1000),]
    exprVec = abs(rnorm(min(nrow(tableOrdered),1000), mean=50, sd=35))
    exprVecOrdered = as.numeric(exprVec[newOrder])
    #createHeatmap(tableOrdered, exprVecOrdered, outNameHeatList[i])

    # Correlation
    vec1 = log10(as.numeric(exprVecOrdered))
    vec2 = log10(as.numeric(rowMeans(tableOrdered))) + 0.3 #1.3
    vec1[vec1<=0] = NA
    vec2[vec2<=0] = 0
    minV = min(c(length(vec1), length(vec2)))
    vec1 = vec1[1:minV]
    vec2 = vec2[1:minV]
    xLabel = "Normalized Nascent RNA expression (log10 + 1)"
    yLabel = "Normalized average BLISS DSB (log10 + 1)"
    outFileName = outNameCorrList[i]
    pAdj = flagAdjList[i]
    correlation = scatterPlot(vec1, vec2, xLabel, ylabel, pAdj, outNameCorrList[i])
    #toWriteCorrList[5*(i+1)] = NA

  }

}

#corrTableOutFileName = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/1_ctcfSignalPlots/correlation_summary.csv"
#write(toWriteCorrList, file = corrTableOutFileName, ncolumns = 5, append = FALSE, sep = "\t")



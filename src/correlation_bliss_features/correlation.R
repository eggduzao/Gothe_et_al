
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

# Input
args <- commandArgs(trailingOnly = TRUE)
graphWidth = as.numeric(args[1])
marginX = as.numeric(args[2])
inputTableFileName = args[3]
outputFileName = args[4]
outputLocation = args[5]
set.seed(13)

###################################################################################################
# Functions
###################################################################################################

# Scatter plot
scatterPlot <- function(vec1, vec2, xLabel, ylabel, outFileName){

  # Parameters
  dataFr = data.frame(X=vec1, Y=vec2)
  colVec = c("gray")
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

  # Calculating correlation
  corrTestSpearman = cor.test(vec1, vec2, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  pValueSpearman = corrTestSpearman$p.value
  corrSpearman = corrTestSpearman$estimate

  # Plotting graph
  pplot = ggplot(dataFr, aes(vec1,vec2)) 
  pplot = pplot + geom_point(alpha=0.1) 
  pplot = pplot + xlab(xLabel) 
  pplot = pplot + ylab(ylabel)
  pplot = pplot + ggtitle(paste('Spearman = ',round(corrSpearman, digits = 4)))
  pplot = pplot + theme_classic()
  pplot = pplot + stat_density_2d(aes(,fill=..level..), bins=11, geom = "polygon")
  pplot = pplot + scale_fill_gradientn(colours = myPalette(11))
  pplot = pplot + theme(plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 5, height = 5)

}

# Line plot
linePlot <- function(vec1, vec2, names, graphWidth, marginX, outFileName){

  # Initilize figure
  pdf(file = outFileName, width = graphWidth, height = 10)
  par(mar = c(marginX,5,1,1), cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, cex.sub = 2.0)

  # Parameters
  vec1 = vec1[!is.na(vec1)]
  vec2 = vec2[!is.na(vec2)]
  pchVec = c(19, 15)
  colVec = c("coral4", "dodgerblue")

  # Graph
  xRange = c(1,length(names))
  yRange = c(-0.21,1.01)
  plot(xRange, yRange, type="n", ylab="Spearman Correlation", xlab=" ", axes = FALSE, cex = 2.0)

  # Axis
  axis(side = 1, at=seq(from = 1, to = length(names), by = 1), labels = names, cex = 2.0, srt=45, las=2)
  axis(side = 2, at=seq(from = -0.2, to = 1.0, by = 0.1), cex = 2.0)

  # Points
  points(vec1, type = "p", pch = pchVec[1], col = colVec[1], cex = 3.0)
  points(vec2, type = "p", pch = pchVec[2], col = colVec[2], cex = 3.0)
  grid()

  # Legend
  legend("topright", inset = c(0.1,0.05), legend = c("BLISS", "Random"), pch = pchVec, col = colVec, cex = 2.0)

  # Control line
  abline(h = vec1[6], col = "lightgoldenrod4", lty = 2, cex = 0.5)

  # Termination
  dev.off()

}

# Calculate corrrelation
spearman <- function(vec1, vec2){

  # Remove NAs
  vec1 = vec1[!is.na(vec1)]
  vec2 = vec2[!is.na(vec2)]

  # Calculating correlation
  corrTestSpearman = cor.test(vec1, vec2, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  pValueSpearman = corrTestSpearman$p.value
  corrSpearman = corrTestSpearman$estimate
  return(c(corrSpearman,pValueSpearman))

}

###################################################################################################
# Execution
###################################################################################################

# Reading features table
inputTable = read.table(inputTableFileName, header = TRUE, sep = "\t")
blissVec = log10(inputTable[,1])

# Calculating correlations
featureVec = c()
randVec = c()
nameVec = c()
for(i in 2:ncol(inputTable)){
  col2 = log10(inputTable[,i])
  nameVec = c(nameVec, colnames(inputTable)[i])
  col3 = sample(blissVec)
  col4 = sample(col2)
  featureVec = c(featureVec, spearman(blissVec, col2)[1])
  randVec = c(randVec, spearman(col3, col4)[1])
  outFileNameCorrelationPlot = paste(outputLocation,colnames(inputTable)[i],".pdf",sep="")
  scatterPlot(blissVec, col2, "DSB Counts", paste(colnames(table)[i]," Counts",sep=""), outFileNameCorrelationPlot)
}

# Reordering table
featureOrder = order(featureVec, decreasing=TRUE)
featureVec = featureVec[featureOrder]
randVec = randVec[featureOrder]
nameVec = nameVec[featureOrder]

# Create plot
linePlot(featureVec, randVec, nameVec, graphWidth, marginX, outputFileName)

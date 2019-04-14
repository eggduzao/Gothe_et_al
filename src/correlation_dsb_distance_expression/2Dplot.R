
###################################################################################################
# Import
###################################################################################################

# Import
rm(list=ls())
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)
library(plotly)
library(scatterplot3d)

# Input
args <- commandArgs(trailingOnly = TRUE)
max_dist = as.numeric(args[1])
inputTableFileName = args[2]
outFileNameDD = args[3]
outFileNameDE = args[4]
outFileNameED = args[5]

###################################################################################################
# Functions
###################################################################################################

# Correlation 2D
corrPlot2D <- function(vecX, vecY, xLab, yLab, initialText, breakVec, y1, y2, outFileName){

  # Initialize plot
  dataFr = data.frame(X = vecX, Y = vecY)
  if(breakVec[1] == 0){labelVec = breakVec / 2}
  else{labelVec = breakVec}

  # Spearman Correlation
  corrTestSpearman = cor.test(vecX, vecY, alternative = "two.sided", method = "spearman", conf.level = 0.95)
  correlation = corrTestSpearman$estimate
  
  # Plotting graph
  pplot = ggplot()
  pplot = pplot + geom_point(data=dataFr, aes(x=X, y=Y), color="gray45", shape=21, size=1, alpha = 0.2)
  pplot = pplot + xlab(xLab) 
  pplot = pplot + ylab(yLab)
  pplot = pplot + ggtitle(paste(initialText,"\nCorrelation = ",round(correlation, digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + scale_x_continuous(breaks = breakVec, labels = labelVec, limits = c(breakVec[1], breakVec[length(breakVec)]))
  pplot = pplot + scale_y_continuous(limits = c(y1,y2))
  pplot = pplot + theme(legend.title=element_text(size=8), plot.title = element_text(hjust = 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 6, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Labels
distLabel = "Distance from the closest loop anchor (Kbp)"
expLabel = "Expression (log10 RPKM)"
dsbLabel = "Average DSBs per gene"

# Reading table
table = as.matrix(read.table(inputTableFileName, header = TRUE))

# Vector X
vectorX = as.numeric(table[,"DISTANCE"])

# Vector Y
vectorY = log10(as.numeric(table[,"EXPRESSION"]))

# Vector Z
vectorZ = as.numeric(table[,"DSB"]) * 100

# Distance vs DSB
initialText = "Distance vs DSBs"
corrPlot2D(vectorX, vectorZ, distLabel, dsbLabel, initialText, seq(0, max_dist, 5), 0, 35, outFileNameDD)

# Distance vs Expression
initialText = "Distance vs Expression"
corrPlot2D(vectorX, vectorY, distLabel, expLabel, initialText, seq(0, max_dist, 5), -2, 5, outFileNameDE)

# Expression vs DSB
initialText = "Expression vs DSBs"
corrPlot2D(vectorY, vectorZ, expLabel, dsbLabel, initialText, c(-2, -1, 0, 1, 2, 3, 4, 5), 0, 35, outFileNameED)



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
outFileName = args[3]

###################################################################################################
# Functions
###################################################################################################

# Correlation 3D
corrPlot3D <- function(vecX, vecY, vecZ, initialText, breakVec, outFileName){

  # Initialize plot
  dataFr = data.frame(X = vecX, Y = vecY, Z = vecZ)

  # Calculating multiple regression
  mreg = lm(Z~X+Y, dataFr)
  mregs = summary(mreg)
  correlation = sqrt(sqrt(sqrt(sqrt(mregs$adj.r.squared))))
  
  # Labels are half the breaks
  labelVec = breakVec / 2
  
  # Plotting graph
  pplot = ggplot()
  pplot = pplot + geom_point(data=dataFr, aes(x=X, y=Z, color=Y), shape=21, size=1, alpha = 0.25)
  pplot = pplot + scale_colour_gradient2(na.value = "darkgreen", low = "darkgreen", mid = "darkgreen", high = "red", midpoint = 0, limits = c(-1.5, 4.5), breaks = c(-2, 5), labels = c(-2, 5), name = "Expression\n(log10 RPKM)", guide = "colourbar")
  pplot = pplot + xlab("Distance from the closest loop anchor (Kbp)") 
  pplot = pplot + ylab("Average DSBs per gene")
  pplot = pplot + ggtitle(paste(initialText,"\nCorrelation = ",round(correlation, digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + scale_x_continuous(breaks = breakVec, labels = labelVec, limits = c(breakVec[1], breakVec[length(breakVec)]))
  pplot = pplot + ylim(0, 35)
  pplot = pplot + theme(legend.title=element_text(size=8), plot.title = element_text(hjust = 0.5))
  pplot = pplot + guides(colour = guide_colourbar(order = 1))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 6, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Reading tabl
table = as.matrix(read.table(inputTableFileName, header = TRUE))

# Vector X
vectorX = as.numeric(table[,"DISTANCE"])

# Vector Y
vectorY = log10(as.numeric(table[,"EXPRESSION"]))

# Vector Z
vectorZ = as.numeric(table[,"DSB"]) * 100

# 3D Correlation
initialText = "Distance vs DSBs vs Expression"
corrPlot3D(vectorX, vectorY, vectorZ, initialText, seq(0, max_dist, 20), outFileName)

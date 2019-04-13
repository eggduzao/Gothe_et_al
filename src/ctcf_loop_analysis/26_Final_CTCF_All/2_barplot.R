
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/2_hotspotBarplot/"
inputFileList = c()
outFileList = c()
yLimList = c()
foldList = c("all", "a", "b")
for(fold in foldList){
  subfoldList = c("active_genes", "all_genes", "inactive_genes", "intergenic", "active_promoters", "all_promoters", "inactive_promoters")
  for(subf in subfoldList){
    cellList = c("K562", "TK6")
    for(cell in cellList){
      inputFileList = c(inputFileList, paste(il,fold,"/",cell,"_",subf,".txt",sep=""))
      outFileList = c(outFileList, paste(il,fold,"/",cell,"_",subf,".pdf",sep=""))
      if(fold == "all"){yLimList = c(yLimList, 37)}
      else if(fold == "a"){yLimList = c(yLimList, 37)}
      else if(fold == "b"){yLimList = c(yLimList, 18)}
    }
  }
}

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecZ, vecC, yLim, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ, vectorc = vectorC)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, fill=vectorz))
  pplot = pplot + geom_bar(stat="identity", position=position_dodge())
  pplot = pplot + theme_classic()
  pplot = pplot + ylim(0, yLim)
  pplot = pplot + xlab(" ")
  pplot = pplot + ylab("Frequency of Regions (%)")
  pplot = pplot + geom_text(aes(label=as.character(paste(vectory,"%",sep=""))), vjust = -0.5, color="black", position = position_dodge(0.9), size=3.5)
  pplot = pplot + geom_text(aes(label=as.character(paste("(",vectorc,")",sep=""))), vjust = -2.0, color="black", position = position_dodge(0.9), size=3.5)
  pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10), 
                        axis.text.y = element_text(size=14), axis.title=element_text(size=16))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 8, height = 6)

}

###################################################################################################
# Execution
###################################################################################################

# Iterating on columns
for(i in 1:length(inputFileList)){

  # Reading table
  inputFileName = inputFileList[i]
  table = read.table(inputFileName, header = TRUE)

  # Vector X
  vectorX = as.character(table[,"X"])

  # Vector Y
  vectorY = as.numeric(table[,"Y"])

  # Vector Z
  vectorZ = as.character(table[,"Z"])

  # Vector C
  vectorC = as.numeric(table[,"C"])

  # Barplot
  yLim = yLimList[i]
  outFileName = outFileList[i]
  barPlot(vectorX, vectorY, vectorZ, vectorC, yLim, outFileName)
}



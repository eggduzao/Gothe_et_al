
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
il = "/media/egg/cheops_agpapan/eduardo/Roukos_Bliss/Results/26_Final_CTCF_All/3_hotspotViolinplot/"
inputFileList = c()
outFileList = c()
foldList = c("all", "a", "b")
for(fold in foldList){
  subfoldList = c("active_genes", "all_genes", "inactive_genes", "intergenic", "active_promoters", "all_promoters", "inactive_promoters")
  for(subf in subfoldList){
    cellList = c("K562", "TK6")
    for(cell in cellList){
      inputFileList = c(inputFileList, paste(il,fold,"/",cell,"_",subf,".txt",sep=""))
      outFileList = c(outFileList, paste(il,fold,"/",cell,"_",subf,".pdf",sep=""))
    }
  }
}

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecZ, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, color=vectorz))
  #pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72))
  pplot = pplot + geom_boxplot(width=0.1, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  pplot = pplot + ylim(-0.5, 5)
  pplot = pplot + xlab(" ")
  pplot = pplot + ylab("Total Signal Intensity (log10)")
  pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10), 
                        axis.text.y = element_text(size=14), axis.title=element_text(size=16), 
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  #pplot = pplot + scale_y_continuous(minor_breaks = seq(-2 , 10, 0.1), breaks = seq(-2, 10, 0.5))
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
  vectorY = log(as.numeric(table[,"Y"]))

  # Vector Z
  vectorZ = as.character(table[,"Z"])

  # Barplot
  outFileName = outFileList[i]
  barPlot(vectorX, vectorY, vectorZ, outFileName)

}



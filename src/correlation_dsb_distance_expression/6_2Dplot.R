
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

# Input Lists
il = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/5_Tables/"
ol = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/6_2D_Plots_DiffCells/"
parameterCellList = c("K562", "TK6", "CD34")
parameterRegionList = c("promoter_gene", "gene", "promoter")
parameterLoopTypeList = c("anchors_with_ctcf", "anchors_with_and_wo_ctcf")
parameterMaxDistList = c("201", "501", "1001")
inputDsbCellList = c()
inputExpCellList = c()
inputDistCellList = c()
outLocList = c()
outputRegion = c()
outputMaxDist = c()
outputCell1List = c()
outputCell2List = c()
outputCell3List = c()

# Input Loop
for(cell1 in parameterCellList){
  for(cell2 in parameterCellList){
    for(cell3 in parameterCellList){
      for(region in parameterRegionList){
        for(loop in parameterLoopTypeList){
          for(maxdist in parameterMaxDistList){
            inputDsbCellList = c(inputDsbCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell1,"_table_filt.txt",sep=""))
            inputExpCellList = c(inputExpCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell2,"_table_filt.txt",sep=""))
            inputDistCellList = c(inputDistCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell3,"_table_filt.txt",sep=""))
            outLocList = c(outLocList, paste(ol,region,"/",loop,"/",maxdist,"/",sep=""))
            outputRegion = c(outputRegion, region)
            outputMaxDist = c(outputMaxDist, maxdist)
            outputCell1List = c(outputCell1List, cell1)
            outputCell2List = c(outputCell2List, cell2)
            outputCell3List = c(outputCell3List, cell3)
          }
        }
      }
    }
  }
}

###################################################################################################
# Functions
###################################################################################################

# Correlation 2D
corrPlot2D <- function(vecX, vecY, vecM, vecC, vecL, xLab, yLab, hjustVec, vjustVec, initialText, mllPointSize, breakVec, y1, y2, outFileName){

  # Initialize plot
  dataFr = data.frame(X = vecX, Y = vecY, M = vecM, C = vecC, L = vecL)
  dataFr1 = subset(dataFr, M == "all")
  dataFr2 = subset(dataFr, M == "mll")
 
  # Labels are half the breaks
  if(xLab == distLabel){
    labelVec = breakVec / 2
  } else {
    labelVec = breakVec
  }

  # Spearman Correlation
  corrTestSpearman = cor.test(vecX, vecY, alternative = "two.sided", method = "spearman", conf.level = 0.95)
  correlation = corrTestSpearman$estimate
  
  # Plotting graph
  pplot = ggplot()
  pplot = pplot + geom_point(data=dataFr1, aes(x=X, y=Y), color="gray45", shape=21, size=1, alpha = 0.2)
  if(mllPointSize == TRUE){
    pplot = pplot + geom_point(data=dataFr2, aes(x=X, y=Y, shape=shp, size=C), shape=21, alpha = 1.0, fill = "blue4", show.legend = F)
    pplot = pplot + geom_text(data=dataFr2, aes(x = X, y = Y, label = L), hjust=hjustVec, vjust=vjustVec, size = 3, color = "blue4", fontface = "bold")
  } else{
    pplot = pplot + geom_point(data=dataFr2, aes(x=X, y=Y, shape=shp), shape=21, alpha = 1.0, fill = "blue4", show.legend = F)
  }
  pplot = pplot + xlab(xLab) 
  pplot = pplot + ylab(yLab)
  pplot = pplot + ggtitle(paste(initialText,"\nCorrelation = ",round(correlation, digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + scale_x_continuous(breaks = breakVec, labels = labelVec)
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

# MLL Point size
mllPointSizeDD = TRUE
mllPointSizeDE = TRUE
mllPointSizeED = FALSE

# Y ranges
y1DD = 0; y2DD = 35
y1DE = -2; y2DE = 5
y1ED = 0; y2ED = 35

# Iterating on all possibilities
for(i in 1:length(outputRegion)){

  # Reading input
  inputDsbFileName = inputDsbCellList[i]
  inputExpFileName = inputExpCellList[i]
  inputDistFileName = inputDistCellList[i]
  outLoc = outLocList[i]
  outputRegion = outputRegion[i]
  maxD = outputMaxDist[i]
  cell1 = outputCell1List[i]
  cell2 = outputCell2List[i]
  cell3 = outputCell3List[i]

  if (cell1 == cell2 & cell1 == cell3 & cell2 == cell3) {
    next
  }

  # Creating output location
  system(paste("mkdir -p ",outLoc,sep=""))

  # Reading table
  tableDsb = as.matrix(read.table(inputDsbFileName, header = TRUE))
  tableExp = as.matrix(read.table(inputExpFileName, header = TRUE))
  tableDist = as.matrix(read.table(inputDistFileName, header = TRUE))

  # Vector X
  vectorXDD = as.numeric(tableDist[,"X"])
  vectorXDE = as.numeric(tableDist[,"X"])
  vectorXED = log10(as.numeric(tableExp[,"Y"]))

  # Vector Y
  vectorYDD = as.numeric(tableDsb[,"Z"]) * 100
  vectorYDE = log10(as.numeric(tableExp[,"Y"]))
  vectorYED = as.numeric(tableDsb[,"Z"]) * 100

  # Vector M
  vectorMDD = as.character(tableDsb[,"M"])
  vectorMDE = as.character(tableDsb[,"M"])
  vectorMED = as.character(tableDsb[,"M"])

  # Vector C
  vectorCDD = as.numeric(tableDsb[,"C"])
  vectorCDE = as.numeric(tableDsb[,"C"])
  vectorCED = as.numeric(tableDsb[,"C"])

  # Vector L
  vectorLDD = c()
  vectorLDE = c()
  vectorLED = c()
  for(k in 1:length(vectorCDD)){
    if(vectorCDD[k] > 1){
      vectorLDD = c(vectorLDD, as.character(vectorCDD[k]))
      vectorLDE = c(vectorLDE, as.character(vectorCDD[k]))
      vectorLED = c(vectorLED, " ")
    } else{
      vectorLDD = c(vectorLDD, " ")
      vectorLDE = c(vectorLDE, " ")
      vectorLED = c(vectorLED, " ")
    }
  }

  # Horizontal and vertical positioning
  mllCounter = length(vectorMDD[vectorMDD == "mll"])
  hjustVec = rep(-0.5, mllCounter)
  vjustVec = rep(-1.0, mllCounter)

  # Initial tests
  initialTextDD = paste("Distance = ",cell3," / DSBs = ",cell1)
  initialTextDE = paste("Distance = ",cell3," / Expression = ",cell2)
  initialTextED = paste("Expression = ",cell2," / DSBs = ",cell1)

  # Break vector
  if(maxD == "201"){
    breakVecDD = c(0, 50, 100, 150, 200)
    breakVecDE = c(0, 50, 100, 150, 200)
    breakVecED = c(-2, -1, 0, 1, 2, 3, 4, 5)
  } else if(maxD == "501"){
    breakVecDD = c(0, 125, 250, 375, 500)
    breakVecDE = c(0, 125, 250, 375, 500)
    breakVecED = c(-2, -1, 0, 1, 2, 3, 4, 5)
  } else{
    breakVecDD = c(0, 250, 500, 750, 1000)
    breakVecDE = c(0, 250, 500, 750, 1000)
    breakVecED = c(-2, -1, 0, 1, 2, 3, 4, 5)
  }

  # Output file name
  outputFileNameDD = paste(outLoc,"DISTvsDSB - dist ",cell3," - dsb ",cell1,".pdf",sep="")
  outputFileNameDE = paste(outLoc,"DISTvsEXPR - dist ",cell3," - expr ",cell2,".pdf",sep="")
  outputFileNameED = paste(outLoc,"EXPRvsDSB - expr ",cell2," - dsb ",cell1,".pdf",sep="")

  # Crop to min value
  minDD = min(length(vectorXDD),length(vectorYDD),length(vectorMDD),length(vectorCDD),length(vectorLDD))
  vectorXDD = vectorXDD[(length(vectorXDD) - minDD + 1):length(vectorXDD)]
  vectorYDD = vectorYDD[(length(vectorYDD) - minDD + 1):length(vectorYDD)]
  vectorMDD = vectorMDD[(length(vectorMDD) - minDD + 1):length(vectorMDD)]
  vectorCDD = vectorCDD[(length(vectorCDD) - minDD + 1):length(vectorCDD)]
  vectorLDD = vectorLDD[(length(vectorLDD) - minDD + 1):length(vectorLDD)]
  minDE = min(length(vectorXDE),length(vectorYDE),length(vectorMDE),length(vectorCDE),length(vectorLDE))
  vectorXDE = vectorXDE[(length(vectorXDE) - minDE + 1):length(vectorXDE)]
  vectorYDE = vectorYDE[(length(vectorYDE) - minDE + 1):length(vectorYDE)]
  vectorMDE = vectorMDE[(length(vectorMDE) - minDE + 1):length(vectorMDE)]
  vectorCDE = vectorCDE[(length(vectorCDE) - minDE + 1):length(vectorCDE)]
  vectorLDE = vectorLDE[(length(vectorLDE) - minDE + 1):length(vectorLDE)]
  minED = min(length(vectorXED),length(vectorYED),length(vectorMED),length(vectorCED),length(vectorLED))
  vectorXED = vectorXED[(length(vectorXED) - minED + 1):length(vectorXED)]
  vectorYED = vectorYED[(length(vectorYED) - minED + 1):length(vectorYED)]
  vectorMED = vectorMED[(length(vectorMED) - minED + 1):length(vectorMED)]
  vectorCED = vectorCED[(length(vectorCED) - minED + 1):length(vectorCED)]
  vectorLED = vectorLED[(length(vectorLED) - minED + 1):length(vectorLED)]

  # 2D Correlation
  corrPlot2D(vectorXDD, vectorYDD, vectorMDD, vectorCDD, vectorLDD, distLabel, dsbLabel, hjustVec, vjustVec, initialTextDD,
             mllPointSizeDD, breakVecDD, y1DD, y2DD, outputFileNameDD)
  corrPlot2D(vectorXDE, vectorYDE, vectorMDE, vectorCDE, vectorLDE, distLabel, expLabel, hjustVec, vjustVec, initialTextDE,
             mllPointSizeDE, breakVecDE, y1DE, y2DE, outputFileNameDE)
  corrPlot2D(vectorXED, vectorYED, vectorMED, vectorCED, vectorLED, expLabel, dsbLabel, hjustVec, vjustVec, initialTextED,
             mllPointSizeED, breakVecED, y1ED, y2ED, outputFileNameED)

}



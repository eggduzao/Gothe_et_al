
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
parameterCellList = c("K562", "TK6", "CD34")
parameterRegionList = c("promoter_gene", "gene", "promoter")
parameterLoopTypeList = c("anchors_with_ctcf", "anchors_with_and_wo_ctcf")
parameterMaxDistList = c("201", "501", "1001")
inputDsbCellList = c()
inputExpCellList = c()
inputDistCellList = c()
outputFileList = c()
outputMaxDist = c()
outputCell1List = c()
outputCell2List = c()
outputCell3List = c()
outSquare = c()

# Input Loop
for(cell1 in parameterCellList){
  for(cell2 in parameterCellList){
    for(cell3 in parameterCellList){
      for(region in parameterRegionList){
        for(loop in parameterLoopTypeList){
          for(maxdist in parameterMaxDistList){

            if(cell1 == cell2 & cell1 == cell3 & cell2 == cell3){
              ol = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/7_3D_Plots_SameCells/"
              outLoc = paste(ol,region,"/",loop,"/",maxdist,"/",sep="")
              system(paste("mkdir -p ",outLoc,sep=""))
              outSquare = c(outSquare, TRUE)
            } else {
              ol = "/Users/egg/Projects/Roukos_Bliss/Results/49_Final_3DPlot/7_3D_Plots_DiffCells/"
              outLoc = paste(ol,region,"/",loop,"/",maxdist,"/",sep="")
              system(paste("mkdir -p ",outLoc,sep=""))
              outSquare = c(outSquare, FALSE)
            }

            inputDistCellList = c(inputDistCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell1,"_table_filt.txt",sep="")) 
            inputExpCellList = c(inputExpCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell2,"_table_filt.txt",sep=""))
            inputDsbCellList = c(inputDsbCellList, paste(il,region,"/",loop,"/",maxdist,"/",cell3,"_table_filt.txt",sep=""))
            name = paste("dist ",cell1," - exp ",cell2," - dsb ",cell3,".pdf",sep="")
            outputFileList = c(outputFileList, paste(outLoc,name,".pdf",sep=""))
            outputCell1List = c(outputCell1List, cell1)
            outputCell2List = c(outputCell2List, cell2)
            outputCell3List = c(outputCell3List, cell3)
            outputMaxDist = c(outputMaxDist, maxdist)

          }
        }
      }
    }
  }
}

###################################################################################################
# Functions
###################################################################################################

# Correlation 3D
corrPlot3D <- function(vecX, vecY, vecZ, vecM, vecC, vecL, hjustVec, vjustVec, initialText, breakVec, osq, outFileName){

  # Initialize plot
  dataFr = data.frame(X = vecX, Y = vecY, Z = vecZ, M = vecM, C = vecC, L = vecL)
  dataFr1 = subset(dataFr, M == "all")
  dataFr2 = subset(dataFr, M == "mll")

  # Calculating multiple regression
  mreg = lm(Z~X+Y, dataFr)
  mregs = summary(mreg)
  if(osq){
    correlation = sqrt(mregs$adj.r.squared)
  } else {
    correlation = mregs$adj.r.squared**2
    if(is.nan(correlation) | is.infinite(correlation)){
      correlation = 0
    }
  }
  
  # Labels are half the breaks
  labelVec = breakVec / 2
  
  # Plotting graph
  pplot = ggplot()
  pplot = pplot + geom_point(data=dataFr1, aes(x=X, y=Z, color=Y), shape=21, size=1, alpha = 0.25)
  pplot = pplot + scale_colour_gradient2(na.value = "darkgreen", low = "darkgreen", mid = "darkgreen", high = "red", midpoint = 0, limits = c(-1.5, 4.5), breaks = c(-2, 5), labels = c(-2, 5), name = "Expression\n(log10 RPKM)", guide = "colourbar") # 3
  pplot = pplot + geom_point(data=dataFr2, aes(x=X, y=Z, shape=shp, size=C), shape=21, alpha = 1.0, fill = "blue4", show.legend = F)
  pplot = pplot + geom_text(data=dataFr2, aes(x = X, y = Z, label = L), hjust=hjustVec, vjust=vjustVec, size = 3, color = "blue4", fontface = "bold")
  pplot = pplot + xlab("Distance from the closest loop anchor (Kbp)") 
  pplot = pplot + ylab("Average DSBs per gene")
  pplot = pplot + ggtitle(paste(initialText,"\nCorrelation = ",round(correlation, digits = 4),sep=''))
  pplot = pplot + theme_classic()
  pplot = pplot + scale_x_continuous(breaks = breakVec, labels = labelVec)
  pplot = pplot + ylim(0, 35)
  pplot = pplot + theme(legend.title=element_text(size=8), plot.title = element_text(hjust = 0.5))
  pplot = pplot + guides(colour = guide_colourbar(order = 1))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 6, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Iterating on all possibilities
for(i in 1:length(outputMaxDist)){

  # Reading input
  inputDistFileName = inputDistCellList[i]
  inputExpFileName = inputExpCellList[i]
  inputDsbFileName = inputDsbCellList[i]
  outputFileName = outputFileList[i]
  cell1 = outputCell1List[i]
  cell2 = outputCell2List[i]
  cell3 = outputCell3List[i]
  maxD = outputMaxDist[i]
  osq = outSquare[i]

  # Reading tabl
  tableDist = as.matrix(read.table(inputDistFileName, header = TRUE))
  tableExp = as.matrix(read.table(inputExpFileName, header = TRUE))
  tableDsb = as.matrix(read.table(inputDsbFileName, header = TRUE))

  # Vector X
  vectorX = as.numeric(tableDist[,"X"])

  # Vector Y
  vectorY = log10(as.numeric(tableExp[,"Y"]))

  # Vector Z
  vectorZ = as.numeric(tableDsb[,"Z"]) * 100

  # Vector M
  vectorM = as.character(tableDsb[,"M"])

  # Vector C
  vectorC = as.numeric(tableDsb[,"C"])

  # Vector L
  vectorL = c()
  for(k in 1:length(vectorC)){
    if(vectorC[k] > 1){
      vectorL = c(vectorL, as.character(vectorC[k]))
    } else{
      vectorL = c(vectorL, " ")
    }
  }

  # Horizontal and vertical positioning
  mllCounter = length(vectorM[vectorM == "mll"])
  hjustVec = rep(-0.5, mllCounter)
  vjustVec = rep(-1.0, mllCounter)

  # Initial tests
  initialText = paste("Dist = ",cell1," / Exp = ",cell2," / DSB = ",cell3)

  # Break vector
  if(maxD == "201"){
    breakVec = c(0, 50, 100, 150, 200)
  } else if(maxD == "501"){
    breakVec = c(0, 125, 250, 375, 500)
  } else{
    breakVec = c(0, 250, 500, 750, 1000)
  }

  # Crop to min value
  minV = min(length(vectorX),length(vectorY),length(vectorZ),length(vectorM),length(vectorC),length(vectorL))
  vectorX = vectorX[(length(vectorX) - minV + 1):length(vectorX)]
  vectorY = vectorY[(length(vectorY) - minV + 1):length(vectorY)]
  vectorZ = vectorZ[(length(vectorZ) - minV + 1):length(vectorZ)]
  vectorM = vectorM[(length(vectorM) - minV + 1):length(vectorM)]
  vectorC = vectorC[(length(vectorC) - minV + 1):length(vectorC)]
  vectorL = vectorL[(length(vectorL) - minV + 1):length(vectorL)]

  # 3D Correlation
  corrPlot3D(vectorX, vectorY, vectorZ, vectorM, vectorC, vectorL, hjustVec, vjustVec, initialText, breakVec, osq, outputFileName)

}


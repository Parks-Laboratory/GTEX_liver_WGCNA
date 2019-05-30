

### 1 loading library


library(WGCNA)
library(dplyr)

options(stringsAsFactors = FALSE);

load ("network_expr_data.Rdata")

femData = datExpr;
femData <- femData * 1.0
femData = data.frame(t(femData))


maleData = datExpr;
maleData <- maleData * 1.0
maleData = data.frame(t(maleData))

## name data

dim(femData)
names(femData)
dim(maleData)
names(maleData)

### 2  work with two sets
nSets = 2

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")


# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(femData)))
multiExpr [[1]]
names(multiExpr[[1]]$data) = femData$gene
rownames(multiExpr[[1]]$data) = names(femData)

multiExpr[[2]] = list(data = as.data.frame(t(maleData)))
multiExpr [[2]]

names(multiExpr[[2]]$data) = maleData$gene
rownames(multiExpr[[2]]$data) = names(maleData)


exprSize = checkSets(multiExpr)
exprSize


####   sample plot tree to find outliers 
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "0-gtex-SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);

dev.off()



# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples


save.image (file = "0-gtex-Consensus-dataInput.RData")


library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load (file = "0-gtex-Consensus-dataInput.RData");


# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(4,6,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();

####  works fine

# Plot the results:
colors = c("black", "red")


# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}


# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "0-scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()


#########  one step network construction and module detection

# original
# net = blockwiseConsensusModules(
#  multiExpr, power = 6, TOMType = "unsigned", minModuleSize = 30, reassignedThreshold = 0,
#  pamRespectsDendro = FALSE, 
#  mergeCutHeight = 0.25, numericLabels = TRUE,
#  minKMEtoStay = 0,
#  saveTOMs = TRUE, verbose = 5)



net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)



consMEs = net$multiMEs
consMEs

moduleLabels = net$colors;
moduleLabels


#  Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
listofgene = list ("color"=moduleColors,"MES"=consMEs)

consTree = net$dendrograms[[1]]; 

# sizeGrWindow(8,6);
# pdf(file = "0-ConsensusDendrogram-auto.pdf", width = 8, height = 6)
# plotDendroAndColors(consTree, moduleColors,
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Consensus gene dendrogram and module colors")
# dev.off()


save.image (file="1-gtex-Consensus-NetworkConstruction-auto.RData")



library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "0-gtex-Consensus-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "1-gtex-Consensus-NetworkConstruction-auto.RData");
lnames

exprSize = checkSets(multiExpr)
nSets = exprSize$nSets

info = data.frame(GeneName = GeneName,
                  femData = femData,
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels));


info <- info [order(info$ModuleColor),]

write.csv(info, file = "1-gtex-consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE)

save.image ("2-gtex-InfoWithModule.RData")



getwd();

setwd("c:/Users/Zongliang/Desktop/PD project/matrix/")

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", 
                   "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
##############################
setwd("C:\\Users\\Zongliang\\Desktop\\PD project\\0.matrixCOPD")
#1.a Loading expression data
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
  femData = read.table("matrix.txt",header = T);
# Take a quick look at what is in the data set:
dim(femData);
names(femData);
matrix=femData[,2:ncol(femData)]
names(matrix);
rownames(matrix)=femData[1:nrow(femData),1];
rownames(matrix)
rowName=rownames(matrix)


##### ECLIPSE and COPD#####
RowName=femData[,1]
DataEclipse=femData[,2:229]
DataCOPD=femData[,230:ncol(femData)]
rownames(DataEclipse)=RowName
rownames(DataCOPD)=RowName

##### 1.check samples ####
gsgEclipse = goodSamplesGenes(DataEclipse, verbose = 3);
gsgCOPD = goodSamplesGenes(DataCOPD, verbose = 3);
gsgEclipse$allOK
gsgCOPD$allOK
DataEclipse=t(data.matrix(DataEclipse))

sampleTreeEclipse = hclust(dist(DataEclipse), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeEclipse, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 100, col = "red");
# Determine cluster under the line
clustEclipse = cutreeStatic(sampleTreeEclipse, cutHeight = 100, minSize = 1)
table(clustEclipse)
# clust 1 contains the samples we want to keep.
keepSamples = (clustEclipse==1)
datExprEclipse = DataEclipse[keepSamples, ]
### soft ####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprEclipse, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(datExprEclipse, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "parkinson",
                       verbose = 3)

##########################33
#1.b Checking data for excessive missing values and identification of outlier microarray samples
#We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(matrix, verbose = 3);
gsg$allOK
matrix=t(data.matrix(matrix))
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTree = hclust(dist(matrix), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf("cluster.pdf")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 100, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 1)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = matrix[keepSamples, ]



nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
###########################
#1.c Loading clinical trait data


##############################
#2.a Automatic network construction and module detection
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
#lnames = load(file = ".RData");
#The variable lnames contains the names of loaded variables.
#lnames
######################3
#2.a.1 Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#2.a.2 One-step network construction and module detection
net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "parkinson",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.3,
                    addGuide = TRUE, guideHang = 0.5)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
table(net$colors)
for (i in 0:max(net$colors)){
  a=which(net$colors==i)
  
  write.table(RowName[a],paste(i,".txt"),row.names=F,col.names=F,quote=F)
}
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Parkinson.RData")

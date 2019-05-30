
library(optparse)
library (WGCNA)
options (stringAsFactors = FALSE)


# load data

datMicroarrays = read.csv ("gtex.txt",header = TRUE, sep = "\t")
datMicroarrays [1:5,1:5]

# rename sample name  and gene name
ArrayName = names (data.frame(datMicroarrays[,-1]))
GeneName = datMicroarrays$GeneName


# exclude the first column
datExpr = data.frame (t(datMicroarrays[,-1]))


names (datExpr) = datMicroarrays [,1]
dimnames (datExpr)[[1]] = names (data.frame (datMicroarrays[,-1]))

datExpr [1:5,1:5]


#  get mean and missing value
meanExpressionByArray = apply (datExpr,1,mean,na.rm=T)

NumberMissingByArray = apply (is.na(data.frame(datExpr)),1,sum)

NumberMissingByArray

# show mean expression per array

sizeGrWindow (10,5)
barplot (meanExpressionByArray,xlab = "sample",ylab = "Mean expression", main =" gtex human liver microarray", names.arg = c(1:117),cex.names = 0.7)


# Keep arrays with less missing values

KeepArray = NumberMissingByArray < 4500
table (KeepArray)
datExpr = datExpr [KeepArray,]

# Y is the clinical trait
# y = y [KeepArray]
ArrayName [KeepArray]


# handling missing data and zero variance in probe profiles

NumberMissingByGene = apply (is.na(data.frame(datExpr)),2,sum)
summary (NumberMissingByGene)
no.presentdatExpr =as.vector (apply(!is.na(as.matrix(datExpr)),2,sum))
table (no.presentdatExpr)
KeepGenes = no.presentdatExpr >= 88
table (KeepGenes)

datExpr = datExpr [,KeepGenes]
GeneName = GeneName [KeepGenes] 

variancedatExpr = as.vector (apply(as.matrix(datExpr),2,var,na.rm=T))
KeepGenes = variancedatExpr > 0
table (KeepGenes)
datExpr = datExpr [,KeepGenes]
GeneName = GeneName [KeepGenes] 


row.names.remove <- c("a6674", "a6567","a6572","a6566")
datExpr <- datExpr[!(row.names(datExpr) %in% row.names.remove), ]


save.image (file="network_expr_data.Rdata")




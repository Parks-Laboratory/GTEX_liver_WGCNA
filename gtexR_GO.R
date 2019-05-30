
BiocManager::install('org.Mm.eg.db')
BiocManager::install('org.Hs.eg.db')
install.packages ("xml2")
BiocManager::install("biomaRt")
install.packages("sqldf")

library (org.Mm.eg.db)
library (org.Hs.eg.db)
library(biomaRt)
library (WGCNA)
library(sqldf)

load (file = "2-gtex-InfoWithModule.RData")

annot = read.csv ("GeneAnotation.csv",sep=",")
annot

probes = names (datExpr)
probes

probes2annot = match (probes,annot$GeneName)

probes2annot

ensembl.genes = annot$id [probes2annot]


# ensembl.genes <- readLines ('ensemble_id.csv')
ensembl.genes

ensembl.genes <- gsub('.{0,2}$', '', ensembl.genes)

ensembl.genes  <- gsub("[^0-9A-Za-z///' ]","#" , ensembl.genes ,ignore.case = TRUE)
ensembl.genes <- gsub("#","", ensembl.genes ,ignore.case = TRUE)
values <- as.vector (ensembl.genes)
values
# values [is.na(values)] <- "ENSG00000270040"
# values


# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="uswest.ensembl.org")
entrezgene = getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
      filters = 'ensembl_gene_id', 
      values = values, 
      mart = mart)

names(entrezgene)[2]<-"genenum"
entrezgene
 
entrezgene<- entrezgene[!duplicated(entrezgene[1]), ]
entrezgene

value_dataframe = data.frame (values)
value_dataframe
names(value_dataframe)[1]<-"ensembl_gene_id"
value_dataframe


# write.csv (entrezgene,"entrezgene.xls")
# write.csv (value_dataframe,"value_dataframe.xls")


combined_df <- sqldf("Select f.*, most.genenum
                      from value_dataframe f 
                      left JOIN (select distinct ensembl_gene_id,genenum from entrezgene) as most 
                      on f.ensembl_gene_id = most.ensembl_gene_id
                      ")

combined_df



# allLLIDs = as.data.frame.matrix(entrezgene) 
allLLIDs <- as.vector (combined_df[[2]])
allLLIDs


GOenr = GOenrichmentAnalysis (moduleColors,allLLIDs,organism="human",nBestP = 10)

tab = GOenr$bestPTerms [[4]]$enrichment

write.table (tab,"0-GOenrichementtable.csv",sep=",",quote= TRUE,row.names=FALSE)






















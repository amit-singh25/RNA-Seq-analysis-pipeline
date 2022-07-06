
```
##############################################################################
Load require Library,or installed.packages('name of the pkg')
#############################################################################

library(DESeq2)
library(biomaRt)
library(gage)
library(gageData)
library(org.Mm.eg.db)
library(reshape2)
library(dplyr)

###############################################################################
Load the meta file, meta file contain tab separated sample infomration 
################################################################################

setwd("~/Desktop/Data")
meta<-read.delim("meta_data.txt",header=T,sep="\t")
meta<-read.delim("meta.txt",header = T,sep="\t")
sampleTable <- data.frame(sampleName = meta$sampleNames, fileName = meta$sampleFiles, condition = meta$sampleCondition)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = ".",design = ~ condition)


####################################################################################
Data visulation by PCA plot
####################################################################################

DESeq2::plotPCA(vsd, intgroup=c("condition"),ntop=1000)
dev.copy2pdf(file='PCA_plot.pdf')
ddsHTSeq <- DESeq(DESeq2Table)
##Comparison conditions between two sample group  
res <- results(ddsHTSeq, contrast=c("condition", "Sample1", "Sample2"))
##Select significant Diffrential regulated gene
resSig <- subset(res, padj < 0.01)

########################################################################################################################
Annotate Ensembl id to gene symbol,gene biotype and entrezgene by biomart here mouse name taken as an example"
########################################################################################################################

resSig$ensembl <- sapply(strsplit(rownames(resSig), split="\\+" ), "[", 1 )
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id","mgi_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = resSig$ensembl,mart = ensembl)
idx <- match(resSig$ensembl, genemap$ensembl_gene_id)
resSig$name <- genemap$mgi_symbol[idx]
resSig$gene <- genemap$entrezgene[idx]
resSig$gene_biotype <- genemap$gene_biotype[idx]
resSig$symbol <- genemap$name_1006[idx]

############################################################################################
Remove noncoding genes from the data
###############################################################################################

idx<-grep("protein_coding", resSig$gene_biotype, fixed = TRUE)
resSig<-resSig[idx,]
resOrdered <- resSig[order(-resSig$log2FoldChange ),]
write.xlsx(data.frame(resOrdered), "DEG_list.xlsx",asTable = TRUE,row.names=T)

#########################################################################################################
GO pathways analysis using GAGE Here all gene expression used unlike only DE gene, here data normalized 
########################################################################################################

GeneCounts <- counts(DESeq2Table)
cnts<-GeneCounts
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))#cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+1)
range(cnts.norm)

#####################################################################################################################################
Make a GO term list from "org.Mm.eg.db" a mouse annotation package from bioconductor but you can also download "Misgdatabase"
or any other sourcse of the data base can be useed such as KEGG, reactome, consesuspathDB etc.
#######################################################################################################################################

mygo <- as.list(org.Mm.egGO2EG)
mygo <- (mygo[!is.na(mygo)])
t <- mget(names(mygo),GOTERM)
names(mygo) <- as.character(lapply(t,Term))

########################################################################################################################
Obtain Significance pathways using gage, gage uses t-test between two condition and listed up and down Go/pathway
#########################################################################################################################
gos <- gage(mymat,gsets = mygo,ref=ref,samp=samp,compare ="paired")
gene_set<-sigGeneSet(gos,cutoff=0.001)
up<-(gene_set$greater)
write.csv(up,file="go_up_regulated.csv")
down<-(gene_set$less)
write.csv(down,file="go_down_regulated.csv")

########################################################################################################################
only with differentially regulated gene based on fold change 
########################################################################################################################
foldchanges = resSig$log2FoldChange
names(foldchanges) = resSig$ensembl
keggres = gage(foldchanges, gsets=mygo, same.dir=TRUE)
lapply(keggres, head)

```

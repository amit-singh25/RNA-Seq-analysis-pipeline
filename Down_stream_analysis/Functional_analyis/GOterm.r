library(gage)
library(org.Mm.eg.db)
library(genefilter)
library(annotate)
library("GO.db")
library("GOstats")
library(KEGG.db)

#####load Deseq count table 
GeneCounts <- counts(DESeq2Table)
cnts<-GeneCounts
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+1)
range(cnts.norm)

#### Name ENSEMBL ID TO GENEID

##### Most of the database identifies gene id, so I converted gene id from ensembl-idusing (org.Mm.eg.db) mouse annotataion pacakge

rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm),org.Mm.egENSEMBL2EG,ifnotfound=NA))
###Replace NA from rowname of the data, basically there some ensembl id which has no enterez id so i removed them in below command

mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])

###### Make a GO term list from "org.Mm.eg.db" a mouse annotation package from bioconductor but you can also download from my fav "Misgdatabase"
mygo <- as.list(org.Mm.egGO2EG)
mygo <- (mygo[!is.na(mygo)])
t <- mget(names(mygo),GOTERM)
names(mygo) <- as.character(lapply(t,Term))



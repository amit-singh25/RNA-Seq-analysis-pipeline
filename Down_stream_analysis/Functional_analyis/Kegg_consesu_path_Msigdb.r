library(gage)
library(gageData)

###load the cound data 

cnts <- counts(DESeq2Table)
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
dim(cnts)
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm=log2(cnts.norm+1)
range(cnts.norm)

rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm),org.Mm.egENSEMBL2EG,ifnotfound=NA))
###Replace NA from rowname of the data, basically there some ensembl id which has no enterez id so i removed them in below command
mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])

####kegg database with entrez gene id
kg.hsa=kegg.gsets("mmu")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]

###significant test

gos <- gage(mymat,gsets = kegg.gs,ref=c(1:3),samp=c(4:6),compare ="paired")
gene_set_1<-sigGeneSet(gos,cutoff=0.1)
out<-(gene_set_1$greater)
write.csv(out,file="kegg_up_path.csv")
write.xlsx(as.data.frame(out), "kegg_up_path.xlsx",asTable = TRUE,row.names=T)
out<-(gene_set_1$less)
write.csv(out,file="kegg_down_goterm.csv")
write.xlsx(as.data.frame(out), "kegg_down_path.xlsx",asTable = TRUE,row.names=T)

#######use other pathways database http://cpdb.molgen.mpg.de/
y<-scan("CPDB_pathways_genes.tab",what="",sep="\n")
ndemo <- vector(mode="list")
for (k in 1:length(y)){
  x=strsplit(y[k],"\t")[[1]]
  ndemo[[x[1]]]<-strsplit(x[length(x)],",")[[1]]
  print(x[1])
}
idx <- grep("disease",names(ndemo))
mouse <- ndemo[-idx]
save(mouse,file="pathway.dat")

####train the data with Msigdb database
####http://software.broadinstitute.org/gsea/msigdb/index.jsp
install.packages("msigdbr")












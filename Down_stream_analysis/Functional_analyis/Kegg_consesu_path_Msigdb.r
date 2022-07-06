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



#################################################-BP############################
go.hs <- go.gsets(species="mouse")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]
################################################################################
fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs,same.dir=TRUE)
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
out_up.bp<-subset(fc.go.bp.p.up,  q.val< 0.5)
write.csv(out_up.bp,file="go_upregulated.csv")
################################################################################
fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
out_down.bp<-subset(fc.go.bp.p.up,  q.val< 0.5)
write.csv(out_down.bp,file="go_upregulated.csv")
################################################################################
############clusterProfilerh##################################################
######profiler
up<- resSig[(resSig$log2FoldChange > 0),]
down<- resSig[(resSig$log2FoldChange < 0),]

yy <-enrichGO(down$entrez[!is.na(down$entrez)], 'org.Mm.eg.db', ont="BP",
              pAdjustMethod = "BH",minGSSize = 15, maxGSSize = 500,
              pvalueCutoff  = 0.01)

pdf('BP_GOterm.pdf',width=10,height=8)
barplot(yy, drop = TRUE)
dev.off()
write.xlsx(data.frame(yy), "BP_GO_term_up.xlsx")
#####if you want to get with symbol name 
#ewp.up <- DOSE::setReadable(yy, org.Mm.eg.db,`keyType`  = "ENTREZID")
###########################################
yy1 <-enrichKEGG(down$entrez[!is.na(down$entrez)], 'mmu',
                  pAdjustMethod = "BH",
                pvalueCutoff  = 0.01)
pdf('KEGGterm.pdf',width=10,height=8) 
barplot(yy1, drop = TRUE)
dotplot(yy1)
write.xlsx(data.frame(yy1), "KEGGterm.xlsx")
dev.off()


p<-"21803|104111|18035|22340|18707|216869"

################################################################################
###############wiki pathways analysis 
# wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus", format = "gmt")
# wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
# wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
# wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
# wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
# wpid2gene
# wpid2name
# ewp <- enricher(resSig$gene[!is.na(resSig$gene)],
#                 TERM2GENE = wpid2gene, 
#                 TERM2NAME = wpid2name,pAdjustMethod = "BH",minGSSize = 10, maxGSSize = 500,qvalueCutoff = 0.05)
# 
# ewp.up <- DOSE::setReadable(ewp, org.Mm.eg.db,`keyType`  = "ENTREZID")
# head(ewp.up)
# dev.copy2pdf(file="wiki_pathwas_day6.pdf")
# write.xlsx(data.frame(ewp.up), "wiki_pathway_day6.xlsx")
##############################################################################

organism="org.Mm.eg.db"
gene_list = sort(foldchanges, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             minGSSize = 3, 
             maxGSSize = 200, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")


gseaKEGG <- gseKEGG(geneList = gene_list, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "mmu", # supported organisms listed below
                    nPerm = 200, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)





bp <- enrichGO(gene_list, ont="BP",keyType = "ENTREZID", OrgDb = 'org.Mm.eg.db')
bp <- pairwise_termsim(bp)

m_df <- msigdbr(species = "Mus musculus")

m_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% 
dplyr::select(gs_name, entrez_gene)
em <- enricher(down$entrez[!is.na(down$entrez)], TERM2GENE=m_t2g)
em

C3_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% 
dplyr::select(gs_name, entrez_gene)
em2 <- GSEA(gene_list, TERM2GENE = C3_t2g)

C5_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
em2 <- GSEA(gene_list, TERM2GENE = C5_t2g)


######creat gene2term
y<-scan("~/Desktop/count/CPDB_pathways_genes.tab",what="",sep="\n")
ndemo <- vector(mode="list")
for (k in 1:length(y)){
  x=strsplit(y[k],"\t")[[1]]
  ndemo[[x[1]]]<-strsplit(x[length(x)],",")[[1]]
  print(x[1])
}
idx <- grep("disease",names(ndemo))
mouse <- ndemo[-idx]
test<-melt(mouse)
colnames(test)<-NULL
g2t<-test[,c(2,1)]
em2 <- GSEA(gene_list, TERM2GENE = g2t)
################################
kg.mmu=kegg.gsets("mmu")
kg.mmu<-kegg.gsets(species = "mmu", id.type = "entrez")
kegg.gs=kg.mmu$kg.sets[kg.mmu$sigmet.idx]
kegg.gs=kg.mmu$kg.sets
idx <- grep("disease",names(kegg.gs))
kegg.gs<- kegg.gs[-idx]
k2g<-melt(kegg.gs)
k2g<-k2g[,c(2,1)]
em2 <- GSEA(gene_list, TERM2GENE = k2g)

###########
keep<-list()
for (tem in rownames(out)){
idx<-melt(as.vector(mget(mygo[["tem"]],org.Mm.egSYMBOL)))
keep[[tem]]<-idx
}

dat<-load("~/Desktop/count/full_data_norm.rda")
akt<-read.delim("Akt_path.txt",header = T,sep=",")
a<-data[data$symbol %in%akt$Pathway_Akt,]
rownames(a)<-a$symbol
a<-as.matrix(a[,c(1:12)])
pdf(file = "Akt_path.pdf",width = 5,height = 8.5)
pheatmap::pheatmap(a[,c(1:3,10:12,4:9)],
                   cluster_rows = TRUE,
                   cluster_cols = F,fontsize_row = 6,
                   #legend_breaks =c(-3,0,5,10,14,
                                    #max(a)),legend_labels = c("-3", "0", "5", "10","14", "Log2cpm\n"),
                   main = "AKT_Path")
dev.off()


wnt<-read.delim("Wnt_pathw.txt",header = T,sep=",")
a<-data[data$symbol %in%wnt$Wnt_Pathw,]
rownames(a)<-a$symbol
a<-as.matrix(a[,c(1:12)])
pdf(file = "wnt_path.pdf",width = 5,height = 8.5)
pheatmap::pheatmap(a[,c(1:3,10:12,4:9)],
                   cluster_rows = TRUE,
                   cluster_cols = F,fontsize_row = 6,
                   #legend_breaks =c(-3,0,5,10,14,
                   #max(a)),legend_labels = c("-3", "0", "5", "10","14", "Log2cpm\n"),
                   main = "Wnt_Path")
dev.off()

#####
idx <- grep("Akt",names(kegg.gs))
keg<- kegg.gs[idx]
akt<-as.character(mget(as.vector(unlist(keg)),org.Mm.egSYMBOL,ifnotfound=NA))
akt<-data.frame(noquote(akt))
idx <- grep("wnt",names(kegg.gs))
keg<- kegg.gs[idx]
wnt<-as.character(mget(as.vector(unlist(keg)),org.Mm.egSYMBOL,ifnotfound=NA))









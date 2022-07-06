
##############
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(patchwork)
library(RaceID)
library(openxlsx)
###################################################################
Load data in to R
####################################################################
BM <- CreateSeuratObject(counts =BM_data, project = "BoneMetastasis", min.cells =3,min.features = 200)
###################################################################
Filter specific cell marker 
###################################################################
BM<-subset( BM, subset = PECAM1>0.5)
BM <-SCTransform(BM)
############################################################
top10 <- head(VariableFeatures(BM), 10)
plot1 <- VariableFeaturePlot(BM)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.copy2pdf(file="variable_plot.pdf")

##########################################################################
Run PCA 
##################################################################
BM <- JackStraw(BM, num.replicate = 100)
JackStrawPlot(BM, dims = 1:50)
dev.copy2pdf(file="JackStra_plot.pdf")
ElbowPlot(BM)
dev.copy2pdf(file="elbow_plot.pdf")

BM  <- RunPCA(BM , features = VariableFeatures(object = BM ))
BM  <- RunUMAP(BM , dims = 1:25, verbose = F)
BM <- FindNeighbors(BM , dims = 1:25, verbose = F)
BM<- FindClusters(BM, resolution = 0.4)


#############################################################################################################
FIND MARKER GENES list and theier count matrix
#############################################################################################################
all.markers <- FindAllMarkers(BM, only.pos = TRUE, min.pct = 0.25, logfc.threshold =0.5,return.thresh = 0.01)
write.xlsx(data.frame(all.markers),"all_cluster_marker_gene.xls",rownames=T)
write.table(as.matrix(GetAssayData(object = BM, slot = "counts")), 'counts.csv', sep = ',', row.names = T, col.names = T, quote = F)

################################################################################################################
PLOT GENE OF INTEREST 
################################################################################################################
features = c( "PECAM1",FLT4")
FeaturePlot(BM, features = features,ncol = 2)
dev.copy2pdf(file="feature_plot.pdf")

VlnPlot(BM,features = features,ncol = 2)
dev.copy2pdf(file="VINPOT.pdf")

DimPlot(BM, reduction = "umap",label = T)
dev.copy2pdf(file="umap.pdf")

pdf(file="gene_umap.pdf",width = 12,height = 6)
p1<-FeaturePlot(object = BM, features,ncol = 2)
p2<-DimPlot(BM, label = T)
CombinePlots(plots = list(p1, p2))
dev.off()
##############################################################################################
Anotate the cell type based on the top marker gene plot them 
#################################################################################
BM <- RenameIdents(object = BM, `0` = "Endothelial cell", 
                                `1` = "Tcell", 
                                `2` = "Macrophages", 
                                `3` = "Luminal cell", 
                                `4` = "B cells",
                                `5` = "Macrophages/MICell")
                              
pdf(file="pca_plot_with_cell_type.pdf",width = 10,height = 10)
DimPlot(BM , reduction = "pca")
dev.off()

pdf(file="label_umap_plot.pdf",width = 10,height = 10)
DimPlot(object = BM , label = TRUE)
dev.off()

pdf(file="label_tsne_plot.pdf",width = 10,height = 10)
TSNEPlot(object = BM)
dev.off()


###########################################################################################
GO/ PATHWAYS ANALYSIS OF SPECIFIC CLUSTER
###########################################################################################

clut<-all.markers[which(all.markers$cluster == "6"),]
all.gen <- subset(clut, p_val_adj< 0.001)
all.gen$geneid <- as.character(mget(rownames(all.gen),org.Hs.egSYMBOL2EG,ifnotfound=NA))
mygo <- as.list(org.Hs.egGO2EG)
mygo <- (mygo[!is.na(mygo)])
t <- mget(names(mygo),GOTERM)
names(mygo) <- as.character(lapply(t,Term))
deseq2.fc=all.gen$avg_log2FC
names(deseq2.fc)=all.gen$geneid
gos <- gage(deseq2.fc,gsets = mygo,ref = NULL, samp = NULL)
gene_set<-sigGeneSet(gos,cutoff=0.01)
up_path<-(gene_set$greater)
down_path<-(gene_set$less)

#############################################################################################
Kegg pathways analysis 
##############################################################################################
kg.hsa=kegg.gsets("hsa")
kg.hsa<-kegg.gsets(species = "hsa", id.type = "entrez")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.gs=kg.hsa$kg.sets
gos <- gage(deseq2.fc,gsets=kegg.gs,ref = NULL, samp = NULL)
up<-(gene_set$greater)
write.csv(up,file="go_upregulated.csv")
down<-(gene_set$less)
write.csv(down,file="go_downregulated.csv")


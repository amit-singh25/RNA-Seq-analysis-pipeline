library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(reshape)
library(tibble)
library(tidyr)
library(dplyr)

ctrl.names <- paste("CTRL_",c(1:378), sep = "")
colnames(ctrl.data) <- ctrl.names

trt.names <- paste("TRT_",c(1:495), sep = "")
colnames(trt.data) <- trt.names

final <- merge(ctrl.data,trt.data, by="row.names", all = TRUE)
final[is.na(final)] <- 0

rownames(final) <- final$Row.names
final$Row.names <- NULL

barplot(colSums(final))
abline(h=1000, col = "red")
#################################################################################
###########qulity control check cell name if necessary remove from the analysis  
###############################################################################
colnames(final[,final[c("Dcn"),]>8])

# name2id <- function(x,id) {
#   ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
#   n <- c()
#   for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
#   }
#   n
# }

# CGenes <- c("Pcna","Mki67","Mir703","Gm44044","Gm22757","Gm4775","Gm17541","Gm8225","Gm8730","Ptma","Actb","Hsp90aa1","Hsp90ab1","Ppia")
# FGenes<- c("Malat1","Xist")

sc <- SCseq(final)
sc <- filterdata(sc,mintotal=500, CGenes = NULL, FGenes = NULL) 
dim(sc@ndata)
sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc, cln = 0, sat = T, clustnr = 30)
plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)
sc <- clustexp(sc,cln=6,sat=FALSE, clustnr = 30)
sc <- findoutliers(sc, probthr = 0.0001)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
sc <- comptsne(sc)
sc <- compfr(sc,knn=10)
plotmap(sc)
dev.copy2pdf(file='tsenclsuetring.plot.pdf')
plotmap(sc,fr=TRUE)
dev.copy2pdf(file='fr_clsutering.plot.pdf')
types <- sub("(\\_|\\.).+","", colnames(sc@ndata))
plotsymbolsmap(sc,types,fr=F,cex = 0.5)
dev.copy2pdf(file='recid_t_sne.plot.pdf')
###heatmpa cluster 
clustheatmap(sc)
dev.copy2pdf(file='cell_cell_corelation.pdf')
###endothelialcell
plotexpmap(sc,"Pecam1",logsc=TRUE,fr=F)
##tcell
plotexpmap(sc,"Cd3g",logsc=TRUE,fr=F)
##bcell
plotexpmap(sc,"Igkc",logsc=TRUE,fr=F)
####
plotexpmap(sc,"Cd74",logsc=TRUE,fr=F)
plotexpmap(sc,"Dcn",logsc=TRUE,fr=F)

#########extracting genes for each clusters
genes <- list()
dg.all <- list()
for (i in c(1:8)) {
  dg <- clustdiffgenes(sc,i,pvalue=.05)
  ##   dg <- dg[grep("Rpl|Rps|Rik|RP|Malat|Jun|Fos|Hsp|Actb|Eef1a1|Ptma", rownames(dg), invert = T),]
  dg <- dg[order(dg$fc, decreasing = T),]
  dg.all[[i]] <- dg
  genes[[i]] <- head(dg,25)
  genes[[i]]$cluster <- rep(i,25)
}
#genes <- unique(unlist(genes))
#final_clut<-ldp
#####save#########################
final_clut<-do.call(rbind.data.frame, genes)
write.xlsx(final_clut,"Race_id_all_cluster_marker.xls",row.names=T)

#plotmarkergenes(sc,genes,cl=c(1:8),samples=types,order.cells=F,aggr = T, cap = 3, cluster_cols = T)
###############################################
###diffentiall express gene in each cluster
###########cluster-1#########################
dg <- clustdiffgenes(sc,6,pvalue=.05)
dg <- dg[order(dg$fc, decreasing = T),]
head(dg,25)
plotmap(sc)
############plot gene biomarker in tsnse plot 
my_gene <- c("Igkc", "Cd3e", "Cd74", "Cd8a", "Pecam1", 
             "Fabp4", "Car3", "Ccl8", "Lyz2",
             "Trbc2", "Dcn", "H2-Ab11", "Folr2")


plotmap_df <- sc@cpart %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  # Bind columns since cells in sc@cpart and tsne should be ordered the same way
  bind_cols(sc@tsne) %>%
  set_colnames(c("cellid", "cluster", "X", "Y")) %>%
  mutate("is_medoid"= cellid %in% sc@medoids, 
         cluster = factor(cluster, levels = 1:length(cluster) ) ) 


# save gene expression data frame
genes_df <- sc@ndata %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% my_gene) %>%
  gather(key = "cellid", value="scaled_expr", -gene) %>%
  left_join(plotmap_df, by="cellid") %>%
  left_join( tibble(cellid=names(sc@counts), count=sc@counts), by="cellid" ) %>%
  # transform scaled expression into transcript counts
  mutate(expr = scaled_expr * min(count) + 0.1) %>%
  arrange(gene, expr)


# For multiple genes
g <- genes_df %>% 
  ggplot(aes(X, Y, color = log2(expr) )) +
  geom_point(size=0.3) + 
  scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100) ) +
  theme_void() +
  facet_wrap(facets = ~gene) +
  coord_fixed(ratio = 1) +
  theme(strip.text.x = element_text(size = 15, face = "bold.italic", colour = "black", angle = 0))
ggsave(filename = paste("ggplotexpmap_", "B_pDC_genes", ".pdf", sep=""), plot = g, device = "pdf")

#types <- sub("(\\_|\\.).+","", colnames(sc@ndata))
#genes <- head(rownames(dg)[dg$fc>1],10)
#plotmarkergenes(sc,genes,samples=types)
#barplot(colSums(sc@ndata))
#############trajectory plot 

ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=15,nmode=T,fr=F)
#ltr <- projback(ltr,pdishuf=500)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)
plotgraph(ltr,scthr=0.6,showCells=FALSE,showTsne=TRUE)
x <- compscore(ltr,scthr=0.6)



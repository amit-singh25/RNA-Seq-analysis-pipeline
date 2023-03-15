####PCA plot with loading arrow.

mymat<-assay(vsd)
p5 <- labdsv::pca(t(mymat),dim=3,cor=F)
mysum <- summary(p5)
plot(p5$score[,1:2],pch=19,col=c(rep("green",12),rep("red",11)),
     xlim=c(-25,40),ylim=c(-25,35),type="n",
     ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))
loads <- as.data.frame(p5$loadings[,1:2])
loads$sym<-mapIds(org.Hs.eg.db, keys=row.names(loads), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
loads$len <- sqrt(loads[,1]^2+loads[,2]^2)
loads <- loads[order(abs(loads[,4]),decreasing=T),]
num <- 4
arrows(0,0,200*loads[1:num,1],200*loads[1:num,2],col="blue",length=0.)
text(200*loads[1:num,1],200*loads[1:num,2],loads[1:num,3],cex=0.5,col="blue")
text(p5$score[,1],p5$score[,2],rownames(p5$score),pos=3)
#test<-loads[1:500,]

#lines(p5$score[1:9,1],p5$score[1:9,2],col="red",type="o",pch=19)
#lines(p5$score[11:20,1],p5$score[11:20,2],col="green",type="o",pch=19)
#legend("topleft",c("Control","Treatment"),fill=c("red","green"),inset=0.02)
#dev.off()


###############distsample ##############################################
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL

###############sampleDists##################################################
pdf(file = "sampleDists.pdf",width = 10, height = 8)
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors)
dev.off()  

####correlation plot ######################################################
pdf(file = "Samplecor.pdf",width = 10, height = 8)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat, method = "pearson")
as.dist(1-vsd_cor, upper = TRUE) %>%
as.matrix %>%
pheatmap::pheatmap(., main = "Pearson correlation",col=colors)
#pheatmap(vsd_cor,display_numbers = TRUE, col=colors)
#pheatmap(vsd_cor, col=colors)
dev.off() 

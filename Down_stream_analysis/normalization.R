library(affy)
library(limma)
library(mogene10sttranscriptcluster.db)
library(sva)
library(genefilter)
#Normalization
#eset.diabetes=justRMA(celfile.path="../raw/Diabetes/Gene expression Diabetes Type 1 and 2")
#pData(eset.diabetes)$type=as.factor(c(rep("Db",5),rep("STZ",4),rep("WT",9)))


#normalizing the Db data
eset.diabetes.Db=justRMA(celfile.path="../raw/Diabetes/Gene expression Diabetes Type 1 and 2",filenames=c("Db_2106.CEL","Db_2107.CEL","Db_2108.CEL","Db_2109.CEL","Db_2110.CEL","WT_2112.CEL","WT_2113.CEL","WT_2115.CEL","WT_2519.CEL","WT_2520.CEL"))
pData(eset.diabetes.Db)$type=as.factor(c(rep("Db",5),rep("WT",5)))

#Output

formated=exprs(eset.diabetes.Db)
formated=as.data.frame(formated)
formated$symbol=as.vector(unlist(mget(rownames(formated),mogene10sttranscriptclusterSYMBOL,ifnotfound=NA)))
formated=formated[!is.na(formated$symbol),]
myiqr=apply(formated[,1:10],1,IQR)
fL=findLargest(rownames(formated),myiqr,"mogene10sttranscriptcluster")
exprs.eset <- formated[fL,]
formated=formated[fL,]
myiqr=apply(formated[,1:10],1,IQR)
formated=formated[order(myiqr,decreasing=T),]

write.table(formated,file=paste(format(Sys.time(), "%Y%m%d"),"normalized_log_expressions_gloms_diabetes_Db.txt",sep="_"),sep="\t",col.names=NA)
save(eset.diabetes.Db,file=paste(format(Sys.time(), "%Y%m%d"),"eset_gloms_Db.dat",sep="_"))
system("rm eset_gloms_Db.dat")
system(paste("ln -s", paste(format(Sys.time(), "%Y%m%d"),"eset_gloms_Db.dat",sep="_"), "eset_gloms_Db.dat"))

pca=prcomp(t(formated[,1:10]))
score.pca=predict(pca)
pdf(paste(format(Sys.time(), "%Y%m%d"),"PCA_gloms_diabetes_Db.pdf",sep="_"))
plot(score.pca[,1],score.pca[,2],col=c(2,2,2,2,2,1,1,1,1,1),xlab="First principal component (40% of variability)", ylab="Second principal component (17% of variability)",main="PCA of the Diabetes type II glomeruli transcriptome")
legend(-20,20,c("WT","Db/Db, Diabetes type II"),col=c(1,2),pch=1)
dev.off()

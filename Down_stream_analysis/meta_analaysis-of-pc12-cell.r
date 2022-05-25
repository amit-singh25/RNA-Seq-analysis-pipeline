#########################
### Approx. with a 5th order poly
# anova pvalue 


#data<-read.table("Adjusted_expression.xls",header=T,sep="\t")
#data = data[c(-117,-11671),]
#data=data[c(51,1:22)]

exp2<-exp2<-read.delim("ngf_cntrl.csv",sep=",",row.names=1)

tpts_ngf_full <- c(0,0.5,1,2,3,4,5,6,8,12,12,12,24,48)
mypoly.ngf <- poly(tpts_ngf_full, 5, raw=TRUE)
# for control
tpts_ngf <- c(0,0.5,1,2,3,4,5,6,8,12,24,48)
tpts_ctl<-c(0,2,4,6,8,12,24,48)

mypoly.ctl <- poly(tpts_ctl, 5, raw=TRUE)
ctl.tpts_df <- data.frame(x=tpts_ngf)

ngf.approx <- as.data.frame(mat.or.vec(nrow(exp2),12))
ctl.approx <- as.data.frame(mat.or.vec(nrow(exp2),12))

for(i in 1:nrow(exp2)){
	print(i)
	# Ctl
	y <- as.numeric(exp2[i,15:22])
	pred <- predict.lm(lm(y ~ mypoly.ctl))
	ctl.approx[i,1:12] <- as.numeric(approx(tpts_ctl,pred,xout=tpts_ngf)$y)
	# NGF
	y <- as.numeric(exp2[i,1:14])
	pred <- predict.lm(lm(y ~ mypoly.ngf))
	ngf.approx[i,1:12] <- pred[c(1:10,13,14)]
}

names(ngf.approx) <- tpts_ngf
names(ctl.approx) <- tpts_ngf
rownames(ngf.approx) <- rownames(exp2)
rownames(ctl.approx) <- rownames(exp2)


########################
## PCA
library(org.Rn.eg.db)
library(labdsv)
 library(Cairo)
  CairoFonts(
             regular="Arial:style=Regular",
             bold="Arial:style=Bold",
             italic="Arial:style=Italic",
             bolditalic="Arial:style=Bold Italic,BoldItalic",
             symbol="Symbol"
     )

mymat <- cbind(ctl.approx[,2:11],ngf.approx[,2:11])
mysym <- mget(rownames(mymat),org.Rn.egSYMBOL,ifnotfound=NA)

p5 <- labdsv::pca(t(mymat),dim=4,cor=F)
mysum <- summary(p5)

CairoPDF(file="PC12_PCA_biplot.pdf",width=10/1.54,height=10/1.54)
plot(p5$score[,1:2],pch=19,col=c(rep("green",12),rep("red",11)),
xlim=c(-25,40),ylim=c(-25,35),type="n",
ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))

loads <- as.data.frame(p5$loadings[,1:2])
loads$sym <- as.character(mget(rownames(loads),org.Rn.egSYMBOL,ifnotfound=NA))
loads$len <- sqrt(loads[,1]^2+loads[,2]^2)
loads <- loads[order(abs(loads[,4]),decreasing=T),]
num <- 40
arrows(0,0,200*loads[1:num,1],200*loads[1:num,2],col="grey",length=0.)
text(200*loads[1:num,1],200*loads[1:num,2],loads[1:num,3],cex=0.5,col="grey")

text(p5$score[,1],p5$score[,2],rownames(p5$score),pos=3)
lines(p5$score[1:10,1],p5$score[1:10,2],col="red",type="o",pch=19)
lines(p5$score[11:20,1],p5$score[11:20,2],col="green",type="o",pch=19)
legend("topleft",c("Control","+NGF"),fill=c("red","green"),inset=0.02)

dev.off()

##############
# Gene Sets
mart <- read.csv("mart_export.txt",header=T,sep=",")
out <- setdiff(mart[,3],rownames(exp2))
idx <- which(mart[,3] %in% out)
mart <- mart[-idx,]

 load("~/Downloads/gsea/Consensus.RData")

  idx <- grep("disease",names(cons2))
 cons2 <- cons2[-idx]
 idx <- lapply(cons2,length)
  idx2 <- which(idx < 5)
  cons2 <- cons2[-idx2]

  # Convert from human to mouse using GeneAnswers
y_rat <- list()
for(i in 1:length(cons2)){
     print(i)
     idx <- which(mart[,1] %in% cons2[[i]])
     y_rat[[i]] <- mart[idx,3]
 }
 names(y_rat) <- names(cons2)
#### 

# Transcription Factors
x <- scan("~/Downloads/gsea/c3.tft.v3.1.entrez.gmt", what="", sep="\n")
y <- strsplit(x, "[[:space:]]+")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)
for(i in 1:length(y)){ y[[i]] <- y[[i]][-1]}
# Convert to rat
y_rat_tf <- list()
for(i in 1:length(y)){
     print(i)
     idx <- which(mart[,1] %in% y[[i]])
     y_rat_tf[[i]] <- mart[idx,3]
 }
 names(y_rat_tf) <- names(y)
idx <- grep("MIR",names(y_rat_tf))
y_rat_tf <- y_rat_tf[-idx]
idx <- grep("UNKNOWN",names(y_rat_tf))
y_rat_tf <- y_rat_tf[-idx]
 # 
 # Gage
 mymat <- cbind(ctl.approx[,2:12],ngf.approx[,2:12])

gos <- gage(mymat,gsets = y_rat,ref=1:11,same.dir=TRUE,compare = "paired")
significant.groups=sigGeneSet(gos,cutoff=0.001)
out <- -log10(as.matrix(significant.groups$greater[,6:16]))
out[out > 3] <- 3
#rownames(out) <- substr(rownames(out),1,60)
pdf("pc12_cons_up.pdf",height=20,width=7)
pheatmap(out,cellwidth=8,cellheight=7,cluster_col=F,cex=0.8,cluster_col_method="maximum")
dev.off()
out <- -log10(as.matrix(significant.groups$less[,6:16]))
out[out > 3] <- 3
pdf("pc12_cons_down.pdf",height=20,width=7)
pheatmap(out,cellwidth=8,cellheight=7,cluster_col=F,cex=0.8,cluster_col_method="maximum")
dev.off()
#####





####plot all regulated gene in one file


CairoPDF(file="all_rgualted_gene_plot.pdf",width=10/1.54,height=10/1.54)
xaxis=c(0,0.5,1,2,3,4,5,6,8,12,24)
par(mfrow=c(3,3))
limits =  c(-2,10)

for(i in 1:17)
{
	ngf.points=NGF[i,1:11]
	con.points=NGF_con[i,1:11]
	inhibitor.points=NGF_inhi[i,1:11]
	
	plot(xaxis,ngf.points,type="l",ylab="Gene Expression Fold (log2)",xlab="Time (hours)",xlim=c(0,25),xaxs = "i",main=NGF$symbol[i],ylim=limits)
	points(xaxis,ngf.points,pch=18,col="blue",lwd = 3)
	lines(xaxis,ngf.points,col="blue",legened=c("NGF","NGF+UO126","Control"))
	
	lines(xaxis,inhibitor.points,col="red")
	points(xaxis,inhibitor.points,pch=18,col="red",lwd = 3)
	
	lines(xaxis,con.points,col="black")
	points(xaxis,con.points,pch=18,col="black",lwd = 3)
	legend(c("NGF","NGF+UO126","Control"),fill=c("blue","red","black"))	
	
}

dev.off()


###top ranked gene 
a<-read.delim("eset.normed_t0.ngf_48h_mds.csv",header=T,sep=",")
b<-read.delim("gene.txt",header=T,sep=",")
ngf<-a[a$symbol %in% b$gene,]
ngf<-ngf[c(2:12,15)]
CairoPDF(file="all_upregulated_gene.pdf",width=22/1.54,height=12/1.54)
xaxis=c(0,0.5,1,2,3,4,5,6,8,12,24)
par(mfrow=c(3,4))
limits =  c(-1,8)
for(i in 1:nrow(ngf)){
ngf.points=ngf[i,1:11]
	
barplot(as.matrix(ngf[i,1:11]),main=ngf$symbol[i],ylim=limits,names.arg=c("0","0.5","1","2","3","4","5","6","8","12","24"),titile="up_regulated.pdf",col = "gray")
	
#barplot(as.matrix(ngf[i,1:11]),ylab="Gene Expression Fold (log2)",xlab="Time (hours)",main=ngf$symbol[i],ylim=limits,names.arg=c("0h","0.5h","1h","2h","3h","4h","5h","6h","8h","12h","24h"),titile="up_regulated.pdf",col = "gray")
}
dev.off()








#######all plot 

library(colorRamps)
cols <-  matlab.like2(1000)

a<-read.delim("eset.normed_t0.ngf_48h_mds.csv",header=T,sep=",")
ngf<-a[c(2:12,15)]

pCairoPDF(file="ngf_timeseries_top_color.pdf",width=10/1.54,height=10/1.54)

tpts_ngf <- c(0,0.5,1,2,3,4,5,6,8,12,24,48)
plot(tpts_ngf,ngf[1,1:11],main='',xlab='Time [h]',
type="n",pch=19,ylab="FE",ylim=c(-3,9),xlim=c(0,24),lwd=3)

for(i in 1:1000){
# Plot points and curves with the respective colors
	points(tpts_ngf, ngf[i,1:11],col=cols[i])
	curve(splinefun(tpts_ngf, ngf[i,1:11], method="monoH.FC")(x),
		  add=TRUE, lwd=0.5,col=cols[i], n=241)
	
}


########all_plot_with_black color

a<-read.delim("eset.normed_t0.ngf_48h_mds.csv",header=T,sep=",")
ngf<-a[c(2:12,15)]
CairoPDF(file="ngf_timeseries_top1000.pdf",width=10/1.54,height=10/1.54)
tpts_ngf <- c(0,0.5,1,2,3,4,5,6,8,12,24)
plot(tpts_ngf,ngf[1,1:11],main='',xlab='Time[h]',
ylab="Fold Expression(Log2)",ylim=c(-3,9),xlim=c(0,25), type="n")
for(i in 1:1000){
ngf.points=ngf[i,1:11]
points(tpts_ngf, ngf[i,1:11])
lines(tpts_ngf,ngf.points, type="l",lwd=0.5,lty=2)}

dev.off()

####3 data set at bar plot 

CairoPDF(file="seconadary_gene.pdf",width=10/1.54,height=10/1.54)
xaxis=c(0,0.5,1,2,3,4,5,6,8,12,24)
par(mfrow=c(3,3))
limits =  c(-3,8)
for(i in 1:nrow(all)){
ngf.po=all[i,c(2:12,16:26,29:39)]	
barplot(rbind(as.matrix(all[i,2:12]),as.matrix(all[i,17:27]),as.matrix(all[i,31:41])),ylab="Gene Expression Fold (log2)",xlab="Time (hours)",main=all$symbol[i],ylim=limits,col=c("darkblue","red","darkgreen"),beside=TRUE,names.arg=c("0","0.5","1","2","3","4","5","6","8","12","24"))
legend("topright", c("NGF", "EGF","UO"), fill=c("darkblue", "red","darkgreen"))
}
dev.off()


######plot line plot for 3 different 
a<-read.delim("ngf_ieg.csv",header=T,sep=",")
b<-read.delim("ngf_uo_ieg.csv",header=T,sep=",")
c<-read.delim("egf_ieg.csv",header=T,sep=",")
pdf(file='all_IEG.pdf', height=10, width=12)
xaxis=c(0,0.5,1,2,3,4,5,6,8,12,24)
par(mfrow=c(3,3))
limits =  c(-3,8)
for(i in 1:nrow(e)){
	ngf.points=e[i,2:12]
	inhibitor.points=e[i,13:23]
	ngf.points=e[i,24:34]
	
#barplot(as.matrix(ngf[i,1:11]),ylab="Gene Expression Fold (log2)",xlab="Time (hours)",main=ngf$anv_padj[i],ylim=limits,label=sprintf("t==%.1e",ngf$anv_padj),names.arg=c("0","0.5","1","2","3","4","5","6","8","12","24"))
#barplot(as.matrix(ngf[i,2:11]),ylab="Gene Expression Fold (log2)",xlab="Time (hours)",main=ngf$symbol[i],ylim=limits,names.arg=c("0","0.5","1","2","3","4","5","6","8","12","24"),titile="EGF_IEG.pdf")
plot(xaxis,ngf.points,type="l",ylab="Gene Expression Fold (log2)",xlab="Time (hours)",xlim=c(0,25),main=a$mysym[i],ylim=limits)))
points(xaxis,ngf.points,pch=18,col="blue")
lines(xaxis,inhibitor.points,col="red")
points(xaxis,inhibitor.points,pch=18,col="red")
lines(xaxis,egf.points,col="black")
points(xaxis,egf.points,pch=18,col="black")
}

dev.off()


#####plot after MEK inhbition 



a<-read.delim("eset.normed_t0.ngf_48h_mds.csv",header=T,sep=",")
b<-read.delim("eset.inhibitor_48h.csv",header=T,sep=",")
d<-cbind(a,b)
e<-d[c(2:12,18:28,15)]
ngf<-read.delim("Inhi_gene.txt",header=T,sep=",")
exp<-e[e$symbol %in% ngf$gene,]
#pdf(file="NGF_INHI.pdf",onefile=T)
par(lwd = 3) 
CairoPDF(file="NGF_INHI.pdf",width=22/2,height=10/1.5)
xaxis=c(0,0.5,1,2,3,4,5,6,8,12,24)
par(mfrow=c(3,3))
limits =  c(-1,8)
for(i in 1:nrow(exp)){
	ngf.points=exp[i,1:11]
	inhibitor.points=exp[i,12:22]
	barplot(rbind(as.matrix(exp[i,1:11]),as.matrix(exp[i,12:22])),main=exp$symbol[i],ylim=limits,col=c("gray25","gray"), las=0.3,beside=TRUE,names.arg=c("0","0.5","1","2","3","4","5","6","8","12","24"))
	}
legend("topright", c("NGF","UO"), fill=c("gray8","gray87"))

dev.off()












#### boolean model 
library("BoolNet")
library(colorRamps)
library(pheatmap)
library(labdsv)
library(Cairo)
CairoFonts(
regular="Arial:style=Regular",
bold="Arial:style=Bold",
italic="Arial:style=Italic",
bolditalic="Arial:style=Bold Italic,BoldItalic",
symbol="Symbol"
)

ngf<-loadNetwork("ngf_bool_time2013.txt")
mycol <- matlab.like(13)
mycol[1] <- "#FFFFFF"
myoder<-c(1:13,60,20,21,16:19,14,15,22:37,61,38:53,62,54:57,58,59)
mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))



### 1st-hour-on 
CairoPDF(file="MAPK_on",width=10/1.54,height=10/1.54)
path1<-getPathToAttractor(ngf,c(1,rep(0,61)))
out_path1<-path1[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8,cluster_col_method="maximum")
dev.off()

#dev.copy2pdf(file="0-1hr.pdf")

### SERPINE1 on

CairoPDF(file="serpine_on.pdf",width=10/1.54,height=10/1.54)
net1=path1[nrow(path1),] 
net1$SERPINE1<-1
path2<-getPathToAttractor(ngf,as.numeric(net1))
out_path2<-path2[,myoder]
pheatmap(t(out_path2)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
#dev.copy2pdf(file="il6_on_1-3hr.pdf")
dev.off()

### NPY_on 
CairoPDF(file="NPY_on.pdf",width=10/1.54,height=10/1.54)
net2=path2[nrow(path2),]
net2$NPY<-1
net2$SERPINE1<-1
path3<-getPathToAttractor(ngf,as.numeric(net2))
out_path3<-path3[,myoder]
pheatmap(t(out_path3)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()

### TNFRSF12a_on 

CairoPDF(file="TNFRSFA_on.pdf",width=10/1.54,height=10/1.54)
net3<-path3[nrow(path3),]
net3$NPY<-1
net3$SERPINE1<-1
net3$TNFRSF12A<-1
path4<-getPathToAttractor(ngf,as.numeric(net3))
out_path4<-path4[,myoder]
pheatmap(t(out_path4)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()


#dev.copy2pdf(file="npy_on.pdf")
#mycols <- c(rep(1,2),rep(2,1),rep(3,1),rep(4,1),rep(5,1),rep(6,1),rep(7,1),rep(8,1),rep(10,1),rep(11,3),rep(12,1))
#getAttractors(ngf, type = "synchronous",method = "random",startStates = list(), genesON = c(), genesOFF = c(), canonical = TRUE, randomChainLength = 10000, avoidSelfLoops = TRUE, geneProbabilities = NULL,returnTable = TRUE) 
#attractors <- getAttractors(ngf, method="chosen",startStates=list(rep(1, length(ngf$genes))))



#######script for inhbition 
CairoPDF(file="MEK_inhbition.pdf",width=10/1.54,height=10/1.54)
ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
mycol <- matlab.like2(12)
mycol[1] <- "#FFFFFF"
myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"MEK",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,48)))
dummy1 = c(rep(1,nrow(path1)))
path1_wi_dummy <- cbind(path1, dummy1)
out_path1<-path1_wi_dummy[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)

dev.off()

####jnk inhbition 

CairoPDF(file="JNK_inhbition.pdf",width=10/1.54,height=10/1.54)
ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
mycol <- matlab.like2(12)
mycol[1] <- "#FFFFFF"
myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"JNK",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,48)))
dummy1 = c(rep(1,nrow(path1)))
path1_wi_dummy <- cbind(path1, dummy1)
out_path1<-path1_wi_dummy[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()


###pi3k inhibtion 


CairoPDF(file="PI3K_inhbition.pdf",width=10/1.54,height=10/1.54)
ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
mycol <- matlab.like2(12)
mycol[1] <- "#FFFFFF"
myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"PI3K",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,48)))
dummy1 = c(rep(1,nrow(path1)))
path1_wi_dummy <- cbind(path1, dummy1)
out_path1<-path1_wi_dummy[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()

###stat3 inhbition 

CairoPDF(file="STAT3_inhbition.pdf",width=10/1.54,height=10/1.54)
ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
mycol <- matlab.like2(12)
mycol[1] <- "#FFFFFF"
myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"STAT3",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,48)))
dummy1 = c(rep(1,nrow(path1)))
path1_wi_dummy <- cbind(path1, dummy1)
out_path1<-path1_wi_dummy[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)

dev.off()


######robust analysis 

ngf<-loadNetwork("ngf_bool_time2013.txt")
attrs <- getAttractors(ngf,canonical=TRUE, method="random",startStates=10000)
perturbationResults <- sapply(1:100, function(i)
{
	perturbedNet <- perturbNetwork(ngf, perturb="functions", method="bitflip")
	perturbedAttrs <- getAttractors(perturbedNet, method="random",startStates=10000,canonical=TRUE)
	attractorIndices <- sapply(attrs$attractors,function(attractor1)
							   {index <- which(sapply(perturbedAttrs$attractors, function(attractor2)
													  {
													  identical(attractor1, attractor2)
													  }))
							   if (length(index) == 0)
							   NA
							   else
							   index
							   })
	return(attractorIndices)
})

numOccurrences <- apply(perturbationResults[1:length(attrs$attractors),,drop=FALSE], 1,function(row)sum(!is.na(row)))
cat("Attractors in original network:\n")
print(attrs)
cat("Number of occurrences of the original attractors",
"in 1000 perturbed copies of the network:\n")
for (i in 1:length(attrs$attractors))
{
	cat("Attractor ",i,": ",numOccurrences[i],"\n",sep="")
}


###error bar plot 


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
##ERK/NGF
ngf_mean<-c(0.04,1.00,0.90,0.80,0.49,0.33,0.26,0.31,0.40,0.20)
uo_mean<-c(0.01,0.47,0.63,0.34,0.17,0.20,0.21,0.23,0.25,0.18)
y.sd<-(ngf_mean)
y1.sd<-(uo_mean)
CairoPDF(file="ngf_uo_inhbition.pdf",width=10/1.54,height=10/1.54)
yy <- matrix(c(ngf_mean,uo_mean),2,10,byrow=TRUE)
ee <- matrix(c(y.sd,y1.sd),2,10,byrow=TRUE)*1.96/10
#barx <- barplot(yy, beside=TRUE,col=c("gray25","darkgray"), ylim=c(0,1.5),axis.lty=1,names.arg=c("0","0.17","0.50","1","2","4","6","8","12","24"), xlab="Time[h]", ylab="Fold Expression",legend=c("NGF","UO"))

barx <- barplot(yy, beside=TRUE,col=c("gray25","darkgray84"), ylim=c(0,1.5),axis.lty=1,names.arg=c("0","0.17","0.50","1","2","4","6","8","12","24"),legend=c("NGF","UO"))

error.bar(barx,yy,ee)
dev.off()
####jnk/ngf
ngf_mean<-c(0.44,0.35,0.22,0.53,1.00,0.51,0.56,0.65,1.05,0.55)
sp_mean<-c(0.14,0.15,0.16,0.28,0.47,0.22,0.26,0.26,0.27,0.25)
y.sd<-(ngf_mean)
y1.sd<-(sp_mean)
CairoPDF(file="ngf_jnk_inhbition.pdf",width=10/1.54,height=10/1.54)
yy <- matrix(c(ngf_mean,uo_mean),2,10,byrow=TRUE)
ee <- matrix(c(y.sd,y1.sd),2,10,byrow=TRUE)*1.96/10
barx <- barplot(yy, beside=TRUE,col=c("gray25","darkgray87"), ylim=c(0,1.5),axis.lty=1,names.arg=c("0","0.17","0.50","1","2","4","6","8","12","24"),legend=c("NGF","SP"))
error.bar(barx,yy,ee)
dev.off()
######AKT/NGF



ngf_mean<-c(0.20,1.00,0.98,0.47,0.39,0.48,0.58,0.82,1.02,0.73)
nly_mean<-c(0.11,0.22,0.22,0.16,0.18,0.34,0.28,0.37,0.36,0.38)
y.sd<-(ngf_mean)
y1.sd<-(nly_mean)
CairoPDF(file="ngf_AKT_inhbition.pdf",width=10/1.54,height=10/1.54)
yy <- matrix(c(ngf_mean,uo_mean),2,10,byrow=TRUE)
ee <- matrix(c(y.sd,y1.sd),2,10,byrow=TRUE)*1.96/10
barx <- barplot(yy, beside=TRUE,col=c("gray25","gray87"), ylim=c(0,1.5),axis.lty=1,names.arg=c("0","0.17","0.50","1","2","4","6","8","12","24"), legend=c("NGF","NLY"))
error.bar(barx,yy,ee)
dev.off()
















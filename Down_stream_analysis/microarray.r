library(limma)
library(mouse4302.db)
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library(KEGG.db)
library(genefilter)

# Read in the Data
d <- read.table("2expression_data_RMA_unfiltered_annotated.txt",header=T,sep="\t",row.names=1)
# Delete all control Probesets
idx <- grep("AFFX",rownames(d))
d <- d[-idx,]

# Sort out non-entrez id
eid <- as.character(mget(rownames(d),mouse4302ENTREZID))
idx <- which(eid == "NA")
d <- d[-idx,]

# Delete Genes having the same Entrez ID according to the interquartile range
iqrs <- apply(d[,3:22],1,IQR)
uniqGenes <- findLargest(rownames(d),iqrs, "mouse4302")
d <- d[uniqGenes, ] 

# Delete the genes having a small IQR
thr <- quantile(d[,3],probs=c(0.1))
f1 <- pOverA(0.25, thr) 
f2 <- function(x) (IQR(x) > 0.1)
ff <- filterfun(f1,f2) 
selected <- genefilter(d[,3:22], ff)
d <- d[selected,]
 
# Name the Columns of the Chips accordingly
mynames <- c("WT 1","WT 2","WT 3","WT 4","WT 5"
,"RapKif 1"
,"RapKif 2"
,"RapKif 3"
,"RapKif 4"
,"RapKif 5"
,"Rap 1"
,"Rap 2"
,"Rap 3"
,"Rap 4"
,"Rap 5"
,"Kif 1"
,"Kif 2"
,"Kif 3"
,"Kif 4"
,"Kif 5")

names(d)[3:22] <- mynames

# Pick out the numeric data
mat <- d[,3:22]
names(mat) <- mynames

# Hierarchical Clustering of the data
 dd <- as.dist(1-cor(d[,3:22],method="pearson"))
 plot(hclust(dd),xlab="",ylab="1-correlation")


# PCA
mat.pca <- prcomp(t(mat), scale = TRUE)
mat.var <- mat.pca$sdev^2
mat.var.norm <- mat.var/sum(mat.var)
cumsum(mat.var)/sum(mat.var)
mat.pred <- predict(mat.pca)
 mycols <- rainbow(4)
 cols <- c(rep(1,5),rep(2,5),rep(3,5),rep(4,5))
 plot(mat.pred, type = "n", xlab = "PC1", ylab = "PC2",xlim=c(-150,140),ylim=c(-150,120))
 text(mat.pred[, 1],mat.pred[, 2], mynames,col=mycols[cols],cex=1.)

 # 3D PCA  
library(rgl)
plot3d(mat.pred[,1],mat.pred[,2],mat.pred[,3],cex=1,type="n",radius=2.5,col=mycols[cols])
text3d(mat.pred[,1],mat.pred[,2],mat.pred[,3]+5,mynames,col=mycols[cols])


 ########################################
 ######### 
 # Limma - get differentially regulated genes
 # Below we pick for different columns to look at different contrasts between
 # wild type, Rap KO, Kif KO, RapKif KO
 mymat <- mat[,c(16:20,8:10)]
 design <- model.matrix(~ 0+factor(c(rep(1,5),rep(2,3)))) 
  mymat <- mat[,c(1:5,8:10)]
 design <- model.matrix(~ 0+factor(c(rep(1,5),rep(2,3)))) 
 
 mymat <- mat[,c(1:5,11:15)]
 design <- model.matrix(~ 0+factor(c(rep(1,5),rep(2,5)))) 
colnames(design) <- c("WT", "Rap") 

mymat <- mat[,c(1:5,16:20)]
 design <- model.matrix(~ 0+factor(c(rep(1,5),rep(2,5)))) 
colnames(design) <- c("WT", "Kif") 

mymat <- mat[,c(8:10,11:15)]
 design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,5)))) 
colnames(design) <- c("RapKif", "Rap") 


mymat <- mat[,c(16:20,11:15)]
 design <- model.matrix(~ 0+factor(c(rep(1,5),rep(2,5)))) 
colnames(design) <- c("Kif", "Rap") 

# WT 
 fit <- lmFit(as.matrix(mymat), design)
 contrast.matrix <- makeContrasts(Kif-Rap, levels=design) 

 fit2 <- contrasts.fit(fit, contrast.matrix) 
 fit2 <- eBayes(fit2)
 res2 <- topTable(fit2, coef=1, adjust="BH",number=nrow(mymat),sort.by="logFC")
res2$Symbol  <- as.character(mget(res2$ID,mouse4302SYMBOL,ifnotfound=NA))
res2$Description  <- as.character(mget(res2$ID,mouse4302GENENAME,ifnotfound=NA))
res2 <- res2[order(res2$logFC,decreasing=TRUE),]
write.xls(res2,file="Kif_vs_Rap.xls")

#####################
 # Gage - this is a gene set enrichment
library(gage)
library(gplots)
library(GeneAnswers)

 # Read in the gene sets from the C2 data and convert the IDs from human to mouse
x <- scan("~/Downloads/gsea/c2.all.v3.0.entrez.gmt", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above
	for(i in 1:length(y)){
	 y[[i]] <- y[[i]][-1]
	}
# Convert from human to mouse using GeneAnswers
y_mouse <- list()
for(i in 1:length(y)){
	print(i)
     homoLL <-      getHomoGeneIDs(y[[i]], species='human', speciesL='mouse', mappingMethod='direct')
     y_mouse[[i]] <- as.character(homoLL)
}
names(y_mouse) <- names(y)

# Get all Kegg Genes and all GO genes
mykegg <- as.list(org.Mm.egPATH2EG)
mykegg <- mykegg[!is.na(mykegg)]
mygo <- as.list(org.Mm.egGO2EG)
mygo <- mygo[!is.na(mygo)]

# Create a matrix with Entrez IDs as row names
# Compare with type to Rap Kif
mymat <- mat[,c(1:5,6:10)]
rownames(mymat) <- as.character(mget(rownames(mymat),mouse4302ENTREZID))

# Use all GO Terms for Gage
gos <- gage(mymat,gsets = mygo,ref=c(1:5))
# Get the significant ones and plot the heatmap
siggo<-sigGeneSet(gos, outname="gse16873.kegg")
out <- siggo$greater
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))

# Plot the heatmap
#par(mar = c(8, 1, 1, 50))
heatmap.2(-log10(out),Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)
out <- siggo$less
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))
heatmap.2(-log10(out[,6:10]),Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)

#################################
# Gage for Rap versus Wildtype
mymat <- mat[,c(1:5,11:15)]
rownames(mymat) <- as.character(mget(rownames(mymat),mouse4302ENTREZID))

gos <- gage(mymat,gsets = mygo,ref=c(1:5))
# Get the significant Terms
siggo<-sigGeneSet(gos, outname="gse16873.kegg")
out <- siggo$greater
#rownames(out) <- as.character(getPathNames(rownames(out)))
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))

# Plot thre result as heatmap
#par(mar = c(8, 1, 1, 50))
myout <- -log10(out[,6:10])
myout[myout > 5] <- 5
pdf(file="RapDownGOTerms.pdf",width=10,height=20)
heatmap.2(myout,Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)
dev.off()

out <- siggo$less
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))
myout <- -log10(out[,6:10])
myout[myout > 5] <- 5
pdf(file="RapDownGOTerms.pdf",width=10,height=20)
heatmap.2(myout,Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)
dev.off()


#################################
# Kif versus Wildtype
mymat <- mat[,c(1:5,16:20)]
rownames(mymat) <- as.character(mget(rownames(mymat),mouse4302ENTREZID))

gos <- gage(mymat,gsets = mygo,ref=c(1:5))
siggo<-sigGeneSet(gos, outname="gse16873.kegg")
out <- siggo$greater
#rownames(out) <- as.character(getPathNames(rownames(out)))
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))

#par(mar = c(8, 1, 1, 50))
myout <- -log10(out[,6:10])
myout[myout > 5] <- 5
pdf(file="RapDownGOTerms.pdf",width=10,height=20)
heatmap.2(myout,Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)
dev.off()

out <- siggo$less
t <- mget(rownames(out),GOTERM)
rownames(out) <- as.character(lapply(t,Term))
myout <- -log10(out[,6:10])
myout[myout > 5] <- 5
pdf(file="RapDownGOTerms.pdf",width=10,height=20)
heatmap.2(myout,Rowv="none",dendrogram="col", margins=c(5,20),cex.lab=10)
dev.off()


###########################################################
# Linear Regression
#################################################################
# Lin regression with rap
# Choose the right columns for the regression, leave out the two bad RapKif Samples
mdr2 <- mat[,c(11:15,1:5,8:10,16:18)]
# X-axis values
xaxr <- c(rep(0,5),rep(1,5),rep(2,3),rep(3,3))
# Allocate some arrays
rs <- rep(0,nrow(mdr))
pv <- rep(0,nrow(mdr))
mdr2$Symbol <- as.character(mget(rownames(mdr),mouse4302SYMBOL))
mdr2$rs <- NA # R squared error
mdr2$pv <- NA  # p-value
mdr2$slope <- NA  # slope 
mdr2$intercept <- NA # intercepts

# Read in the limma analysis for Rap vs. Kif and choose only those that 
# are differentially regulated between Rap and Kif
kif <- read.table("Rap - Kif.csv",sep="|",header=T,row.names=2)

# Loop through all genes
for(i in 1:nrow(mdr)){
	print(i)
	# Linear regression
	fm <- lm(as.numeric(mdr2[i,1:16]) ~ xaxr)
	# Get output and calculated the output values
	out <- summary(fm)
	mdr2$rs[i] <- out$r.squared
	f.stat <- out$fstatistic
	mdr2$pv[i] <- 1-pf(f.stat["value"],f.stat["numdf"],f.stat["dendf"])
	mdr2$intercept[i] <-fm$coeff[1]
	mdr2$slope[i] <-fm$coeff[2]
}

##########################
## Only consider mdr
##########################
isna <- which(is.na(mdr2$rs))
mdr.o <- mdr2[-isna,]
# Take only those genes differentially regulated and upregulated in kif 
idx    <- which(kif$adj.P.Val < 0.05 & abs(kif$logFC) > 0.1)
idx <- which(rownames(mdr.o) %in% rownames(kif)[idx])

# Calc the standard deviation and mean of the R squared values
mysd   <- sd(mdr.o[idx,"rs"])
mymean <- mean(mdr.o[idx,"rs"])
# Take 2* standard defiation as cutoff
cutoff <- mymean + 1*mysd
idx <- which(mdr.o$rs > cutoff)
# If the rs is greater than this, take it as progress, otherwise as causal gene
mdr.o$progress <- 0
mdr.o$progress[idx] <- 1
# Add the connectivity data
cons <- read.table("cons_pkd.txt",sep="|",header=T,row.names=1,nrows=-1)
# Merge the arrays
mdr.c <- merge(mdr2,cons[,],by.x="row.names",by.y="row.names")
mdr.c <- mdr.c[,c(1:26,49)]
# Sortthe data
mdr.c <- mdr.c[order(mdr.c[,"conn"],decreasing=T),]
# Prepare the aoutput for protression genes
idx <- which(mdr.c$progress == 1)
prog2 <- mdr.c[idx,c(1,18:27)]
prog2 <- prog2[order(prog2$rs,abs(prog2$slope),decreasing=T),]
prog2$ensembl <- as.character(mget(as.character(prog2$Row.names),mouse4302ENSEMBL))
prog2$chr <- as.character(mget(as.character(prog2$Row.names),mouse4302CHR))
prog2$chrloc <- as.character(mget(as.character(prog2$Row.names),mouse4302CHRLOC))
prog2$chrlocend <- as.character(mget(as.character(prog2$Row.names),mouse4302CHRLOCEND))
write.xls(prog2[,c(1:6,11:15)],file="Progression_for_Anna2.xls")
#write.xls(prog[,c(1:6,11)],file="Progression_for_Anna.xls")

##################################
## Use the trap map and see how many regulators we find among the progression/
## The causal genes
trap <- read.delim("~/trap_network.txt",header=T,sep="\t")
trap <- trap[,1:3]
regs <- as.character(unique(trap[,1]))

idx <- which(mdr.o$progress == 1)
prog_symols <- as.character(mget(rownames(mdr.o)[idx],mouse4302SYMBOL))
idx <- which(prog_symols %in% regs)

idx2 <- which(trap[,1] %in% prog_symols[idx])
trap[idx2,]
#################################
# Read in the causal genes
causal <- read.delim("Causal_for_Anna2.xls",header=T,sep="\t")
idxc <- which(causal$Symbol.x %in% regs)

idx3 <- which(trap[,1] %in% causal$Symbol.x[idxc])
trap[idx3,]
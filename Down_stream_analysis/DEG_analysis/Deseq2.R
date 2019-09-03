rm(list=ls())
library(Cairo)
CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol")
library(DESeq2)
library(biomaRt)
library(openxlsx)

setwd("~/mouse_data/count")
#meta<-read.delim("meta.txt",header=T,sep="\t")
sampleFiles<- c("O1.htseq_count.txt",
                "O2.htseq_count.txt",
                "O4.htseq_count.txt",
                "P2.htseq_count.txt",
                "P4.htseq_count.txt",
                "A2.htseq_count.txt")

####name need to chnage A2 to P7 
#sampleNames <- c("A2","A5","A6","CT1","CT2","CT4","I1","12","I6","I7","O1","O2","04","P2","P4","P7")
sampleNames <- c("O1","O2","O4","P2","P4","P7")
#sampleCondition <- c("A","A","A","CT","CT","CT","I","I","I","I","O","O","O","P","P","P")
sampleCondition <- c("O","O","O","P","P","P")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = ".",design = ~ condition)

####PCA plot

Cairo(file="sample_PCA_plot",width=8,height=7,pointsize=12,type="pdf",units="in",dpi=96)
vsd <- varianceStabilizingTransformation(DESeq2Table, blind=TRUE)
Pvars <- rowVars(assay(vsd))
ntop = 1000
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,length(Pvars)))]
PCA <- prcomp(t(assay(vsd)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
PC3 = PCA$x[,3], PC4 = PCA$x[,4],
condition = colData(vsd)$condition)
(qplot(PC1, PC2, data = dataGG, color =  condition, size = I(1.5))+ggtitle("PCA plot")
+ theme(plot.title = element_text(hjust = 0.5))
+labs(x = paste0("PC1, ", round(percentVar[1],4),"%"),y = paste0("PC2, ", round(percentVar[2],4),"%"))
+xlim(c(-120,120)))

dev.off()
#####os simply use Deseq2
DESeq2::plotPCA(vsd, intgroup=c("condition"),ntop=1000)
dev.copy2pdf(file='PCA_plot.pdf')

#####Now compute DEG

ddsHTSeq <- DESeq(DESeq2Table)
##Comparison conditions between OvsLP
res <- results(ddsHTSeq, contrast=c("condition", "O", "P"))
##take adj.Pvalue 0.05 as significant DE
resSig <- subset(res, padj < 0.01)
##### Annotate Ensembl id to gene symbol,gene biotype and entrezgene by biomart"
resSig$ensembl <- sapply(strsplit(rownames(resSig), split="\\+" ), "[", 1 )
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id","mgi_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = resSig$ensembl,mart = ensembl)
idx <- match(resSig$ensembl, genemap$ensembl_gene_id)
resSig$name <- genemap$mgi_symbol[idx]
resSig$gene <- genemap$entrezgene[idx]
resSig$gene_biotype <- genemap$gene_biotype[idx]
resSig$symbol <- genemap$name_1006[idx]
#####Remove noncoding genes from the data
idx<-grep("protein_coding", resSig$gene_biotype, fixed = TRUE)
resSig<-resSig[idx,]
resOrdered <- resSig[order(-resSig$log2FoldChange ),]
write.xlsx(data.frame(resOrdered), "DEG_list.xlsx",asTable = TRUE,row.names=T)
####










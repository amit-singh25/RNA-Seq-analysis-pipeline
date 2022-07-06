library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "/Users/amit/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-current/Genes/genes.gtf"
FASTA <- FaFile("~/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/BowtieIndex/genome.fa")
#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="GRCm38.71", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
gr<-seqlevels(reducedGTF)[1:25]
reducedGTF<- keepSeqlevels(reducedGTF,gr)

#Open the fasta file
#FASTA <- FaFile(FASTAfile)
#open(FASTA)
#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")
write.table(output, file="GC_lengths.tsv", sep="\t")

#### other methods

#getGeneLengthAndGCContent(id, org, mode=c("biomart", "org.db"))

final<-read.delim("final_hela_cell.txt",header=T,sep="\t",row.names = 1)
sizeFactors.subset=colSums(final)
dds <-final[rowSums(final) > 10, ]
idx <- match (rownames (dds),rownames (output))
gc<-as.data.frame(output[idx,])

stopifnot(all(rownames(dds) == rownames(gc)))
stopifnot(colnames(dds) == names(sizeFactors.subset))

dds<- cqn(dds, lengths=gc$Length,  x=gc$GC, sizeFactors = sizeFactors.subset,verbose = TRUE)

par(mfrow=c(1,2))
cqnplot(cqn.subset, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(cqn.subset, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
dds <- as.data.frame(cqn.subset$y + cqn.subset$offset)

#######


colnames(dds)<-c("hela_cell_1","hela_cell2")
dds$mean<-rowMeans(dds)
fn<-ecdf(dds$mean)
dds$percnt<-fn(dds$mean)*100
dds<-as.data.frame(cbind(dds,gc))
dds$ensembl<-rownames(dds)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = dds$ensembl,mart = ensembl)
idx <- match(dds$ensembl, genemap$ensembl_gene_id)
dds$symbol <- genemap$hgnc_symbol[idx]
dds$entrezgene_id <- genemap$entrezgene[idx]

dds<- dds[order(dds$mean,decreasing=T),]


######

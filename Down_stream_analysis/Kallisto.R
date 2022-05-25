library(GenomicFeatures)
library(Rsamtools)
library(tximport)
library(readr)
library(DESeq2)
#setwd("~/Desktop/mouse_long_code/")
txdb <- makeTxDbFromGFF(file="~/Desktop/mouse_long_code/genecode_genome/gencode/gencode.vM16.long_noncoding_RNAs.gtf",format="gtf")
#saveDb(txdb, file="txdb.sqlite")
#loadDb("txdb.sqlite")
##extracting information from txdb
#g <- genes(txdb) # GRanges, just start to end, no exon/intron information
#tx <- transcripts(txdb) # GRanges, similar to genes()
#e <- exons(txdb) # GRanges for each exon
#ebg <- exonsBy(txdb, by="gene") # exons grouped in a GRangesList by gene
#ebt <- exonsBy(txdb, by="tx") # similar but by transcript
# then get the transcript sequence
#genome<-FaFile("~/Desktop/mouse_long_code/genecode_genome/gencode.vM16.lncRNA_transcripts.fa")
#txSeq <- extractTranscriptSeqs(genome, ebt)
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
a<- gsub("\\..*","",tx2gene[,1])
b<- gsub("\\..*","",tx2gene[,2])
c<-cbind(a,b)
colnames(c)=colnames(tx2gene)
tx2gene <- as.data.frame(c)
write.table(tx2gene,file="tx2gene.csv",sep=",",quote = F)

####
#a<-read.delim("gene_name.txt",sep="|",header = F)
#tx2gene<-a[,c(2,6)]
#colnames(tx2gene)<-c("TXNAME","GENEID")
#write.table(tx2gene,file="tx2gene.csv",sep=",",quote = F)

#######clean the data 
####make all tsv file to one folder here we fisrt cut first colum thed cut all samples that are in 5th cloumn tpm, 
#i have nine sample so first it take first line which is name then gene lnegth then tpm which is 5th cloumn  
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30,35,40,45 > transcript_tpms_all_samples.tsv
###then it write the folder name like below 
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
######
perl -ple 's/(ENSMUST[0-9.]+)\|(ENSMUSG[0-9.]+)\|.*?\t(\d+\t[0-9.]+\t[0-9.]+\t[0-9.]+)/$1\t$2\t$3/g' transcript_tpms_all_samples.tsv > new.tsv
perl -ple 's/target_id/target_id\tensembl_id/g' new.tsv > final.tsv
######get the gene name and and transcript id from fastafile from following comand
grep "^>" gencode.vM16.lncRNA_transcripts.fa | sed -e "s/>//g;s/ /\t/g" | cut -f 1,4,7 >tex.txt 

kallisto.dir<-"/Users/amit/Desktop/mouse_long_code/gencode_out"
samples<- c("ERR1141924","ERR1141925",	"ERR1141926",	"ERR1141927",	"ERR1141928",	"ERR1141929",
            "ERR1141930","ERR1141931","ERR1141932")
kallisto.files <- file.path(kallisto.dir,samples,  "abundance.tsv")
names(kallisto.files)<- c("ERR1141924","ERR1141925",	"ERR1141926",	"ERR1141927",	"ERR1141928",	"ERR1141929",
                          "ERR1141930","ERR1141931","ERR1141932")
all(file.exists(kallisto.files))
txi <- tximport(kallisto.files, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion = TRUE)
#txi.kallisto.tsv <- tximport(kallisto.files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
#txi <- tximport(kallisto.files, type="kallisto", txIn=TRUE, txOut=FALSE, tx2gene=TranscriptToGeneMap, countsFromAbundance="no")
txi.tx <- tximport(kallisto.files, type = "kallisto", txOut = TRUE, tx2gene = tx2gene,ignoreTxVersion = TRUE)
#condition <- factor(c("E","E","E","H","H","H","L","L","L"))

sampleTable <- data.frame(condition = factor(rep(c("E", "H","L"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
counts <- counts(dds,norm=TRUE)
counts <- as.data.frame(counts)
res <- results(dds, contrast=c("condition", "E", "L"))
resSig <- subset(res, padj < 0.05)


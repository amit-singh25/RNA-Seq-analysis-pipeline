setwd("/Users/amit/Desktop/mouse_data/count/final_plot/ovsp")

GeneCounts <- counts(DESeq2Table)
expr <- DGEList(counts=GeneCounts)
expr <- calcNormFactors(expr)
norm <- cpm(expr, log = TRUE)
norm$ensembl <- rownames(norm)
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id","mgi_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = norm$ensembl,mart = ensembl)
idx <- match(norm$ensembl, genemap$ensembl_gene_id)
norm$name <- genemap$mgi_symbol[idx]
norm$gene <- genemap$entrezgene[idx]
norm$symbol <- genemap$name_1006[idx]
write.xlsx(data.frame(norm),"norm_data_all_sample.xlsx",asTable = TRUE,row.names=T)
write.csv(data.frame(norm),file="norm_data_all_sample.csv")
###convert tds file to bedgraph file by using igvtools

igvtools tdftobedgraph GSM1686146_H3K9me3_wt_CCGTCC_102214_s.tdf GSM1686146_H3K9me3_wt_CCGTCC_102214_s.bedgraph
###then conver nc12 file to proper visulazion format
setwd("/Users/amit/Desktop/michael_chipseq/")
###for single 
track<-import("GSM1686146_H3K9me3_wt_CCGTCC_102214_s.bedgraph")
b<-c("I", "II" ,"III","IV", "V", "VI","VII","supercont12.8" ,"supercont12.9","supercont12.10","supercont12.11", "supercont12.12","supercont12.13","supercont12.14", "supercont12.15", "supercont12.16",
 "supercont12.17","supercont12.18","supercont12.19", "supercont12.20")
seqlevels(track)<-b
rtracklayer::export(track, "GSM1686146_H3K9me3_wt_CCGTCC_102214_s_clear.bedgraph")
H3K9me3<-coverage(track,weight = mcols(track)[[1]])
save(H3K9me3,file="H3K9me3.rda")


###
track<-import("GSM1686139_H3K27me3_wt_CGTACG_102214_s.bedgraph")
b<-c("I", "II" ,"III","IV", "V", "VI","VII","supercont12.8" ,"supercont12.9","supercont12.10","supercont12.11", "supercont12.12","supercont12.13","supercont12.14", "supercont12.15", "supercont12.16","supercont12.17","supercont12.18","supercont12.19", "supercont12.20")
seqlevels(track)<-b
rtracklayer::export(track, "GSM1686139_H3K27me3_wt_CGTACG_102214_s.bedgraph_clear.bedgraph")
H3K27me3<-coverage(track,weight = mcols(track)[[1]])
save(H3K27me3,file="H3K27me3.rda")

########load serplo2 data get alignment
bowtie2 -p 4 -x /Users/amit/Desktop/genome/neurospora_crassa_bowtie_index -sensitive-local -U C3K04ACXX_S2_35_sequence.fastq -S C3K04ACXX_S2_35_sequence.sam
samtools view -bS C3K04ACXX_S2_35_sequence.sam > C3K04ACXX_S2_35_sequence.bam
samtools sort C3K04ACXX_S2_35_sequence.bam >C3K04ACXX_S2_35_sequence_sort.bam
samtools index C3K04ACXX_S2_35_sequence_sort.bam

####
 htseq-count -f bam -m intersection-nonempty -r pos --type=exon --idattr=gene_id --stranded=reverse C3K04ACXX_S2_35_sequence_sort.bam /Users/amit/Desktop/genome/Neurospora_crassa.NC12.34.gtf >serplo2_35min.txt

gtf<-import("~/Desktop/genome/Neurospora_crassa.NC12.34.gtf")
###extract start codon
Neurospora.fasta<-FaFile("~/Desktop/genome/Neurospora_crassa.NC12.dna.toplevel.fa")
a<-read.delim("final_serpol2_data.txt",header=T,sep="\t",)
final<-a[ order( -a[,"Count"] ), ]
final<-final[1:150,]
rownames(final) <- c()

#a<-c("NCU02265","NCU03967")
gene<-gtf[((mcols(gtf)$gene_id %in% final$Gene_id))]

##gene1<-gene[((mcols(gene)$type =="type"))]
gene1<-gene[((mcols(gene)$type =="gene"))]
#bed<-flank(gene1,width = 1500,both = TRUE)
bed<-gene1 + 1500
names(bed)<-mcols(bed)[[5]]


#seq<-flank(gene,width = 1500,both = TRUE)
#seq<-flank(bed, 2)
#selected_gene<-getSeq(x = Neurospora.fasta, param = seq)
#names(selected_gene)<-mcols(seq)[[5]]
#writeFasta(selected_gene, file="sleceted_gene_seq.fa")
#writeFasta(selected_gene, file="sleceted_gene_seq.txt")



######plot with peak range object
#Rle.H3K27me3<-import("/Users/amit/Desktop/michael_chipseq/GSM1686139_H3K27me3_wt_CGTACG_102214_s.bedgraph_clear.bedgraph",as="Rle")
H3K27me3.Profiles<-H3K27me3
#H3K27me3.Profiles<-S4Vectors::runmean(Rle.H3K27me3,101,endrule ="constant")
H3K27me3.Profiles<-H3K27me3.Profiles[bed]
for (i in 1:length(H3K27me3.Profiles)){
    if(as.vector(strand(bed))[i]=="-"){
        H3K27me3.Profiles[[i]]<-rev(H3K27me3.Profiles[[i]])
    }
}

#Rle.H3K9me3<-import("/Users/amit/Desktop/michael_chipseq/  GSM1686146_H3K9me3_wt_CCGTCC_102214_s_clear.bedgraph",as="Rle")
H3K9me3.Profiles<-H3K9me3
#H3K9me3.Profiles<-S4Vectors::runmean(Rle.H3K9me3,101,endrule ="constant")
H3K9me3.Profiles<-H3K9me3.Profiles[bed]
for (i in 1:length(H3K9me3.Profiles)){
    if(as.vector(strand(bed))[i]=="-"){
        H3K9me3.Profiles[[i]]<-rev(H3K9me3.Profiles[[i]])
    }
}


####plot
#csp1.Profiles.smooth<-S4Vectors::runmean(csp1.Profiles,101,endrule ="constant")
setwd("/Users/amit/Desktop/michael_chipseq/images")
for(i in 1:length( H3K9me3.Profiles)) {
    png(file = paste(names(bed)[i],".png",sep=""))
    dev
#for(i in 1:length( H3K9me3.Profiles)) {
#    png(file = paste(names(bed)[i],".png",sep=""))
#    plot( H3K9me3.Profiles[[i]],type="l",xlab="Base",ylab="Reads",main=names(bed)[i],
#    axes=F,lwd=5,ylim=c(0,max(c( H3K9me3.Profiles[[i]],(H3K27me3.Profiles[[i]]*10)))))
#    axis(side=1,at=c("0","1500", length(H3K9me3.Profiles[[i]])-1500, length(H3K9me3.Profiles[[i]])),labels=c("-1500","TSS","ATG","+1500"))
#    axis(side=2)
#    #axis(side=4,at=c("0"),labels=c("0"))
#    lines(H3K27me3.Profiles[[i]]*10,type="l",col="red",lwd=3)
#    legend("topright",c("H3K9me3.Profiles","H3K27me3.Profiles \n 10x magnified"),fill=c("black","red"))
#
#    #textxy(df$dese2, df$qpcr, labs=df$gene_name, cex=1)
#
#    dev.off()
#}

#####
library( hwriter )
page <- openPage( "methylation_profile.html",
head = paste( sep="\n",
"<script>",
"   function set_image( name ) {",
"      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
"   }",
"</script>" ) )
cat(file=page,
'<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
hwrite(final, border=NULL, page=page,
onmouseover = sprintf( "set_image( '%s' );", final$Gene_id ) )
cat( file=page,
'</div></td><td style="vertical-align:top"><img id="plot" width="600px"></td></tr></table>' )
closePage(page)
browseURL( "methylation_profile.html" )

#####
<script type="text/javascript">
var row = document.getElementsByTagName('tr');
var main_url = 'http://fungidb.org/fungidb/app/record/gene/';
for(i = 2; i < row.length; i++)
{
    td_last = row[i].firstElementChild;
    gene_name = td_last.textContent
    row[i].firstElementChild.textContent = '';
    gene_url = main_url + gene_name;
    var a_link = document.createElement('a');
    a_link.href = gene_url;
    a_link.text = gene_name;
    row[i].lastElementChild.appendChild(a_link);
}
</script>
</body>


#for(i in 1:length(H3K9me3.Profiles)) {
#   png(file = paste(mcols(bed)[[4]][i],".png",sep=""))
#   plot(csp1.Profiles[[i]],type="l",xlab="Base",ylab="Reads",main=mcols(bed)[[4]][i])
#   if((length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i])<=length(csp1.Profiles[[i]]) &
#   (length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i])>0){
#       points(x=(length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i]),
#       y=csp1.Profiles[[i]][(length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i])],
#       col="red",pch=16)
#       lines(x=c((length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i]),
#       (length(csp1.Profiles[[i]])/2)-(mcols(bed)[[3]][i])),
#       y=c(-50000000,5000000000000),
#       col="red",lty=3)
#   }
#   #axis(side = 1,at = (length(csp1.Profiles.smooth[[i]])/2)-(mcols(bed)[[3]][i]),labels = "TSS")
#   dev.off()
#}



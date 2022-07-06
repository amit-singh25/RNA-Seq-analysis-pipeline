library("DNAshapeR")
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library('ShortRead')


gtf<-import("~/Desktop/genome/Neurospora_crassa.NC12.34.gtf",as="Rle")

D<-import("~/mnt/N.some_data/alignment_A/1_sortA.bigWig",as="Rle")
D1<-import("~/mnt/N.some_data/alignment_A/2_sortA.bigWig",as="Rle")
D2<-import("~/mnt/N.some_data/alignment_A/3_sortA.bigWig",as="Rle")
D3<-import("~/mnt/N.some_data/alignment_A/4_sortA.bigWig",as="Rle")
D4<-import("~/mnt/N.some_data/alignment_A/5_sortA.bigWig",as="Rle")
D5<-import("~/mnt/N.some_data/alignment_A/6_sortA.bigWig",as="Rle")
D6<-import("~/mnt/N.some_data/alignment_A/7_sortA.bigWig",as="Rle")
D7<-import("~/mnt/N.some_data/alignment_A/8_sortA.bigWig",as="Rle")
D8<-import("~/mnt/N.some_data/alignment_A/9_sortA.bigWig",as="Rle")
D9<-import("~/mnt/N.some_data/alignment_A/10_sortA.bigWig",as="Rle")

gr2=GRanges(seqname="VII",range=IRanges(start=3131440, end=c(3132000)))

D.Profiles<-D[gr2]
D1.Profiles<-D1[gr2]
D2.Profiles<-D2[gr2]
D3.Profiles<-D3[gr2]
D4.Profiles<-D4[gr2]
D5.Profiles<-D5[gr2]
D6.Profiles<-D6[gr2]
D7.Profiles<-D7[gr2]
D8.Profiles<-D8[gr2]
D9.Profiles<-D9[gr2]

final<-cbind(mean(D.Profiles[[1]]),mean(D1.Profiles[[1]]),
      mean(D2.Profiles[[1]]),mean(D3.Profiles[[1]]),
      mean(D4.Profiles[[1]]),mean(D5.Profiles[[1]]),
      mean(D6.Profiles[[1]]),mean(D7.Profiles[[1]]),
      mean(D8.Profiles[[1]]),mean(D9.Profiles[[1]]))
colnames(final)<-c("Dark","1","5","10","20","30","45","60","120","240")         
barplot(final, horiz = T, ylab = "Light Pluse (Min)",xlab = "Mean Covarge (FRQ) (VII:3131440-3132000)",names.arg=names(final),las=1)         
setwd("~/Desktop/")

#####

gr3=GRanges(seqname="VII",range=IRanges(start=3108000, end=c(3128000)))

D.Profiles<-D[gr3]
D1.Profiles<-D1[gr3]
D2.Profiles<-D2[gr3]
D3.Profiles<-D3[gr3]
D4.Profiles<-D4[gr3]
D5.Profiles<-D5[gr3]
D6.Profiles<-D6[gr3]
D7.Profiles<-D7[gr3]
D8.Profiles<-D8[gr3]
D9.Profiles<-D9[gr3]

final1<-cbind(mean(D.Profiles[[1]]),mean(D1.Profiles[[1]]),
             mean(D2.Profiles[[1]]),mean(D3.Profiles[[1]]),
             mean(D4.Profiles[[1]]),mean(D5.Profiles[[1]]),
             mean(D6.Profiles[[1]]),mean(D7.Profiles[[1]]),
             mean(D8.Profiles[[1]]),mean(D9.Profiles[[1]]))

norm<-final/final1
colnames(final1)<-c("Dark","1","5","10","20","30","45","60","120","240")    
barplot(norm, horiz = T, ylab = "Light Pluse (Min)",xlab="Normalized Mean Covarge (FRQ) (VII:3131440-3132000)",names.arg=names(norm),las=1,main="Normalized Region VII:[3108000-3128000]") 



###future check
write.table(gr2, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)
bedtools intersect -wa -a Neurospora_crassa.NC12.34.gtf -b foo.bed 

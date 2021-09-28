library(dplyr)
library(ggplus)
library(tidyr)
#####plot all ressignificant gene in one pdf file 
pdf('light_inducibale_gene.pdf')
p<- data.frame(colData(fullData),t(counts(fullData, normalized=TRUE)[rownames(ressig),]))[-58,]%>%
  gather( "gene", "ncount", -(1:5) ) %>%
  ggplot + geom_line(aes( x=time, y=log2(ncount+1), col=samplelabel)) +facet_wrap(~gene)
p<- p + xlab( "Time[h]") + ylab("Log2(count)+1") + theme_bw()
p<-p+geom_point(aes(x=time,y=log2(ncount+1),col=samplelabel))
facet_multiple( plot = p, facets = 'gene', ncol = 4, nrow = 4)
dev.off()

######plot single gene to check 
data.frame( colData(fullData),t(counts(fullData, normalized=TRUE)["NCU02713", , drop=FALSE]))[-58,]%>%
  gather( "gene", "ncount", -(1:5) )%>%
  ggplot +geom_line(aes(x=time, y=log2(ncount+1), col=samplelabel)) + facet_wrap(~gene)
p+ xlab("Time[h]") + ylab("Log2(count)+1") + theme_bw()
p<-
###### plot in one gene in one png file, removed one sample for fulldata
for( i in row.names(norm))
 {
  p<-data.frame(colData(fullData),t(counts(fullData, normalized=TRUE)[i, , drop=FALSE]) )[-58,] %>%
    gather( "gene", "ncount", -(1:5) )%>%
    ggplot +geom_line(aes(x=time, y=log2(ncount+1), col=samplelabel)) + facet_wrap(~gene)
  p<- p+ xlab("Time[h]") + ylab("Log2(count)+1")
  p<-p+theme_bw()+coord_fixed(ratio=.5)
  p<-p+geom_point(aes(x=time,y=log2(ncount+1),col=samplelabel))
  p<-p+ scale_x_continuous(breaks = c(0,2, 4, 6,8,10,12,14,16,18,20,22)) 
  p<-p+ ylim(0,16)
 p<-p+theme(axis.text = element_text (face="bold", size=14), axis.title=element_text(size=18,face="bold"))
 p<-p+theme(legend.title = element_text(size=16, face="bold"),legend.text = element_text(size = 16, face = "bold"))
 p<-p+theme(strip.text = element_text(size=25,face="bold"))
 #theme(legend.title = element_text(colour="blue", size=16, face="bold"))
ggsave(paste('~/Desktop/Light-dark/images/',i , '.png', sep=''), width=20, height=15, p)
  print(i)
}
### scale_y_continuous(limits = c(-0.3,0.3))
####single plot
p<-data.frame( colData(fullData),t(counts(fullData, normalized=TRUE)["NCU02356", , drop=FALSE]))[-58,]%>%
  gather( "gene", "ncount", -(1:5) )%>%
  ggplot +geom_line(aes(x=time, y=log2(ncount+1), col=samplelabel)) + facet_wrap(~gene)
p+ xlab("Time[h]") + ylab("Log2(count)+1") + theme_bw()


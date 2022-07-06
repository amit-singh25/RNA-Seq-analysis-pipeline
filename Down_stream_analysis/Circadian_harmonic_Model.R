####norm data of dl-wt 
#data<- log2(norm[,c(2:14)]+1)
data<-  log2(norm[,c(2:14)]+1)
#x<-colData(fullData)$time[c(2:14)]
x<-colData(fullData)$time[c(2:6,9:14)]
#x = seq( 0,, by=2)
 approx <- as.data.frame( mat.or.vec( nrow( data), 14))
 pval<- as.vector( c() )
 amplitute <- as.vector( c() )
  phase <- as.vector( c() )
 for( i in 1:nrow( data ) ){
 y <- as.numeric( data[i,1:11 ] )
fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
fit0 <- lm( y ~ x )

#fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
#fit0 <-glm.nb( y ~ x )
amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
#phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
phase[i]=-1*atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
model<-anova( fit0, fit1 )
#model<-anova( fit0, fit1, test="Chisq" )
#pval[i]<-model$`Pr(>F)`[2]
pval[i]<-model$`Pr(>F)`[2]
#pval[i]<-model[["Pr(Chi)"]][2]
approx <- as.data.frame( cbind( data, pval, amplitute, phase ) )
 approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
 #approx$p.adj_by<-p.adjust( approx$pval, method = "BY",n = length( approx$pval ) )
 #approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
 #approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )
 print(i)
 
 }
###
#sig<- subset(approx,p.adj< 0.1)

sig<-approx[ order(approx[,"p.adj"] ), ]
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="final_with_out_11_with_out_0_log.rda")

#### creat aplot with html output
sig<-sig[,c(25:30)]
page <- openPage( "ryth_gene_final.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
bgcolor1=matrix( c( '#FAEBD7'), nr=10, nc=6 )
bgcolor2=matrix( c( '#E0EEEE'), nr=793, nc=6 )
hwrite(sig, border=NULL, page=page, bgcolor=rbind( bgcolor1,bgcolor2),
onmouseover = sprintf( "set_image( '%s' );",sig$ens.id))
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="700px"></td></tr></table>' )
closePage(page)
browseURL( "ryth_gene_final.html" )

###### useing deseq2
dd <- fullData[,fullData$samplelabel=="dd-wt"]
dd$phase <- as.integer( dd$time) / 22 * 2*pi
#colData( dd)$condition<-as.factor( c( "0h","2h","4h", "6h","8h","10h","12h","14h","16h","18h","20h","22h") )
colData( dd)$samplelabel<-as.factor(( dd)$samplelabel)
design = model.matrix( ~ dd$condition*( dd$time + cos( dd$phase) + sin( dd$phase) ))
#design = model.matrix( ~ dd$samplelabel*( dd$time + cos( dd$phase) + sin( dd$phase) ))
#y<- DGEList( counts = assay( dd))
#design <- ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi )
#table(rowSums(y$counts==0)==9)
keep <- rowSums( cpm( y) > 0.5) >= 2
table( keep)
#keep
y<-y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors( y)
v <- voom( y, design, plot=TRUE)
efit <- eBayes( vfit)
#######polar plot 
df<-melt(dat, id=c("gene","phase"))
ggplot(data=df,aes(x=df$gene,y=df$pahse,gorp=varia))+
      geom_bar(stat="identity")+
  coord_polar()+
  scale_fill_brewer(palette="Greens")+xlab("")+ylab("")
library(ggplot2)
##########################
mat <- norm[,c(29:52)]
#####take avrergae 
newmat <- as.data.frame( mat.or.vec( nrow( mat), 12))
newmat[,1] = apply(mat[,c(1,13)],1,mean)
newmat[,2] = apply(mat[,c(2,14)],1,mean)
newmat[,3] = apply(mat[,c(3,15)],1,mean)
newmat[,4] = apply(mat[,c(4,16)],1,mean)
newmat[,5] = apply(mat[,c(5,17)],1,mean)
newmat[,6] = apply(mat[,c(6,18)],1,mean)
newmat[,7] = apply(mat[,c(7,19)],1,mean)
newmat[,8] = apply(mat[,c(8,20)],1,mean)
newmat[,9] = apply(mat[,c(9,21)],1,mean)
newmat[,10] = apply(mat[,c(10,22)],1,mean)
newmat[,11] = apply(mat[,c(11,23)],1,mean)
newmat[,12] = apply(mat[,c(12,24)],1,mean)
rownames(newmat)<-rownames(mat)
colnames(newmat)<-c("0hr","2hr","4hr","6hr","8hr","10hr","12hr","14hr","16hr","18hr","20hr","22hr")
save(newmat,file="dark_rep_mean.rda")
#newmat<-log2(newmat[,c(1:12)]+1)
x<-colData(fullData)$time[c(29:40)]
#data<-newmat[,c(1:12)]
##
approx <- as.data.frame( mat.or.vec( nrow( newmat), 15))
pval<- as.vector( c() )
amplitute <- as.vector( c() )
phase <- as.vector( c() )
for( i in 1:nrow( newmat ))  {
y <- as.numeric( newmat[i,1:12 ] )
fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
fit0 <- lm( y ~ x )
#fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
#fit0 <-glm.nb( y ~ x )
amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
#phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
phase[i]=atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
model<-anova( fit0, fit1 )
#model<-anova( fit0, fit1, test="Chisq" )
#pval[i]<-model$`Pr(>F)`[2]
pval[i]<-model$`Pr(>F)`[2]
#pval[i]<-model[["Pr(Chi)"]][2]
print(i)
}
approx <- as.data.frame( cbind(newmat, pval, amplitute, phase ) )
#approx$p.adj<-p.adjust( pval, method = "BH",n = length( pval ) )
approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
#approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )

#######
sig<- subset(approx,p.adj< 0.1)
sig<-approx[ order(approx[,"p.adj"] ), ]
dim(sig)
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="ryth_dd_rep_with_out_0_log.rda")


#######
index=(apply(newmat,1,sum)>20) 
row_sub = apply(newmat, 1, function(row) all(row !=0 ))
data<-newmat[index,]
approx <- as.data.frame( mat.or.vec( nrow(data), 15))
pval<- as.vector( c() )
amplitute <- as.vector( c() )
phase <- as.vector( c() )
x<-colData(fullData)$time[c(29:40)]
for( i in 1:nrow( data ))  {
  y <- as.numeric( data[i,1:12 ] )
  #fit1 <- lm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
  #fit0 <- lm( y ~ x )
fit1<-glm.nb( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ) )
fit0 <-glm.nb( y ~ x )
amplitute[i]<- sqrt( coef(fit1 )[3]^2 + coef( fit1 )[4]^2 )
#phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
phase[i]=atan2( coef( fit1 )[4], coef( fit1 )[3] )*22/( 2*pi )
model<-anova( fit0, fit1 )
model<-anova( fit0, fit1, test="Chisq" )
pval[i]<-model$`Pr(>F)`[2]  #pval[i]<-model$`Pr(>F)`[2]
  #pval[i]<-model[["Pr(Chi)"]][2]
#  approx <- as.data.frame( cbind( data, pval, amplitute, phase ) )
 # approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval ) )
  #approx$p.adj_by<-p.adjust( approx$pval, method = "BY",n = length( approx$pval ) )
  #approx$p.adj_homel<-p.adjust( approx$pval, method = "hommel",n = length( approx$pval ) )
  #approx$p.adj_homel<-p.adjust( approx$pval, method = "bonferroni",n = length( approx$pval ) )
  print(i)
}

#######
sig<- subset(approx,p.adj< 0.1)
sig<-approx[ order(approx[,"p.adj"] ), ]
keep<-list()
for (i in rownames( sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[ idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
#df<-df[!duplicated(df), ]
sig<-as.data.frame( cbind( sig,df[!duplicated(df), ]))
save(sig,file="ryth_dd_rep_with__out_0_log.rda")









  


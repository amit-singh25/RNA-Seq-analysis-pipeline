####norm data of dl-wt 

data<-  norm[,c(41:52)] 
data<- log2( norm[,c(41:52)]+1 )
x<-colData(fullData)$time[41:52]
#x = seq( 0,, by=2)
approx <- as.data.frame( mat.or.vec( nrow( data), 15))
pval<- as.vector(c())
amplitute <- as.vector(c())
phase <- as.vector(c())
for( i in 1:nrow( data) ){
  y <- as.numeric( data[ i,1:12 ])
  fit1 <-glm( y ~ x + cos( x/7.33*2*pi ) + sin( x/7.33*2*pi ))
  fit2<-glm( y ~ x + cos( x/11*2*pi ) + sin( x/11*2*pi ))
  fit3<-glm( y ~ x + cos( x/22*2*pi ) + sin( x/22*2*pi ))
  #fit0 <- glm( y ~ x )
  model<-anova(fit1,fit2,fit3,test="Chisq")
  pval[i]<-model$`Pr(>Chi)`[2]
  amplitute[i]<- sqrt( coef(fit1)[3]^2 + coef(fit1)[4]^2 )
  phase[i]<- atan2( coef(fit1)[4], coef( fit1)[3] )
  approx <- as.data.frame( cbind( data,pval,amplitute,phase) )
  approx$p.adj<-p.adjust( approx$pval, method = "BH",n = length( approx$pval))
  print(i)
}
###
#sig<- subset( approx,p.adj< 0.1)
#sig<-sig[ order(sig[,"p.adj"] ), ]
#save(sig,file="sig_12hr_adj.rda")
#########
sig<- subset( approx,pval< 0.1)
#sig<-sig[ order( sig[,"pval"] ), ]
#save(sig,file="sig_12hr_pval.rda")
keep<-list()
for (i in rownames(sig)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[idx,]
  keep[[i]]<-idx
} 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
#df<-df[!duplicated(df), ]
sig<-as.data.frame(cbind(sig,df[!duplicated(df), ]))

####
sig<-sig[,c(13:18)]
page <- openPage( "light_ind.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
#bgcolor1=matrix(c('#FAEBD7'),nr=510,nc=6)
#bgcolor2=matrix(c('#E0EEEE'),nr=285,nc=6)
hwrite(ressig, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );",ressig$ens.id))

cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="800px"></td></tr></table>' )
closePage(page)
browseURL( "light_ind.html" )

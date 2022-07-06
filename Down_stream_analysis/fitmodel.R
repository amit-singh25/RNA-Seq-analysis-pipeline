fitModel <- function(gene) {
  df <- subset(
    data.frame(
      colData( fullData ),
      ncount = counts( fullData, normalized=TRUE )[ gene, ] ),
    fullData$samplelabel == "dl-wt" )
  mu<- as.numeric(mean(df$ncount[3:6]))
  muD <- as.numeric(mean(df$ncount[11:14]))
  nls( log2(ncount+1) ~ ifelse( time > 11, mu + A*exp(-k*(time-11)), muD ),
       algorithm="port",
       #start = c( mu=300, A=2^5, k=0.5, muD=300 ),
       start = c( mu=mu, A=2^5, k=0.5, muD=muD ), 
       lower = c( mu=1, A=2^0, k=0.001, muD=0.001 ), 
       upper = c( mu=1e8, A=2^20, k=50, muD=1e8 ), 
       data=df )
}

getCurve <- function( fit ) {
  fittedline <- data.frame( time = seq( 0, 22, length.out=1000 ) )
  fittedline$ncount <- predict( fit, fittedline )
  fittedline
}
####extract fitted value for plot 
store<-list()
cof<-list()
for(i in rownames(up)){
  print(i)
  fit<-fitModel(i)
  line<- getCurve(fit)  
  store[[i]]<-line
  cof[[i]]<-coef(fit)
}
####get data frame for the html output
up <- up[ order( -up[,"log2FoldChange"] ), ]
final<-cbind(as.data.frame(up),t(as.data.frame(cof)))
anno<-read.delim("Ncrassa.genes.FungiDB.txt",header=T,sep="\t")
anno<-anno[,1:2]
keep<-list()
for (i in rownames(final)){
  idx<-grep(i,anno$X.Gene.ID.)
  idx<-anno[idx,]
  keep[[i]]<-idx
 } 
df<-melt(keep,id=c("X.Gene.ID.","X.Gene.Name.or.Symbol."))
df<-df[,1:2]
colnames(df)<-c("ens.id","symbol")
df<-df[!duplicated(df), ]
final<-cbind(final,df)

######for single plot 
p<-data.frame (colData(fullData),t(counts(fullData, normalized=TRUE )["NCU00552", , drop=FALSE] ) ) %>%
  mutate ( tl = ifelse( time > 11, time - 11, 1e5 ) )  %>%
  gather ( "gene", "ncount", starts_with("NCU") )

fit <- nls(
  ncount ~ ifelse( time>11, mu + A*exp(-k*(time-11)), muD ),
  algorithm="port",
  start = c( mu=300, A=2^5, k=0.5, muD=300 ), 
  lower = c( mu=1, A=2^0, k=0.01, muD=1 ), 
  upper = c( mu=1e5, A=2^20, k=50, muD=1e5 ), 
  data=subset(p, samplelabel=="dl-wt" & time>0 ), trace=TRUE )

fittedline <- data.frame( time = seq( 0, 22, length.out=1000 ) )
fittedline$tl <- ifelse( fittedline$time > 11, fittedline$time - 11, 1e5 )
fittedline$ncount <- predict( fit, fittedline )
ggplot(data) + 
  geom_point(aes( x=time, y=ncount), data=subset( p, samplelabel=="dl-wt" ) ) +
  geom_line(aes( x=time, y=ncount), data=fittedline ) +
  scale_y_continuous(trans="log2")





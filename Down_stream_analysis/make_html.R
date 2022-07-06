library( hwriter )

########## here the output of html file name is open. images is the folder where all the png files are.
page <- openPage( "open.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:500px; overflow-y:scroll">' )
#####
hwrite(ressig, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", ressig$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="200px"></td></tr></table>' )
closePage(page)
browseURL( "open.html" )

#########another way of writeing html file.
#load( "final.rda" )
final$log2FoldChange <- sprintf( "%.2f", final$log2FoldChange )
plotfiles <- paste0( rownames(final), ".png" )
plotfiles2<-paste0( rownames(final), ".jpeg" )
final$plot <-hwriteImage(plotfiles, height="70px", table=FALSE, link = plotfiles )
final$plot2<-hwriteImage(plotfiles2, height="70px", table=FALSE, link = plotfiles2 )
maintable <- hwrite( final, 
                     table.style = "border-collapse: collapse;",
                     cellpadding = 5,row.bgcolor='#ffdc98',
                     center=TRUE,onmouseover="this.bgColor='#ffaaaa'", 
                     onmouseout="this.bgColor='white'", bgcolor='white')

hwrite(maintable, page="table.html")

#####hetamp2
heatmap.2(as.matrix( dd[,c(1:12)]), col=greenred(30), trace="none",
                    Colv=FALSE, Rowv=FALSE, dendrogram = "none",scale="row",
          labRow = FALSE,cexRow=0.075,main="Significant rythmic genes in dark phase",keysize=1)






########ressig plot 
ressig$log2FoldChange <- sprintf( "%.2f", ressig$log2FoldChange )
plotfiles <- paste0( rownames(ressig), ".png" )
ressig$plot <-hwriteImage(plotfiles, height="70px", table=FALSE, link = plotfiles )
#final$plot2<-hwriteImage(plotfiles2, height="70px", table=FALSE, link = plotfiles2 )
maintable <- hwrite( ressig, 
                     table.style = "border-collapse: collapse;",
                     cellpadding = 5,row.bgcolor='#ffdc98',
                     center=TRUE,onmouseover="this.bgColor='#ffaaaa'", 
                     onmouseout="this.bgColor='white'", bgcolor='white')
hwrite(maintable, page="differential regulated gene table.html")
####write below script 
<script type="text/javascript">
	var row = document.getElementsByTagName('tr');
	var main_url = 'http://fungidb.org/fungidb/app/record/gene/';
	for(i = 2; i < row.length; i++)
	{
		td_last = row[i].lastElementChild;
		gene_name = td_last.textContent
		row[i].lastElementChild.textContent = '';
		gene_url = main_url + gene_name;
		var a_link = document.createElement('a');
		a_link.href = gene_url;
		a_link.text = gene_name;
		row[i].lastElementChild.appendChild(a_link);
	}
</script>
</body>


########extract html file in r 
library(XML)
library(RCurl)
library(xlsx)
extrextract information from webpage 
#http://fungidb.org/fungidb/app/record/gene/
http://nptel.ac.in/courses/109104115/
##click right mouse and look for pageview source, look for conetnt that you look for this 
###use internal node,so what ever html there it will extract 
extracthtml<-htmlTreeParse('http://nptel.ac.in/courses/109104115/',useInternalNodes = TRUE)
###we need to find the extract that is use full, here we have to get the tag, we have to specify the begnining the tag and end tag  
content<-getNodeSet(extracthtml,"//ul//li")
content
x<-length(content)
content<-sapply(content,xmlValue)

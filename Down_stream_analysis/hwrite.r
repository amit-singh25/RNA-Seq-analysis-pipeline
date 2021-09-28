library( hwriter )
ressig$log2FoldChange <- sprintf( "%.2f", ressig$log2FoldChange )
ressig$padj <- sprintf( "%.2e", ressig$padj )
ressig<-ressig[,c(2,6:8)]
page <- openPage( "ligt_induced_gene.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:700px; overflow-y:scroll">' )


#####
hwrite(ressig, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", ressig$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="7re00px"></td></tr></table>' )
closePage(page)
browseURL( "ligt_induced_gene.html" )

###########12hr dark plot 
ressig_rep_3$log2FoldChange <- sprintf( "%.2f", ressig_rep_3$log2FoldChange )
ressig_rep_3$padj <- sprintf( "%.2e", ressig_rep_3$padj )
ressig_rep_3<-ressig_rep_3[,c(2,6:8)]

page <- openPage( "12hr_rep.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
hwrite(ressig_rep_3, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", ressig_rep_3$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="700px"></td></tr></table>' )
closePage(page)
browseURL( "12hr_rep.html" )








sig<-sig[,c(13:18)]


#####################

page <- openPage( "ryth_gene_adj.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
 hwrite(sig, border=NULL, page=page, bgcolor=matrix(c('#aaaaff'),nr=10,nc=6))
       onmouseover = sprintf( "set_image( '%s' );", sig$ens.id)

cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="800px"></td></tr></table>' )
closePage(page)
browseURL( "ryth_gene_adj.html" )









ressig_rep_2$log2FoldChange <- sprintf( "%.2f", ressig_rep_2$log2FoldChange )
ressig_rep_2$padj <- sprintf( "%.2e", ressig_rep_2$padj )
ressig_rep_2<-ressig_rep_2[,c(2,6:8)]
page <- openPage( "12hr_rep.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
hwrite(ressig_rep_2, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", ressig_rep_2$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="800px"></td></tr></table>' )
closePage(page)
browseURL( "10hr_rep.html" )


#######ven list gene 

page <- openPage( "ven_list.html",
                  head = paste( sep="\n",
                                "<script>",
                                "   function set_image( name ) {",
                                "      document.getElementById( 'plot' ).setAttribute( 'src', 'images/' + name + '.png' );",
                                "   }",
                                "</script>" ) )
cat(file=page,
    '<table><tr><td style="vertical-align:top"><div style="height:800px; overflow-y:scroll">' )
#####
hwrite(df, border=NULL, page=page,
       onmouseover = sprintf( "set_image( '%s' );", df$ens.id  ) )
cat( file=page,
     '</div></td><td style="vertical-align:top"><img id="plot" width="800px"></td></tr></table>' )
closePage(page)
browseURL( "ven_list.html" )























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
browseURL( "table.html" )

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

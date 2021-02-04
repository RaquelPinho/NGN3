target <-crispr_set3$target
chrom <- as.character(gd2$chrm)
# This range is for the amplicon
start<- gsub(".*:|-.*", "", as.character(seqlevels(target)))
end<- gsub(".*:|.*-", "", as.character(seqlevels(target)))
data_target <- data.frame(chrom, start, end)
tgt <- makeGRangesFromDataFrame(data_target)

txng<- GenomicFeatures::makeTxDbFromGFF("NGN3_2.gtf", format = "gtf")


# Creating overlap ranges for the plot of the location of target
amplicon<-"ref_seq.bed"
gd_amp <- rtracklayer::import(amplicon)
starts<- start(gd_amp) + start(gd2)
ends <- start(gd_amp) + end(gd2)
ir_guide <- IRanges(starts,ends)

star_a<-rep(1,length(gd_amp))
end_a<-width(gd_amp)
start_t<-start(gd2)
end_t <- end(gd2)
dataplot <- data.frame(star_a,end_a,start_t,end_t )


ggplot2::ggplot(dataplot) + geom_point(aes_(x = quote(seqnames(gd2)), y = quote(ys), group = quote(ys), shape = quote(shp)), size = 2) + geom_line(data = lns, aes_(x = quote(tloc), y = quote(ys), group = quote(ys))) + scale_shape_identity() 
                       
                                                                                                                                                                     
                   color = "black", aes_(xmin = quote(start), xmax = quote(end), 
                                         ymin = quote(ymin), ymax = quote(ymax)))
p <- p + geom_rect(data = target_df, aes_(xmin = quote(xmin), 
                                          xmax = quote(xmax), ymin = quote(ymin), ymax = quote(ymax)), 
                   colour = target.colour, fill = NA, size = target.size)
if (!isFALSE(plot.title)) {
  p <- p + ggtitle(plot.title)
}
p <- p + theme_minimal() + theme(axis.text.x = element_text(size = gene.text.size), 
                                 axis.text.y = element_blank(), panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", 
                                                                                                    colour = NA), panel.background = element_rect(fill = "white", 
                                                                                                                                                  colour = NA), panel.spacing = panel.spacing, text = element_text(size = gene.text.size), 
                                 axis.ticks.y = element_blank()) + ylab(NULL) + xlab(NULL)
return(p)
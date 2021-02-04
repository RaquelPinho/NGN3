t was found in the following places
namespace:CrispRVariants
with value

function (txdb, target, target.colour = "red", target.size = 1, 
          gene.text.size = 10, panel.spacing = grid::unit(c(0.1, 0.1, 
                                                            0.1, 0.1), "lines"), plot.title = NULL, all.transcripts = TRUE) 
{
  genomicfeatures <- requireNamespace("GenomicFeatures")
  stopifnot(isTRUE(genomicfeatures))
  genes <- .getOverlappingGenes(txdb, target, all.transcripts = all.transcripts)
  if ("grob" %in% class(genes)) 
    return(genes)
  result <- .makeGeneSegments(genes, txdb, target)
  all_exs <- result$all_exs
  gene_spans <- result$gene_spans
  min_st <- min(all_exs$start)
  max_end <- max(all_exs$end)
  tcks <- unname(quantile(min_st:max_end, seq(1, 100, by = 2) * 
                            0.01))
  tcks <- lapply(as(gene_spans, "GRangesList"), function(sp) {
    tcks[tcks > start(sp) & tcks < end(sp)]
  })
  tck_lns <- lapply(tcks, length)
  tcks <- data.frame(tloc = unlist(tcks), ys = rep(1:length(tcks), 
                                                   lapply(tcks, length)))
  lns <- data.frame(tloc = c(start(gene_spans), end(gene_spans)), 
                    ys = rep(seq_along(gene_spans), 2))
  all_exs$ymax <- all_exs$ts + 0.3
  all_exs$ymin <- all_exs$ts - 0.3
  is_utr <- all_exs$type == "utr"
  all_exs$ymax[is_utr] <- all_exs$ts[is_utr] + 0.2
  all_exs$ymin[is_utr] <- all_exs$ts[is_utr] - 0.2
  target_df <- data.frame(xmin = start(target), xmax = end(target), 
                          ymin = 0, ymax = ceiling(max(all_exs$ymax)))
  if (is.null(plot.title)) {
    plot.title <- paste(unique(genes$GENEID), sep = ";")
  }
  strands <- rep(as.character(strand(gene_spans)), tck_lns)
  strands[strands == "-"] <- 60
  strands[strands == "+"] <- 62
  tcks$shp <- as.integer(strands)
  p <- ggplot2::ggplot(tcks) + geom_point(aes_(x = quote(tloc), 
                                               y = quote(ys), group = quote(ys), shape = quote(shp)), 
                                          size = 2) + geom_line(data = lns, aes_(x = quote(tloc), 
                                                                                 y = quote(ys), group = quote(ys))) + scale_shape_identity()
  p <- p + geom_rect(data = all_exs, fill = "black", 
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
}
<bytecode: 0x0000020e2aee1940>
  <environment: namespace:CrispRVariants>
  > 
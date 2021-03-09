###################
### functions #####
###################

#' This function creates the annotation plot with genomic features of the regions in a especified GRange file
#' @param txlocigff is a TxDb with the annotation of the genomic regions amp_ranges
#' @param amp_ranges is a GRange containing the Loci information of interest with seqlevels matching txlocigff
#' @param Locus is a especific locus to be plotted
#' @param crisprset is a crisprset object containing the allele variants information of the locus
#' @param gd2 is a GRange file containing the loci in the crisprset and amp_ranges with the names matching crisprset
#' @return annotation plot 
plotAnnotation_Mod <- function(txlocigff, Locus, crisprset, amp_ranges, gd2) {

target.colour <- "red"
guide.colour <- "darkred"
target.size <- 1
gene.text.size <- 10
panel.spacing <- grid::unit(c(0.1, 0.1, 0.1, 0.1), "lines")
plot.title <- NULL

target<- crisprset$target
idx <- which(amp_ranges$ID == Locus)
genes <- CrispRVariants:::.getOverlappingGenes(txlocigff, amp_ranges[idx])
# this results in a data.frame with columns : "GENEID","CDSCHROM", "CDSSTART", "CDSEND", "EXONSTART", "EXONEND", "TXID", "TXNAME","TXSTRAND". 
results <- CrispRVariants:::.makeGeneSegments(genes, txlocigff, amp_ranges[idx]) 
# this results in a list composed of all_exs (a data.frame) and gene_spans (a GRange object)
all_exs <- results$all_exs
gene_spans <- results$gene_spans
min_st <- min(all_exs$start)
max_end <- max(all_exs$end)
tcks <- unname(quantile(min_st:max_end, seq(1, 100, by = 2) *   0.01))
tcks <- lapply(as(gene_spans, "GRangesList"), function(sp) {
  tcks[tcks > start(sp) & tcks < end(sp)]
})
tck_lns <- lapply(tcks, length)
# just the length of every tck vector that will probably always be 50 because we a doing sequence from 1-100 by 2 in the steps above.
tcks <- data.frame(tloc = unlist(tcks), ys = rep(1:length(tcks),  lapply(tcks, length)))
# Now tck will be the data.frame used to plot
lns <- data.frame(tloc = c(start(gene_spans), end(gene_spans)), 
                  ys = rep(seq_along(gene_spans), 2))
# seq_along creats a vector with the seq from 1 to the number of ranges in gene spans.
all_exs$ymax <- all_exs$ts + 0.3
all_exs$ymin <- all_exs$ts - 0.3
is_utr <- all_exs$type == "utr"
all_exs$ymax[is_utr] <- all_exs$ts[is_utr] + 0.2
all_exs$ymin[is_utr] <- all_exs$ts[is_utr] - 0.2
all_exs$colour <-ifelse(all_exs$type == "exon", "black","white")
target_df <- data.frame(xmin = start(amp_ranges[idx]), xmax = end(amp_ranges[idx]), 
                        ymin = 0, ymax = ceiling(max(all_exs$ymax)))
# ceiling just round up the number for the smalest integer possible.
# This I need to change because it has gd2
if (is.null(plot.title)) {
  plot.title <- paste(unique(amp_ranges$ID[idx], gd2$genes_names[gd2$names == amp_ranges$ID[idx]]))
}
strands <- rep(as.character(strand(gene_spans)), tck_lns)
strands[strands == "-"] <- 60
strands[strands == "+"] <- 62
tcks$shp <- as.integer(strands)
# this will give the direction of the shape in the plot depending on the strand direction.
p <- ggplot2::ggplot(tcks) + geom_point(aes_(x = quote(tloc), 
                                             y = quote(ys), group = quote(ys), shape = quote(shp)), size = 2) + geom_line(data = lns, aes_(x = quote(tloc),  y = quote(ys), group = quote(ys))) + scale_shape_identity()
# this creates  a line from the start of the gene to the end of the gene, I need to put scale_shape_identity() so it doesn't give me an error. it says that it says that the data has already been scaled.

p <- p + geom_rect(data = all_exs, aes_(xmin = quote(start), xmax = quote(end), ymin = quote(ymin), ymax = quote(ymax)), colour = "black", fill = all_exs$colour)
# This put black rectangles in the features of the gene (utr, exons)  

#p <- p + geom_rect(data = target_df, aes_(xmin = quote(xmin),  xmax = quote(xmax), ymin = quote(ymin), ymax = quote(ymax)),    colour = target.colour, fill = NA, size = target.size)
p <- p + geom_rect(data = target_df[1,], aes_(xmin = quote(xmin),  xmax = quote(xmax), ymin = quote(ymin), ymax = quote(ymax)),    colour = target.colour, fill = NA, size =1)
# This will add the amplicon region
p <- p + geom_rect(data = target_df[1,], aes_(xmin = quote(xmin + start(target)),  xmax = quote(xmin + end(target)), ymin = quote(ymin), ymax = quote(ymax)),    colour = "darkred", fill = "darkred", size =1)

# This will add the region that is being visualized

if (!isFALSE(plot.title)) {
  p <- p + ggtitle(plot.title)
}

p <- p + theme_minimal() + theme(axis.text.x = element_text(size = gene.text.size), axis.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white",  colour = NA), panel.background = element_rect(fill = "white",   colour = NA), panel.spacing = panel.spacing, text = element_text(size = gene.text.size), axis.ticks.y = element_blank()) + ylab(NULL) + xlab(NULL)


return(p)

}


#' This function combine the annotation, alignment and freqHeatmap plots in one plot
#' @param top.plot is the annotation plot
#' @param right.plot is the freqHeatmap 
#' @param left.plot is the alignement plot
#' @return merged plot 
arrangePlots <- function(top.plot, left.plot, right.plot, fig.height = NULL,
                         col.wdth.ratio  = c(2, 1), row.ht.ratio = c(1,6),
                         left.plot.margin = grid::unit(c(0.1,0,3,0.2), "lines")){
  
  # Set the size ratio of the top and bottom rows
  plot_hts <- if (is.null(fig.height)){ row.ht.ratio
  }else { fig.height/sum(row.ht.ratio)*row.ht.ratio }
  
  # Remove y-axis labels from right plot
  right.plot <- right.plot + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
  
  # Adjust margins of left.plot
  left.plot <- left.plot + theme(plot.margin = left.plot.margin)
  
  # Convert plots to grobs, lock plot heights
  p2 <- ggplot2::ggplotGrob(left.plot)
  p3 <- ggplot2::ggplotGrob(right.plot)
  p3$heights <- p2$heights
  
  # Return arranged plots
  return(gridExtra::grid.arrange(top.plot,
                                 gridExtra::arrangeGrob(p2, p3, ncol = 2, widths = col.wdth.ratio),
                                 nrow = 2, heights = plot_hts, newpage = FALSE))
}



#' This function returns the mutation eficiencies and the variant frequency matrix.
#' @param crisprset is a crisprset object
#' @return mut_eff is a matrix with the mutation efficency calculated with / without/ minus the read couts / 
#' frequency from the control samples (what is removed in the minus is the greater frequency when more than one control sample).
#' @return vcc is the mutation frequency matrix
mut.efficiency <- function(crisprset){
	# get the variant counts 
	vc_all <-CrispRVariants::variantCounts(crisprset)

	# calculate variant frequencies 
	total_reads <-apply(vc_all,2,sum)
	vcc <- sweep(vc_all, 2,total_reads , '/')

	# Calculation of mutation efficiency including control reads
	rows_SNV <- grep("SNV", rownames(vc_all))
	mut_reads <- apply(vc_all[-rows_SNV,],2,sum)
	mut_eff <- mut_reads / total_reads

	# Calculation of mutation efficiency excluding control reads
	controls <- grep("CL", colnames(vc_all))
	reads_control <- which(vc_all[,controls[1]] != 0 | vc_all[,controls[2]] != 0)
	rm_reads <- unique(c(reads_control, rows_SNV))
	mut_reads_nc <- apply(vc_all[-rm_reads,], 2, sum)
	mut_eff_nc <- mut_reads_nc / total_reads

	# Calculating mutation efficiency normilized by the control reads 
	# Get higher frequencies between controls 
	vcc_mfmc <- vcc
	maxfreq<- apply(vcc[-1, controls],1, max)
	mut_freq_minus_control <- as.vector(sweep(vcc[-1,],1, maxfreq,"-"))
	dims<- dim(mut_freq_minus_control)
	mut_freq_minus_control <- ifelse(mut_freq_minus_control < 0, vcc[-1,], mut_freq_minus_control)
	dim(mut_freq_minus_control)<- dims
	vcc_mfmc[-1,] <- mut_freq_minus_control
	mut_eff_mc <- apply(vcc_mfmc[-rows_SNV,],2,sum)

	#Creating the summary (mean , min and max efficiency)
	mut_eff_dt <- as.data.frame(rbind(mut_eff, mut_eff_nc))
	mut_eff_dt <- as.data.frame(rbind(mut_eff_dt, mut_eff_mc))
	sum_m <- mut_eff_dt %>% t() %>% as.data.frame() %>% summarise_all(.funs = list(
	         mean = mean,
	         min = min,
	         max = max), na.rm=TRUE)  %>% matrix(3,3)
	colnames(sum_m) <- c( "Mean", "Min", "Max")
	rownames(sum_m) <- rownames(mut_eff_dt)
	mut_eff <- transform(merge(mut_eff_dt, sum_m, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
	rownames(mut_eff) <- c("Mut.efficiency.minus.control", "Mut.efficiency", "Mut.efficiency.wo.controlreads")
	
		return(list(mut_eff = mut_eff, vcc = vcc ))
}


#' Get the rownames of readMatrix above a given threshold
#' @param col a column from readMatrix
#' @param dataConservationThreshold a threshold in terms of lost data
#' @return a list of rownames from readMatrix
filter.alleles <- function(col, .dataConservationThreshold = 0.90){
  col <- sort(col)
  csum.col <- cumsum(col)
  alleles <- names(col)[which(csum.col >= (1 - .dataConservationThreshold) * sum(col))]
  data_alleles_kept <-col[alleles]
  min <- min(data_alleles_kept)
  # Checking if there is repetitive minimal value (this indicates if there were exclusion and retention of alleles with same number of reads - To remove this bias on having to choose which allele to remove I will put the minimum as the next lower value )
  if(length(data_alleles_kept[data_alleles_kept == min]) > 1){
    data_alleles_kept <- data_alleles_kept[data_alleles_kept > min]
    min <- min(data_alleles_kept)
    alleles <- names(data_alleles_kept)
    }
    return(alleles)
}


#' Get the rownames of readMatrix above a given threshold
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return a list of rownames from readMatrix from all samples
filter.allAlleles <- function(readMatrix, dataConservationThreshold = 0.90){
  filtered.rownames <- apply(readMatrix, 2, filter.alleles, .dataConservationThreshold = dataConservationThreshold)
  unique(unlist(filtered.rownames))
}


#' Get the rownames of readMatrix above a given threshold
#' @param readMatrix a matrix with number of reads per allele per sample
#' @return plots of conservation in number of reads or percentage
plotConservation <- function(readMatrix, thresholds = seq(0.80, 0.99, by = 0.01)){
  #percentage of unique alleles conserved (kept) when using the threshold
  actualDataConservation <- sapply(thresholds, function(t) length(filter.allAlleles(readMatrix, dataConservationThreshold = t)) / nrow(readMatrix))
  #number of alleles kept when using the threshold
  allelesConserved <- sapply(thresholds, function(t) length(filter.allAlleles(readMatrix, dataConservationThreshold = t)))
  df <- data.frame(thresholds = thresholds, dataConservation = actualDataConservation, alleleConservation = allelesConserved)
  df_m <-reshape2::melt(df, id.vars = "thresholds")
  ggplot(df_m, aes(x = thresholds, y = value, group = variable)) +
   geom_point() +
    facet_grid(variable ~ ., scales = "free", labeller = label_parsed, switch = "y") + geom_vline(xintercept = 0.9, linetype = "dashed", color = "red" ) +
   ggpubr::theme_pubr() + labs(title = "Allele conservation per read-threshold", x = "Threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5))
 
 }

#' Get the rownames of readMatrix above a given threshold
#' @param col a colunm from readMatrix with number of reads per allele in a given sample
#' @return plots of conservation in number of reads or percentage
plotConservationSample <- function(col, thresholds = seq(0.80, 0.99, by = 0.01)){
  #percentage of unique alleles conserved (kept) when using the threshold
  col<- col[col > 0]
  actualDataConservation <- sapply(thresholds, function(t) length(filter.alleles(col, .dataConservationThreshold = t)) / length(col))
  #number of alleles kept when using the threshold
  allelesConserved <- sapply(thresholds, function(t) length(filter.alleles(col, .dataConservationThreshold = t)))
  df <- data.frame(thresholds = thresholds, dataConservation = actualDataConservation, alleleConservation = allelesConserved)
  df_m <-reshape2::melt(df, id.vars = "thresholds")
  ggplot(df_m, aes(x = thresholds, y = value, group = variable)) +
   geom_point() +
    facet_grid(variable ~ ., scales = "free", labeller = labeller("dataConservation", "alleleConservation"), switch = "y") + geom_vline(xintercept = 0.9, linetype = "dashed", color = "red" ) +
   ggpubr::theme_pubr() + labs(title = "Allele conservation per read-threshold", x = "Threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5), panel.spacing = unit(2,"lines"), axis.title.x = element_text(margin = margin(t = 100, r = 0, b = 0, l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 100, b = 0, l = 0)))
 
 }

#' Get the rownames of readMatrix above a given threshold
#' @param col the name column from readMatrix or the column of readMatrix
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return plots of cumulative density with trheshold of removed alleles 
#' @return table of mutation efficiency
Mut.eff.threshold <- function(col, .dataConservationThreshold = 0.9){
  ref_name <- names(col[1])
  names.kept <- filter.alleles(col, .dataConservationThreshold) 
  col_k <- col[names.kept]
  mut.eff <- sum(col_k[which(names(col_k)!= ref_name)])/sum(col_k)
  return(mut.eff)
}
#' Get the rownames of readMatrix above a given threshold
#' @param col column from readMatrix or the column of readMatrix
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return plots of cumulative density with trheshold of removed alleles 
#' @return table of mutation efficiency
Plot.Cum.dens <- function(col, .dataConservationThreshold = 0.9) {
  names.kept <- filter.alleles(col, .dataConservationThreshold) 
  col_k <- col[names.kept]
  min <- min(col_k)
  col <- as.data.frame(col)
  colnames(col)<- "Counts"
  p <- ggplot(col, aes(Counts)) + stat_ecdf(geom =  "point") + ylim(min = 0.9, max = 1) + xlim(min = 0, max = (max(col)+500)) + geom_vline(xintercept = min, color = "red", linetype = "dashed" )+
  ggpubr::theme_pubr() + xlab("Read counts") + ylab("Cumulative density") + labs(title = "Allele conservation")
 
  return(ggplotly(p))
}


#' Get the mutation efficiency of a specific sample and general stats of readMatrix when using a given threshold and removing control samples
#' @param colname the name column from readMatrix or the column of readMatrix
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return table of mutation efficiency
#removing controls from calculation
Mut.eff.all.threshold <- function(colname, readMatrix, .dataConservationThreshold = 0.9) {
	mut.eff.list <- apply(readMatrix,2, Mut.eff.threshold)
	cl<-grep("CL",colnames(readMatrix))
    min<- min(mut.eff.list[-cl])
    max <- max(mut.eff.list[-cl])
    mean <- mean(mut.eff.list[-cl])
    sd <- sd(mut.eff.list[-cl])
    col <- readMatrix[,colname]
    mut.eff.col <- Mut.eff.threshold(col)
    df <- data.frame(mu.eff.sample = mut.eff.col, min.geral = min, mean.geral = mean, max.geral = max, sd.geral = sd)
    return(df)
}

#' Get the heatmap of samples missing in each locus 
#' @param listMatrix the list of variant counts for each locus
#' @return heatmap of samples missing in data
MatrixSamples <- function(listMatrix) {
  #color
  palette <- wesanderson::wes_palette("GrandBudapest1",3)
  # get the samples in each locus
  list_samples <- lapply(listMatrix, colnames)
  # get total number of samples
  total_number_samples <- max(unlist(lapply(list_samples,length)))
  # get complete list of samples
  total_name_samples <- list_samples[[which(unlist(lapply(list_samples, length)) == total_number_samples)]]
  # create matrix for heatmap
  binary_list <- lapply(list_samples, FUN = function(x) { x <- total_name_samples %in% x
    x <-ifelse(x == TRUE,1,0)
    })
  
  binary_list <- as.matrix(bind_rows(binary_list))
  
  rownames(binary_list) <- total_name_samples
  
  p <- superheat(t(binary_list), order.rows = order(1:ncol(binary_list), decreasing = T),  heat.pal = palette, heat.lim = c(0,1), heat.na.col = "#b3e2cd",left.label.col = "lavenderblush1", left.label.size = 0.1, left.label.text.size = 4, bottom.label.size = 0.2, bottom.label.text.size = 4, bottom.label.text.angle = 90, bottom.label.col = "lavenderblush1", title = "Samples missing per locus" , row.title = "Locus", column.title = "Samples", legend = F )
        return(p)
   
  }

#' Get the heatmap of samples missing in each locus 
#' @param listMatrix the list of variant counts for each locus
#' @return heatmap of samples missing in data
MissingSamples <- function(listMatrix) {
  # get the samples in each locus
  list_samples <- lapply(listMatrix, colnames)
  # get total number of samples
  total_number_samples <- max(unlist(lapply(list_samples,length)))
  # get complete list of samples
  total_name_samples <- list_samples[[which(unlist(lapply(list_samples, length)) == total_number_samples)]]
  # create matrix for heatmap
  missing_list <- lapply(list_samples, FUN = function(x) { x <- total_name_samples %in% x
          x <- total_name_samples[!x]
          x <- paste(x, collapse = ", ")
              })
  missing_df <- data.frame(Locus = names(missing_list), Missing_samples = unlist(missing_list, use.names = FALSE))
  return(missing_df)
}

#' add information of missing samples into read counts -which means adding 0 to missing samples entry in total read counts per sample list
#' @param listMatrix the list of variant counts (read counts per allele per sample) for each locus
#' @return missing samples entry in total read counts per sample list including missing samples
SumReadsPerSample <- function(listMatrix){
# Getting the total read counts per sample per locus
total_read_counts_per_samples_per_locus<- lapply(listMatrix, FUN = function(x) {apply(x,2, sum)})
# Adding missing samples data
n_samples <- sapply(total_read_counts_per_samples_per_locus, length)
n_samples_max <- max(n_samples) 
sample_names <- sapply(total_read_counts_per_samples_per_locus[which(n_samples == n_samples_max)[1]],names)
total_read_counts_complete <- lapply(total_read_counts_per_samples_per_locus, FUN = function(x) {
    mis_samples <- !(sample_names %in% names(x))
    names(mis_samples) <- sample_names
    c_sum  <- c()
    for(i in 1:length(sample_names)) {
    c_sum[i] <- ifelse(mis_samples[i] == FALSE , x[names(mis_samples[i])],0)
      }
    names(c_sum)<- sample_names
    return(c_sum)
  })
return(total_read_counts_complete)
}

#' Get the data table of basic statistics for read counts 
#' @param listMatrix the list of variant counts (read counts per allele per sample) for each locus
#' @return basic stats without removing missing samples - that means considering those 0
GetStatsReads <- function(listMatrix){
# Getting the total read counts per sample per locus
total_read_counts_complete <- SumReadsPerSample(listMatrix)
# removing samples with 0 reads 
  total_read_counts_complete <- lapply(total_read_counts_complete, FUN = function(le) {
    no_reads <- which(le == 0)
    if(!(purrr::is_empty(no_reads))) {
    le <- le[- no_reads]
    }
     return(le)
  } )
# Getting total reads per locus
total_reads_per_locus <- sapply(total_read_counts_complete,sum)
# Getting total reads 
total_reads_general <- sum(unlist(total_reads_per_locus))
# Getting stats summary information
# per locus
sum_per_locus <- list(psych::describe(total_reads_per_locus))

# per sample
sum_per_sample <- lapply(total_read_counts_complete,psych::describe)

# Combining information 
sum_info <- append(sum_per_sample, sum_per_locus)
names(sum_info) <- c(names(sum_per_sample), "Total")

# Creating a data frame with information of interest
options(scipen = 999)
dt <- as.matrix(dplyr::bind_rows(sum_info))
rownames(dt) <- names(sum_info)
dt <- dt[,c(2,3,4,5,8,9,13)]
# returning info data
return(dt)
}

#' Get barplot of the mean read counts per group 
#' @param listMatrix the list of variant counts (read counts per allele per sample) for each locus
#' @param dataGroup vector with metadata of samples. Should be in the same order of the colnames in listMatrix
#' @return barplot of mean of read counts per group in a ggplot object
PlotReadPerGroup <- function(listMatrix,groupData) {
     #Get the total read counts per Loci per sample
     total_read_counts_complete <- SumReadsPerSample(listMatrix)
     #Contruct data.frame for plot /containing total reads per sample and the sample group information
     group_info <- lapply(seq_along(total_read_counts_complete), FUN = function(i) {
       x <- data.frame(
         Counts = as.numeric(total_read_counts_complete[[i]]), 
         Group = groupData,
         Loci = names(total_read_counts_complete)[i])
     })
     group_info <- do.call(rbind, group_info)
     dt <- group_info %>% 
       group_by(Loci) %>%
       group_modify(~ group_by(., Group) %>%
         summarize(mean = mean(Counts), sd = sd(Counts))
       )

    #Create plot 
     palette <- wesanderson::wes_palette("GrandBudapest1",4)
     p <- ggplot(dt) +  
     aes(Group, mean, fill = Group) + 
     scale_fill_manual(values= palette)+
     geom_col(position = "dodge") + 
     facet_wrap(~ Loci, scales = "free_y") + 
     labs( x = "Group", y= "Mean", title = "Mean of read counts per group") +
     ggpubr::theme_pubr() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     theme(legend.position = c(0.85, 0.1), legend.text = element_text(size =12))+
     scale_x_discrete(labels=c("Blastocysts" = "Blasts", "Blastocysts_M" = "Blasts_M", "Control" = "Control", "Partenotes_M" = "Parte_M"))
     return(p)
}




#' Get the data information of the alleles kept for each sample a given threshold
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return a list of numeric data from the alleles kept for each sample
filter.Alleles.allSamples <- function(readMatrix, dataConservationThreshold = 0.90){
  filtered.rownames <- apply(readMatrix, 2, filter.alleles, .dataConservationThreshold = dataConservationThreshold)
  filtered.data <- lapply(seq_along(filtered.rownames), FUN = function(i){
    filtered.data <- readMatrix[filtered.rownames[[i]], names(filtered.rownames)[i]]
    names(filtered.data) <- filtered.rownames[[i]]
    filtered.data
    })
  names(filtered.data) <- names(filtered.rownames)
  return(filtered.data)
}

#' Get a melted data frame of the alleles per sample to plot allele conservation general 
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param dataConservationThreshold a threshold in terms of conserved data, that is the proportion of reads to be kept.
#' @return a list of numeric data from the alleles kept for each sample
meltedAlleleInfo <- function(readMatrix, dataConservationThreshold = 0.90){
  filtered.data <- filter.Alleles.allSamples(readMatrix,dataConservationThreshold)
  melted.data <- lapply(seq_along(filtered.data), FUN = function(i) {
    x <- data.frame(
         Counts = as.numeric(filtered.data[[i]]), 
         Sample = names(filtered.data)[i])
  })
  melted.data <- do.call(rbind,melted.data)
  melted.data$Alleles <- unlist(unname(sapply(filtered.data,names)))
  return(melted.data)
}

#' PLot the minimal read counts from samples from readMatrix above a given threshold
#' @param readMatrix a matrix with number of reads per allele per sample
#' @return plots of minimal number of reads per thrshold per sample
plotMin_Allsamples <- function(readMatrix, thresholds = seq(0.80, 0.99, by = 0.01)){
  # list of filtered data by threshold
  filtered.data.list <- lapply(thresholds, function(t){
    filtered.data <- meltedAlleleInfo(readMatrix, t)
      })
  names(filtered.data.list) <- thresholds
  minCount_perSample <- lapply(thresholds, function(t){
         x<-  as.data.frame(filtered.data.list[[as.character(t)]] %>% group_by(., Sample) %>% summarise(min = min(Counts)))
         x$Threshold <- t 
         x
  })
  minCount_perSample <- do.call(rbind, minCount_perSample)
  
  p <- ggplot() +geom_point(data = minCount_perSample, aes(x = Threshold, y = min, fill = Sample, color = Sample )) + ggpubr::theme_pubr()
  return(p)
 
 }

#' Adjust the distance of labels of x-axis and y-axis of a ggplotly obj
#' @param gg a ggplotly object
#' @param x the adjstment of the distance of the x-axis label to the x-axis (more negative, more distant)
#' @param y the adjstment of the distance of the y-axis label to the y-axis (more negative, more distant)
#' @return plot with adjusted margins for the labels
layout_ggplotly <- function(gg, x = -0.06, y = -0.05){
  # The 1 and 2 goes into the list that contains the options for the x and y axis labels respectively
  #adjusting margins of x axis label
   gg$x$layout$annotations[[1]]$y <- -0.06
  #adjusting margins of y axis label
   gg$x$layout$annotations[[2]]$x <- -0.05
   gg
}

#' Get the vector with type of mutation of each allele.
#' @param allele_name a character vector with the allele name (information)
#' @return a character vector with the mutation type
getAlleleType <- function(allele_name){
        if(grepl('SNV', allele_name)){
          return("SNV")
        } else if(grepl("(?=.*D)(?=.*I)",allele_name, perl = T)){
          return("Mixed")
        } else if(grepl('I', allele_name)){
          return("Insertion")
        } else if(grepl('D', allele_name)){
          return("Deletion")
        } else if(grepl('no', allele_name)){
          return("Reference") 
        } else{
          return("Other")
        } 
      }

#' Create a matrix for circlize data, plus metadata for samples and metadata for variants
#' @param readMatrix of the locus of interest containing read counts per allele per sample
#' @param dataConservationThreshold threshold for allele conservation (is called by the function meltedAlleleInfo)
#' @param groupInfo vector containing the group information on the same order as the col names in readMatrix
#'@return a list containing a matrix, a metadata for the samples and a metadata for the variants
getDTcirclize <- function(readMatrix,groupInfo, dataConservationThreshold = 0.9) {
  #Get allele information after filtering in a data frame of sample per allele
  filtered.data <- meltedAlleleInfo(readMatrix, dataConservationThreshold)
  filtered.data <- reshape2::dcast(filtered.data, Sample ~ Alleles, value.var = 'Counts')
  # Replace NA by 0
  filtered.data[is.na(filtered.data)] <-0
  # Replace read counts by frequency
  total_reads <- apply(filtered.data[,-1],1,sum)
  filtered.data[,-1] <- sweep(filtered.data[,-1], 1,total_reads , '/')
  # Put samples as row names
  rownames(filtered.data) <- filtered.data %>% dplyr::pull(Sample)
  filtered.data <- filtered.data[,-1]
  filtered.data <- as.matrix(filtered.data)
  # Making the sample metadata 
  group_data <- groupInfo
  pal_complete <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush1","peachpuff1", "rosybrown1", "#B3e2cd", "darkslategray3","lightseagreen", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  Meta_samples <- data.frame(Samples = rownames(filtered.data),Total = total_reads, Group = groupInfo[metadata$Samples %in% rownames(filtered.data)]) 
  color<- sapply(Meta_samples$Group,FUN= function(x) {
                        ifelse(x == 'Blastocysts',"#E6A0C4", 
                        ifelse(x == 'Blastocysts_M',"#C6CDF7",
                        ifelse(x == 'Control', "#D8A499","#7294D4")))
                      })
  Meta_samples$Color <-color
  # Meta for the variants 
  pal_variantes <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush3","#b3e2cd","lightseagreen")
  type <- sapply(colnames(filtered.data), getAlleleType)
  color_v <- sapply(type, FUN = function(type) {switch(type, 
    'SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")}) 

  Meta_variants <- data.frame(Variants = colnames(filtered.data), Type = type, Color = color_v)
  # return data in a list
  return(list(matrix = filtered.data, MetadataSample = Meta_samples, MetadataVariants = Meta_variants))
}

#' Get a melted data frame of the alleles per sample to plot allele conservation general 
#' @param readMatrix a matrix with number of reads per allele per sample
#' @param freq.threshold a frequency threshold for filtereing the alleles.
#' @return a list of numeric data from the alleles kept for each sample
meltedAlleleInfoFreqFiltered <- function(readMatrix, freq.threshold = 0.0625){
  filtered.data <- frequencyFilterPerSample(readMatrix,freq.threshold)
  melted.data <- lapply(seq_along(filtered.data), FUN = function(i) {
    x <- data.frame(
         Frequency = as.numeric(filtered.data[[i]]), 
         Sample = names(filtered.data)[i])
  })
  melted.data <- do.call(rbind,melted.data)
  melted.data$Alleles <- unlist(unname(sapply(filtered.data,names)))
  return(melted.data)
}

#' Create a matrix for circlize data, plus metadata for samples and metadata for variants
#' @param readMatrix of the locus of interest containing read counts per allele per sample
#' @param freq.threshold threshold for allele conservation (is called by the function meltedAlleleInfo)
#' @param group character string designating which group to use : "Blastocysts", "Blastocysts_M", "Parternotes_M"
#' @return a list containing a matrix, a metadata for the blastocyst samples and a metadata for the variants
getDTcirclizeByGroupFreqFiltering <- function(readMatrix, freq.threshold = 0.0625, group = "Blastocysts" ) {
  #Get allele information after filtering in a data frame of sample per allele
  filtered.data <- meltedAlleleInfoFreqFiltered(readMatrix, freq.threshold =0.0625)
  if(group == "Blastocysts" ){
    filtered.data <- filtered.data[-(grep("M", filtered.data$Sample)),]
  }
  if(group == "Blastocysts_M" ){
    filtered.data <- filtered.data[c(grep("BLM", filtered.data$Sample),grep("BLCL", filtered.data$Sample)),]
  }
  if(group == "Partenotes_M" ){
    filtered.data <- filtered.data[c(grep("PTM", filtered.data$Sample),grep("PTMCL", filtered.data$Sample)),]
  }
  filtered.data <- reshape2::dcast(filtered.data, Sample ~ Alleles, value.var = 'Frequency')
  # Replace NA by 0
  filtered.data[is.na(filtered.data)] <-0
  # Put samples as row names
  rownames(filtered.data) <- filtered.data %>% dplyr::pull(Sample)
  filtered.data_m <- filtered.data[,-1]
  if(is.vector(filtered.data_m)){
    sprintf(" WARNING: There is only one allele (%s) present in the samples",colnames(filtered.data)[2])
  }
  if(!is.vector(filtered.data_m)) {
    filtered.data <- as.matrix(filtered.data_m)
  }
  # Making the sample metadata 
  group_data <- metadata$Group[match(rownames(filtered.data),metadata$Samples)]
  pal_complete <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush1","peachpuff1", "rosybrown1", "#B3e2cd", "darkslategray3","lightseagreen", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  alleles_list <- frequencyFilterPerSample(readMatrix[,match(rownames(filtered.data),colnames(readMatrix))])
  if(is.list(alleles_list)){
    total_reads <- lapply(seq_along(alleles_list), FUN = function(i) {
      total_reads <- readMatrix[names(alleles_list[[i]]),names(alleles_list)[i]]
      total_reads <- sum(total_reads)
    })
  }
  if(is.vector(alleles_list)){
    total_reads <- lapply(seq_along(alleles_list),FUN = function(i) {
      total_reads <- readMatrix[colnames(filtered.data)[2],names(alleles_list[i])]
        })
  }
  total_reads <- unlist(total_reads)
  Meta_samples <- data.frame(Samples = rownames(filtered.data),Total = total_reads, Group = group_data) 
  color<- sapply(Meta_samples$Group,FUN= function(x) {
    ifelse(x == 'Blastocysts',"#E6A0C4", 
           ifelse(x == 'Blastocysts_M',"#C6CDF7",
                  ifelse(x == 'Control', "#D8A499","#7294D4")))
  })
  Meta_samples$Color <-color
  # Meta for the variants 
  pal_variantes <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush3","#b3e2cd","lightseagreen")
  type <- sapply(colnames(filtered.data), getAlleleType)
  color_v <- sapply(type, FUN = function(type) {switch(type, 
                                                       'SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")}) 
  
  Meta_variants <- data.frame(Variants = colnames(filtered.data), Type = type, Color = color_v)
  #reducing the name of the SNV variants
  snv_v <- grep("SNV",colnames(filtered.data))
  colnames(filtered.data)[snv_v] <- paste0("SNV",seq(1:length(snv_v)))
  Meta_variants$VariantsName <- colnames(filtered.data)
  
  # return data in a list
  return(list(matrix = filtered.data, MetadataSample = Meta_samples, MetadataVariants = Meta_variants))
}
  



#' Creating chord diagram for the variant frequency per sample
#' @param readMatrix of the locus of interest containing read counts per allele per sample
#' @param dataConservationThreshold threshold for allele conservation (is called by the function meltedAlleleInfo)
#' @param groupInfo vector containing the group information on the same order as the col names in readMatrix
#'@return a chord diagram plot
VizCirclize <- function(readMatrix,groupInfo, dataConservationThreshold = 0.9 , show_no_variants = FALSE,  freq.threshold = 0.001, scaled = FALSE){
 
  # formatting the data
  circ_data <- getDTcirclize(readMatrix, groupInfo, dataConservationThreshold)
  ## getting the data 
  m <- circ_data[[1]]
  meta_samples <- circ_data[[2]]
  meta_variants <- circ_data[[3]]
  # organizing and filtering data
  ## ordering mutations
  order_m <- order(meta_variants$Type)
  m <- m[,order_m]
  ## set the no variant presence
  if(show_no_variants == FALSE) {
    no_v <- which(colnames(m) == "no variant")
    m[,no_v] <- 0.001
  } 
  
  # color set up
  grid.col <- meta_samples$Color[meta_samples$Samples %in% rownames(m)]
  grid.col <- c(grid.col, meta_variants$Color[meta_variants$Variants %in% colnames(m)])
  col_grad <-circlize::colorRamp2(range(m), c("darkslategray3","#FD6467")) 
  col_grad_m <- col_grad(m)
  col_grad_m[m <= freq.threshold] <- "#00000000" 
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_mutations <- c('SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")
  
  # initiate plot
  circlize::circos.clear()
  par(mar = rep(0, 4), cex=.75)
  circlize::circos.par(start.degree = -90)
  
  circlize::chordDiagram(x = m, directional = 1, 
               transparency = 0.3,
               grid.col = grid.col,
               link.sort = TRUE,
               link.decreasing = TRUE,
               symmetric = FALSE,
               diffHeight = 0,
               link.lwd = 1,
               col = col_grad_m,
               scale = TRUE,
               annotationTrack = NULL,
               preAllocateTracks = list(
                 list(track.height = 0.03, track.margin = c(0, 0)),
                 list(track.height = 2.25 * circos.par("track.height")),
                 list(track.height = 0.02, track.margin = c(0, 0)),
                 list(track.height = 0.02, track.margin = c(0, 0)),
                 list(track.height = 0.02, track.margin = c(0, 0))
               ),
               big.gap = 30)
  
  for (sample in row.names(m)){
    circlize::highlight.sector(sample,
                     track.index = 2,
                     text = paste(sample, meta_samples$Total[ meta_samples$Samples == sample], sep = "    "),
                     cex = 1,
                     col = NA, 
                     border = NA, 
                     facing = "reverse.clockwise",
                     niceFacing = FALSE
    )
  }
  
  for (allele in colnames(m)){
    circlize::highlight.sector(allele,
                     track.index = 2,
                     text = allele,
                     cex = 1,
                     col = NA, 
                     border = NA, 
                     facing = "clockwise",
                     niceFacing = TRUE
    )
  }
  
  for (i in seq_along(levels(meta_samples$Group))){
    Sample_group <- levels(meta_samples$Group)[i]
    ind <- which(meta_samples$Group == Sample_group)
    circlize::highlight.sector(rownames(m)[ind],
                     track.index = 2,
                     text = Sample_group, 
                     col = NA, 
                     border = pal_groups[i], 
                     facing = "bending.inside",
                     niceFacing = FALSE,
                     text.vjust = "40mm",
                     cex = 1
    )
    circlize::highlight.sector(rownames(m)[ind], track.index = c(3,4,5), col = pal_groups[i])
  }
  
  for (i in seq_along(levels(meta_variants$Type))){
    Mut_type <- levels(meta_variants$Type)[i]
    ind <- which(meta_variants$Type[order(meta_variants$Type)] == Mut_type)
    circlize::highlight.sector(colnames(m)[ind],
                     track.index = 2,
                     text = ifelse(Mut_type %in% c("Other", "Reference"), "",Mut_type),
                     col = NA, 
                     border = pal_mutations[Mut_type], 
                     facing = "bending.inside",
                     niceFacing = FALSE,
                     text.vjust = "40mm",
                     cex = 1
    )
    circlize::highlight.sector(colnames(m)[ind], track.index = c(3,4,5), col = pal_mutations[Mut_type])
  }
}


## Calculating allele frequencies
## Convertes a list of read count matrix into a list of frequency matrix
#'@param readMatrix is a matrix of read counts by alleles per sample
#'@return a freqMatrix 
Convert2Freq <- function(readMatrix) {
      total_reads <-apply(readMatrix,2,sum)
      vcc <- sweep(readMatrix, 2,total_reads , '/')
    }
## Filtering out all the alleles with frequency lower than freq.threshold considering all samples
frequencyFilterAllSamples <- function(readMatrix, freq.threshold = 0.0625) {
  freqMatrix <- Convert2Freq(readMatrix)
  allele_list <- lapply(seq_along(colnames(freqMatrix)), FUN = function(i) {
    x <-freqMatrix[,i]
    x <- x[x> freq.threshold]})
  allele_names_all <- unique(unlist(lapply(allele_list,names)))
  return(allele_names_all)
}

## list of alleles with frequency lower than the freq.threshold per sample
frequencyFilterPerSample <- function(readMatrix, freq.threshold = 0.0625) {
   freqMatrix <- Convert2Freq(readMatrix)
   allele_list <- apply(freqMatrix,2, FUN = function(x) {x[x> freq.threshold]})
   return(allele_list)
 }


#' Create a matrix for circlize data, plus metadata for samples and metadata for variants
#' @param readMatrix of the locus of interest containing read counts per allele per sample
#' @param freq.threshold threshold for allele frequency (is called by the function meltedAlleleInfo.2
#' @param groupInfo vector containing the group information on the same order as the col names in readMatrix
#'@return a list containing a matrix, a metadata for the samples and a metadata for the variants
getDTcirclize.2 <- function(readMatrix,groupInfo, dataConservationThreshold = 0.9) {
  #Get allele information after filtering in a data frame of sample per allele
  filtered.data <- meltedAlleleInfo(readMatrix, dataConservationThreshold)
  filtered.data <- reshape2::dcast(filtered.data, Sample ~ Alleles, value.var = 'Counts')
  # Replace NA by 0
  filtered.data[is.na(filtered.data)] <-0
  # Replace read counts by frequency
  total_reads <- apply(filtered.data[,-1],1,sum)
  filtered.data[,-1] <- sweep(filtered.data[,-1], 1,total_reads , '/')
  # Put samples as row names
  rownames(filtered.data) <- filtered.data %>% dplyr::pull(Sample)
  filtered.data <- filtered.data[,-1]
  filtered.data <- as.matrix(filtered.data)
  # Making the sample metadata 
  group_data <- groupInfo
  pal_complete <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush1","peachpuff1", "rosybrown1", "#B3e2cd", "darkslategray3","lightseagreen", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  Meta_samples <- data.frame(Samples = rownames(filtered.data),Total = total_reads, Group = groupInfo[metadata$Samples %in% rownames(filtered.data)]) 
  color<- sapply(Meta_samples$Group,FUN= function(x) {
                        ifelse(x == 'Blastocysts',"#E6A0C4", 
                        ifelse(x == 'Blastocysts_M',"#C6CDF7",
                        ifelse(x == 'Control', "#D8A499","#7294D4")))
                      })
  Meta_samples$Color <-color
  # Meta for the variants 
  pal_variantes <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush3","#b3e2cd","lightseagreen")
  type <- sapply(colnames(filtered.data), getAlleleType)
  color_v <- sapply(type, FUN = function(type) {switch(type, 
    'SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")}) 

  Meta_variants <- data.frame(Variants = colnames(filtered.data), Type = type, Color = color_v)
  # return data in a list
  return(list(matrix = filtered.data, MetadataSample = Meta_samples, MetadataVariants = Meta_variants))
}
      
#' Creating chord diagram for the variant frequency per sample
#' @param readMatrix of the locus of interest containing read counts per allele per sample
#' @param freq.filtering threshold for allele filtering (is called by the function meltedAlleleInfoFreqFiltered)
#' @param group character string designating which group to use : "Blastocysts", "Blastocysts_M", "Parternotes_M"
#'@return a chord diagram plot
VizCirclizeByGroupFreqfiltered <- function(readMatrix, freq.threshold = 0.0625, show_no_variants = FALSE, scaled = FALSE, group = "Blastocysts"){
  
  # formatting the data
  circ_data <- getDTcirclizeByGroupFreqFiltering(readMatrix, freq.threshold = freq.threshold, group = group )
  ## getting the data 
  m <- circ_data[[1]]
  if(!is.numeric(m)){
    stop(sprintf("Only one allele present: %s!",colnames(m)[[2]]))
  }
  meta_samples <- circ_data[[2]]
  meta_variants <- circ_data[[3]]
  # organizing and filtering data
  ## ordering mutations
  order_m <- order(meta_variants$Type)
  m <- m[,order_m]
  if(ncol(m) == 2){
    m2 <- t(m) 
  } else {m2 <- m}
  ## set the no variant presence
  if(show_no_variants == FALSE) {
    no_v <- which(colnames(m) == "no variant")
    m[,no_v] <- 0.001
  } 
  
  # color set up
  grid.col <- meta_samples$Color[match(rownames(m), meta_samples$Samples)]
  grid.col <- c(grid.col, as.character(meta_variants$Color[match(colnames(m),meta_variants$VariantsName)]))
  col_grad <-circlize::colorRamp2(range(m), c("darkslategray3","#FD6467")) 
  col_grad_m <- col_grad(m)
  col_grad_m[m <= freq.threshold] <- "#00000000" 
  if(ncol(m) == 2){
    col_grad_m <-t (col_grad_m)
    }
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_mutations <- c('SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")
  
  # initiate plot
  circlize::circos.clear()
  par(mar = rep(0, 4), cex=.75)
  circlize::circos.par(start.degree = ifelse(ncol(m)==2, 90,-90))
  
  circlize::chordDiagram(m2, directional = ifelse(ncol(m)==2, -1,1), 
                         transparency = 0,
                         grid.col = grid.col,
                         link.sort = TRUE,
                         link.decreasing = TRUE,
                         symmetric = FALSE,
                         diffHeight = 0,
                         link.lwd = 1,
                         col = col_grad_m,
                         scale = TRUE,
                         annotationTrack = NULL,
                         preAllocateTracks = list(
                           list(track.height = 0.03, track.margin = c(0, 0)),
                           list(track.height = 2.25 * circos.par("track.height")),
                           list(track.height = 0.02, track.margin = c(0, 0)),
                           list(track.height = 0.02, track.margin = c(0, 0)),
                           list(track.height = 0.02, track.margin = c(0, 0))
                         ),
                         big.gap = 30)
  
  for (sample in row.names(m)){
    circlize::highlight.sector(sample,
                               track.index = 2,
                               text = paste(sample, meta_samples$Total[ meta_samples$Samples == sample], sep = "    "),
                               cex = 1,
                               col = NA, 
                               border = NA, 
                               facing = "reverse.clockwise",
                               niceFacing = FALSE
    )
  }
  
  for (allele in colnames(m)){
    circlize::highlight.sector(allele,
                               track.index = 2,
                               text = allele,
                               cex = 1,
                               col = NA, 
                               border = NA, 
                               facing = "clockwise",
                               niceFacing = TRUE
    )
  }
  
  for (i in c(which(levels(meta_samples$Group) %in% meta_samples$Group))){
    Sample_group <- levels(meta_samples$Group)[i]
    ind <- which(meta_samples$Group == Sample_group)
    circlize::highlight.sector(rownames(m)[ind],
                               track.index = 2,
                               text = Sample_group, 
                               col = NA, 
                               border = pal_groups[i], 
                               facing = "bending.inside",
                               niceFacing = FALSE,
                               text.vjust = "40mm",
                               cex = 1
    )
    circlize::highlight.sector(rownames(m)[ind], track.index = c(3,4,5), col = pal_groups[i])
  }
  
  for (i in c(which(levels(meta_variants$Type) %in% meta_variants$Type))){
    Mut_type <- levels(meta_variants$Type)[i]
    ind <- which(meta_variants$Type[order(meta_variants$Type)] == Mut_type)
    circlize::highlight.sector(colnames(m)[ind],
                               track.index = 2,
                               text = ifelse(Mut_type %in% c("Other", "Reference"), "",Mut_type),
                               col = NA, 
                               border = pal_mutations[Mut_type], 
                               facing = "bending.inside",
                               niceFacing = FALSE,
                               text.vjust = "40mm",
                               cex = 1
    )
    circlize::highlight.sector(colnames(m)[ind], track.index = c(3,4,5), col = pal_mutations[Mut_type])
  }
}

#' Create a percent stack barplot of the allele frequencies by sample.
#' @param readMatrix is a matrix of sample by allele read count 
#' @param group is a string indicating the group to be used ("Blastocysts","Blastocysts_M", "Partenotes_M" and "All")
#' @return a ggplot barplot of frequencies per sample
Barplot.freq<- function(readMatrix, group = 'All', freq.threshold = 0.0625) {
  ## getting alleles with frequencies greater than 6.25%
  melted.filtered.data <- meltedAlleleInfoFreqFiltered(readMatrix,freq.threshold)
  
  #Calculating frequency after the background is removed
  Freq_after_filter <- lapply(unique(melted.filtered.data$Sample),FUN = function(x) {
    sample <- x
    counts <- sum(readMatrix[melted.filtered.data$Allele[melted.filtered.data$Sample == sample], x])
    total_reads <- sum(readMatrix[,x])
    Freq_after_filter <- total_reads*melted.filtered.data$Frequency[melted.filtered.data$Sample == sample]/counts
  })
  melted.filtered.data$Freq_after_filter <- unlist(Freq_after_filter)
  #Selecting group
  if(group == "Blastocysts"){
    melted.filtered.data <- melted.filtered.data[-grep("M",melted.filtered.data$Sample),]
  }
  if(group == "Blastocysts_M"){
    melted.filtered.data <- melted.filtered.data[grep("BLM",melted.filtered.data$Sample),]
  }
  if(group == "Partenotes_M"){
    melted.filtered.data <- melted.filtered.data[grep("P",melted.filtered.data$Sample),]
  }
  # reordering the samples
  melted.filtered.data$Sample<-factor(melted.filtered.data$Sample, levels = rev(levels(melted.filtered.data$Sample)))
  # Put frequencies as percentages
  melted.filtered.data$Percent <- signif(melted.filtered.data$Freq_after_filter,2)*100
  # plot
  p =  ggplot(melted.filtered.data, aes(x= Sample, y = Percent, fill = Alleles)) + 
  geom_bar(position = "fill", stat =  "identity")+
  scale_y_continuous(labels = scales::percent_format()) + geom_text(data = melted.filtered.data, aes(x= Sample, y = Freq_after_filter, label = paste0(Percent,"%")), position = position_stack(vjust = 0.5 )) + theme_pubr() + coord_flip()
  return(p)
  }

#' This function creates a list informing the position of complementary to the guide in the target/off-target sites
#' @param references is a DNAstring set made from the fasta files of the reference sequences
#' @param locimetadata is a data.frame object containing the information of the loci
#' @param gd2 is a GRange file containing the loci in the crisprset and amp_ranges with the names matching crisprset
#' @return a list of the  target sequence in the target/off-target sites, target location (grange) in respect to the reference fasta file, and the pam location, position of the first nucleotide in the PAM in respect to the reference file
Target_position <- function(gd,locimetadata,references) {
  # extracting target location an pam direction
  strand <- as.character(strand(gd))
  target_seq <- lapply(seq_along(gd$names),FUN = function(i) {
    name <- gd$names[i]
    target_seq <-locimetadata$Seq_loci[locimetadata$ID == name]
    loc_PAM <- unlist(gregexpr('PAM',target_seq)[1])
    target_seq <- gsub("PAM","", target_seq)
    return(list(target_seq = target_seq, PAM = loc_PAM))
  })
  names(target_seq) <- gd$names
  
  # extracting references strings from DnaStringSet
  references_seq <- lapply(seq_along(gd$names), FUN = function(i){
    references_seq <- toString(references[[i]])
  })
  # extracting target and PAM location
  targets_loc <-lapply(seq_along(gd$names), FUN = function(i){
    if(strand[i] == "-") {
      ref_seq <- reverseComplement(references[i])
      ref_seq <- toString(ref_seq)
      target_region <- unlist(gregexpr(target_seq[[i]][1],ref_seq))
    } else {
    target_region <- unlist(gregexpr(target_seq[[i]][1],references_seq[[i]][1]))
    }
    if(target_seq[[i]][2] == 1) {
      PAM_region <- target_region
      target_region <- PAM_region + 3
      target.loc <- target_region + 3
      target_region <- IRanges::IRanges(start = target_region, end = target_region + nchar(target_seq[[i]][1]))
    } else {
      PAM_region <- target_region + nchar(target_seq[[i]][1]) - 3
      target.loc <- PAM_region - 5
      target_region <- IRanges::IRanges(start = target_region, end = (target_region + nchar(target_seq[[i]][1]) - 4))
    }
    if(strand[i] == "-") {
      target.loc <- nchar(references_seq[[i]][1]) - target.loc
      PAM_region <- nchar(references_seq[[i]][1]) - PAM_region - 1
      target_region <- IRanges::IRanges( start =  (nchar(references_seq[[i]][1]) - end(target_region)+ 1), end = (nchar(references_seq[[i]][1]) - start(target_region) +1))
    }
    target_seq <- unlist(target_seq[[i]][1])
    return(list(target_seq = target_seq, target_loc = target.loc, PAM = PAM_region, target_region = target_region))
  })
  names(targets_loc) <- gd$names
  return(targets_loc)
}

#' Get the number of alleles in the sample filtered with specific frequency greater then threshold 
#' @param readmatrix a matrix with number of reads per allele per sample
#' @param sample character string with the name of the sample as it is shown in readmatrix
#' @return plots of conservation in number of reads or percentage 
plotFreqConSample <- function(readmatrix, sample ="BL01", freq.threshold = seq(0.01, 0.20, by = 0.01)){
  #getting the total alleles before filtering
  col<- readmatrix[,sample]
  #percentage of unique alleles conserved (kept) when using the threshold
  col<- col[col > 0]
  
  actualDataConservation <- sapply(freq.threshold, function(f) length(frequencyFilterPerSample(readmatrix, freq.threshold = f)[[sample]]) / length(col))
  #number of alleles kept when using the threshold
  allelesConserved <- sapply(freq.threshold, function(f) length(frequencyFilterPerSample(readmatrix, freq.threshold = f)[[sample]]))
  df <- data.frame(freq.threshold = freq.threshold, dataConservation = actualDataConservation, alleleConservation = allelesConserved)
  df_m <-reshape2::melt(df, id.vars = "freq.threshold")
  ggplot(df_m, aes(x = freq.threshold, y = value, group = variable)) +
    geom_point() +
    facet_grid(variable ~ ., scales = "free", labeller = labeller("dataConservation", "alleleConservation"), switch = "y") + geom_vline(xintercept = 0.0625, linetype = "dashed", color = "red" ) +
  ggpubr::theme_pubr() + labs(title = "Allele conservation afeter frequency filtering", x = "Frequency threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5), panel.spacing = unit(2,"lines"), axis.title.x = element_text(margin = margin(t = 100, r = 0, b = 0, l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 100, b = 0, l = 0)))
}

#' Get the number of alleles in the sample filtered with specific frequency greater then threshold in different regions of the gene 
#' @param MatrixList a list of matrix ( upstream, target and downstream) with number of reads per allele per sample 
#' @param sample character string with the name of the sample as it is shown in readmatrix
#' @return plots of conservation in number of reads or percentage 
plotFreqConSample_error_ev <- function(MatrixList, sample ="BL01", freq.threshold = seq(0.01, 0.20, by = 0.01)){
  # getting the data for the plots
  list.plot <-lapply(MatrixList, function(m) plotFreqConSample(readmatrix = m, sample = sample))
  list.data<- lapply(seq_along(list.plot), function(i) {
    x <- list.plot[[i]]$data
    x$Region <- rep(names(list.plot)[i], nrow(x))
    x
  })
  #combinig data from three regions
  data <- do.call(rbind,list.data)
  # plotting
  title <- sprintf("Allele conservation after frequency filtering in %s", sample)
  ggplot(data, aes(x = freq.threshold, y = value, group = variable, color = Region)) +
    geom_point(size = 6) +
    facet_grid(variable ~ ., scales = "free", labeller = labeller("dataConservation", "alleleConservation"), switch = "y") + geom_vline(xintercept = 0.0625, linetype = "dashed", color = "red" ) + 
    scale_x_continuous(breaks = freq.threshold[c(TRUE,FALSE)]) +
    ggpubr::theme_pubr() + labs(title = title, x = "Frequency threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5), panel.spacing = unit(2,"lines"), axis.title.x = element_text(margin = margin(t = 100, r = 0, b = 0, l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 100, b = 0, l = 0)))
                                                                                                                                                 
}

#' Get the number of alleles in the sample filtered with a allele conservation trheshold in different regions of the gene 
#' @param MatrixList a list of matrix ( upstream, target and downstream) with number of reads per allele per sample 
#' @param sample character string with the name of the sample as it is shown in readmatrix
#' @return plots of conservation in number of reads or percentage 
plotConThrSample_error_ev <- function(MatrixList, sample ="BL01", thresholds = seq(0.8, 0.99, by = 0.01)){
  # geting colunm
  col <- lapply(MatrixList, function(m) m[,sample])
  col <- lapply(col, function(cc) cc[cc>0])
  names(col) <- names(MatrixList)
  # getting the data for the plot
  list.plot <-lapply(col, function(c) plotConservationSample(col = c))
  list.data<- lapply(seq_along(list.plot), function(i) {
    x <- list.plot[[i]]$data
    x$Region <- rep(names(list.plot)[i], nrow(x))
    x
  })
  #combinig data from three regions
  data <- do.call(rbind,list.data)
  # plotting
  title <- sprintf("Allele conservation after filtering in %s", sample)
  ggplot(data, aes(x = thresholds, y = value, group = variable, color = Region)) +
    geom_point(size =6) +
    facet_grid(variable ~ ., scales = "free", labeller = labeller("dataConservation", "alleleConservation"), switch = "y") + geom_vline(xintercept = 0.9, linetype = "dashed", color = "red" ) + 
    scale_x_continuous(breaks = thresholds[c(TRUE,FALSE)]) +
    ggpubr::theme_pubr() + labs(title = title, x = "Threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),
   plot.subtitle = element_text(hjust = 0.5), panel.spacing = unit(2,"lines"), axis.title.x = element_text(margin = margin(t = 100, r = 0, b = 0, l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 100, b = 0, l = 0)))
 }

#' Plot the allele conservation across all samples given threshold in all regions of the gene
#' @param listMatrix a list of matrix with number of reads per allele per sample
#' @return plots of conservation in number of reads or percentage
plotConservation_error_ev <- function(listMatrix, thresholds = seq(0.80, 0.99, by = 0.01)){
      #get data for all regions
      list.plot <- lapply(listMatrix, function(m) plotConservation(m))
      list.data<- lapply(seq_along(list.plot), function(i) {
        x <- list.plot[[i]]$data
        x$Region <- rep(names(list.plot)[i], nrow(x))
        x
      })
      #combinig data from three regions
      data <- do.call(rbind,list.data)
      data$Region <- factor(data$Region, levels=c("downstream","target","upstream"))
      
      #plots 
      
      ggplot(data, aes(x = thresholds, y = value, group = variable, color = Region)) +
        geom_point(size=6) +
        scale_x_continuous(breaks = thresholds[c(TRUE,FALSE)]) +
        facet_grid(variable ~ ., scales = "free", labeller = label_parsed, switch = "y") + geom_vline(xintercept = 0.9, linetype = "dashed", color = "red" ) +
        ggpubr::theme_pubr() + labs(title = "Allele conservation per read-threshold", x = "Threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5))
      
    }

    
#' Plot the allele conservation across all samples given threshold in all regions of the gene
#' @param listMatrix a list of matrix with number of reads per allele per sample
#' @return plots of conservation in number of reads or percentage
plotFreqCon_error_ev <- function(listMatrix,  freq.threshold = seq(0.01, 0.20, by = 0.01)){
      #get data for all regions all thresholds
      n_alleles <- lapply(listMatrix, function(m) sapply(freq.threshold, function(f) length(frequencyFilterAllSamples(m, freq.threshold = f))))
      per_alleles <- lapply(n_alleles, function(v) (v/v[1])*100)
      list.data <- lapply(seq_along(n_alleles), function(i){
        x <- data.frame(freq.threshold = freq.threshold,dataConservation = per_alleles[[i]], alleleConservation = n_alleles[[i]],Region = names(n_alleles)[i] )
        x <- reshape2::melt(x,id.vars= c("freq.threshold", "Region") )
        x 
      })
      
      #creating data.frame
      data <- do.call(rbind,list.data)
      data$Region <- factor(data$Region, levels=c("downstream","target","upstream"))

      #plots 
      
      ggplot(data, aes(x = freq.threshold, y = value, group = variable, color = Region)) +
        geom_point(size=6) +
        scale_x_continuous(breaks = freq.threshold[c(TRUE,FALSE)]) +
        facet_grid(variable ~ ., scales = "free", labeller = label_parsed, switch = "y") + geom_vline(xintercept = 0.0625, linetype = "dashed", color = "red" ) +
        ggpubr::theme_pubr() + labs(title = "Allele conservation after filtering", x = "Freq.threshold", y = "Conservation") + theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
        }


# Function to get number of alleles when minimal (max = 2) number of alleles is achieved in control regions.
#' @param MatrixList is a list of readmatrix of counts per allele per sample for upstream, downstream, and target regions
#' @param Sample character string with the  sample name as shown on readmatrix
#' @return list of alleles present in each of the regions with frequency greater than frequency when background is removed from control regions 
filt_up_down <- function(MatrixList, Sample) {
  ## Setting freq threshold to up 30%
  freq.threshold <- seq(0.01,0.30, by = 0.01)
  #Getting number of alleles per sample in the different trhesholds 
  alleles_list <- lapply(MatrixList, FUN= function(m){
  alleles_list <- sapply(freq.threshold, function(f) length(frequencyFilterPerSample(m, freq.threshold = f)[[Sample]]))
    })
  alleles_list <- do.call(cbind,alleles_list)
  alleles_list <- as.data.frame(alleles_list)
  alleles_list$freq.threshold <- freq.threshold
  ## Getting threshold where the minimal number of alleles are achieved on control regions
  cond <- alleles_list$upstream == min(alleles_list$upstream) & alleles_list$downstream == min(alleles_list$downstream)
  min_freq <- min(alleles_list$freq.threshold[cond])
  ## filtering sample by the min_freq
  alleles_sample <- lapply(MatrixList, function(m) {
    m <- Convert2Freq(m)
    alleles <- apply(m, 2 , function(mm) mm[mm > min_freq])
    alleles <-alleles[[Sample]]
  })
  names(alleles_sample)<- names(vc_all_list)
  return(alleles_sample)
}

# Function to get min_freq to be used in each sample.
#' @param MatrixList is a list of readmatrix of counts per allele per sample for upstream, downstream, and target regions
#' @param Sample character string with the  sample name as shown on readmatrix
#' @return list of alleles present in each of the regions with frequency greater than frequency when background is removed from control regions 
min_freq <- function(MatrixList, Sample) {
  ## Setting freq threshold to up 30%
  freq.threshold <- seq(0.01,0.30, by = 0.01)
  #Getting number of alleles per sample in the different trhesholds 
  alleles_list <- lapply(MatrixList, FUN= function(m){
    alleles_list <- sapply(freq.threshold, function(f) length(frequencyFilterPerSample(m, freq.threshold = f)[[Sample]]))
  })
  alleles_list <- do.call(cbind,alleles_list)
  alleles_list <- as.data.frame(alleles_list)
  alleles_list$freq.threshold <- freq.threshold
  ## Getting threshold where the minimal number of alleles are achieved on control regions
  cond <- alleles_list$upstream == min(alleles_list$upstream) & alleles_list$downstream == min(alleles_list$downstream)
  min_freq <- min(alleles_list$freq.threshold[cond])
  return(min_freq)
}

#' Create a filtered matrix for circlize data
#' @param locus of the locus of interest, it need to be one of the locus contained in the reads_vc_all_list_min_freq_filtered object
getReadMatrixMinFreq <- function(locus) {
  filtered.data <- readRDS("./reads_vc_all_list_min_freq_filtered.RDS")[[locus]]
  filtered.data <- lapply(seq_along(filtered.data), function(i) {
      d <- filtered.data[[i]]
      d <- d[,c(1,3)]
      row.names(d) <- NULL
      d$Sample <- names(filtered.data[i])
      d} )
  filtered.data <- do.call(rbind,filtered.data)
  colnames_mat <- unique(filtered.data$Sample)
  rownames_mat <- unique(filtered.data$Variants)
  mat <- sapply(seq_along(colnames_mat), function(i) {
    mat <- sapply(seq_along(rownames_mat), function(j) {
    mat <- ifelse(length(filtered.data$Reads[filtered.data$Sample == colnames_mat[i] & filtered.data$Variants == rownames_mat[j]])!= 0, filtered.data$Reads[filtered.data$Sample == colnames_mat[i] & filtered.data$Variants == rownames_mat[j]], 0)}
    )
  })
  colnames(mat) <- colnames_mat
  rownames(mat)  <- rownames_mat
  return(mat)
  
  }
#' Create a matrix for circlize data, plus metadata for samples and metadata for variants
#' @param locus of the locus of interest, it need to be one of the locus contained in the reads_vc_all_list_min_freq_filtered object
#' @param min.freq named vector of frequencies threshold (is called by the function meltedAlleleInfo)
#' @param metadata dataframe containing the group information for the samples
#' @param group character string designating which group to use : "Blastocysts", "Blastocysts_M", "Parternotes_M"
#' @return a list containing a matrix, a metadata for the blastocyst samples and a metadata for the variants 
getDTcirclizeByGroupMinFreq <- function(locus, metadata,  min.freq, group = "Blastocysts" ) {
  #Get allele information after filtering in a data frame of sample per allele
  
  fil.mat <- getReadMatrixMinFreq(locus)
  fil.mat <- fil.mat[,which(colnames(fil.mat) %in% metadata$Samples[metadata$Group %in% c(group, "Control")])] 
  # Replace NA by 0
  fil.mat[is.na(fil.mat)] <-0
  if(is.vector(fil.mat[apply(fil.mat, 1, function(x) !all(x==0)),])){
    unique_var <- rownames(fil.mat)[which(fil.mat != 0, arr.ind =  T )[1]]
   sprintf(" WARNING: There is only one allele (%s) present in the samples.",unique_var)
  }
  fil.mat <- fil.mat[apply(fil.mat[,-1], 1, function(x) !all(x==0)),]
 
  # Replace NA by 0
  fil.mat[is.na(fil.mat)] <-0
  # Put samples as row names
  
  if(!is.vector(fil.mat)) {
    fil.mat <- as.matrix(fil.mat)
  
  # Making the sample metadata 
  group_data <- metadata$Group[match(colnames(fil.mat),metadata$Samples)]
  pal_complete <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush1","peachpuff1", "rosybrown1", "#B3e2cd", "darkslategray3","lightseagreen", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  total_reads <- apply(fil.mat,2,sum)
  names(total_reads) <- colnames(fil.mat)
  Meta_samples <- data.frame(Samples = colnames(fil.mat),Total = total_reads, Group = group_data) 
  color<- sapply(Meta_samples$Group,FUN= function(x) {
    ifelse(x == 'Blastocysts',"#E6A0C4", 
           ifelse(x == 'Blastocysts_M',"#C6CDF7",
                  ifelse(x == 'Control', "#D8A499","#7294D4")))
  })
  Meta_samples$Color <-color
  # Meta for the variants 
  pal_variantes <- c("#F1BB7B", "#FD6467", "#5B1A18","lavenderblush3","#b3e2cd","lightseagreen")
  type <- sapply(rownames(fil.mat), getAlleleType)
  color_v <- sapply(type, FUN = function(type) {switch(type, 
                                                       'SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")}) 
  
  Meta_variants <- data.frame(Variants = rownames(fil.mat), Type = type, Color = color_v)
  #reducing the name of the SNV variants
  snv_v <- grep("SNV",rownames(fil.mat))
  rownames(fil.mat)[snv_v] <- paste0("SNV",seq(1:length(snv_v)))
  Meta_variants$VariantsName <- rownames(fil.mat)
  } else {
  Meta_samples <- data.frame(Samples = names(fil.mat),Total = fil.mat, Group = metadata$Group[metadata$Samples %in% names(fil.mat)]) 
  Meta_variants <- c(unique_var , getAlleleType(unique_var))
  sprintf("WARNING: There is only one allele (%s) present in the samples.",unique_var)
  }

  # return data in a list
  return(list(matrix = fil.mat, MetadataSample = Meta_samples, MetadataVariants = Meta_variants))

}
  
#' Creating chord diagram for the variant frequency per sample
#' @param locus of the locus of interest, it need to be one of the locus contained in the reads_vc_all_list_min_freq_filtered object
#' @param group character string designating which group to use : "Blastocysts", "Blastocysts_M", "Parternotes_M"
#'@return a chord diagram plot
VizCirclizeByGroupMinFreq <- function(locus, metadata,  min.freq, show_no_variants = FALSE, scaled = FALSE, group = "Blastocysts"){
  
  # formatting the data
  circ_data <- getDTcirclizeByGroupMinFreq(locus = locus, metadata = metadata, min.freq = min.freq, group = group)
  ## getting the data 
  if(!is.matrix(circ_data[[1]])){
    stop(sprintf("Only one allele present: %s!",circ_data[[3]][1]))
  }
  sum_m <- apply(circ_data[[1]],2,sum)
  m  <- sweep(circ_data[[1]], 2,sum_m , '/')
  m <- t(m)
  m
  meta_samples <- circ_data[[2]]
  meta_variants <- circ_data[[3]]
  # organizing and filtering data
  ## ordering mutations
  order_m <- order(meta_variants$Type)
  m <- m[,order_m]
  if(ncol(m) == 2){
    m2 <- t(m) 
  } else {m2 <- m}
  ## set the no variant presence
  if(show_no_variants == FALSE) {
    no_v <- which(colnames(m) == "no variant")
    m[,no_v] <- 0.001
  } 
  
  # color set up
  grid.col <- meta_samples$Color[match(rownames(m), meta_samples$Samples)]
  grid.col <- c(grid.col, as.character(meta_variants$Color[match(colnames(m),meta_variants$VariantsName)]))
  col_grad <-circlize::colorRamp2(range(m), c("darkslategray3","#FD6467")) 
  col_grad_m <- col_grad(m)
  col_grad_m[m <= 0.001] <- "#00000000" 
  if(ncol(m) == 2){
    col_grad_m <-t (col_grad_m)
    }
  pal_groups <- c("#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")
  pal_mutations <- c('SNV' = "#F1BB7B", 'Insertion' = "#FD6467", 'Deletion' = "#5B1A18", 'Mixed' = "lavenderblush3", 'Reference' = "#b3e2cd", 'Other' = "lightseagreen")
  
  # initiate plot
  circlize::circos.clear()
  par(mar = rep(0, 4), cex=.75)
  circlize::circos.par(start.degree = ifelse(ncol(m)==2, 90,-90))
  
  circlize::chordDiagram(m2, directional = ifelse(ncol(m)==2, -1,1), 
                         transparency = 0,
                         grid.col = grid.col,
                         link.sort = TRUE,
                         link.decreasing = TRUE,
                         symmetric = FALSE,
                         diffHeight = 0,
                         link.lwd = 1,
                         col = col_grad_m,
                         scale = TRUE,
                         annotationTrack = NULL,
                         preAllocateTracks = list(
                           list(track.height = 0.03, track.margin = c(0, 0)),
                           list(track.height = 2.25 * circos.par("track.height")),
                           list(track.height = 0.02, track.margin = c(0, 0)),
                           list(track.height = 0.02, track.margin = c(0, 0)),
                           list(track.height = 0.02, track.margin = c(0, 0))
                         ),
                         big.gap = 30)
  
  for (sample in row.names(m)){
    circlize::highlight.sector(sample,
                               track.index = 2,
                               text = paste(sample, meta_samples$Total[ meta_samples$Samples == sample], sep = "    "),
                               cex = 1,
                               col = NA, 
                               border = NA, 
                               facing = "reverse.clockwise",
                               niceFacing = FALSE
    )
  }
  
  for (allele in colnames(m)){
    circlize::highlight.sector(allele,
                               track.index = 2,
                               text = allele,
                               cex = 1,
                               col = NA, 
                               border = NA, 
                               facing = "clockwise",
                               niceFacing = TRUE
    )
  }
  
  for (i in c(which(levels(meta_samples$Group) %in% meta_samples$Group))){
    Sample_group <- levels(meta_samples$Group)[i]
    ind <- which(meta_samples$Group == Sample_group)
    circlize::highlight.sector(rownames(m)[ind],
                               track.index = 2,
                               text = Sample_group, 
                               col = NA, 
                               border = pal_groups[i], 
                               facing = "bending.inside",
                               niceFacing = FALSE,
                               text.vjust = "40mm",
                               cex = 1
    )
    circlize::highlight.sector(rownames(m)[ind], track.index = c(3,4,5), col = pal_groups[i])
  }
  
  for (i in c(which(levels(meta_variants$Type) %in% meta_variants$Type))){
    Mut_type <- levels(meta_variants$Type)[i]
    ind <- which(meta_variants$Type[order(meta_variants$Type)] == Mut_type)
    circlize::highlight.sector(colnames(m)[ind],
                               track.index = 2,
                               text = ifelse(Mut_type %in% c("Other", "Reference"), "",Mut_type),
                               col = NA, 
                               border = pal_mutations[Mut_type], 
                               facing = "bending.inside",
                               niceFacing = FALSE,
                               text.vjust = "40mm",
                               cex = 1
    )
    circlize::highlight.sector(colnames(m)[ind], track.index = c(3,4,5), col = pal_mutations[Mut_type])
  }
}


#' Create a percent stack barplot of the allele frequencies by sample.
#' @param locus is the locus of interest 
#' @param group is a string indicating the group to be used ("Blastocysts","Blastocysts_M", "Partenotes_M" and "All")
#' @return a ggplot barplot of frequencies per sample
Barplot.minfreq<- function(locus, group = 'All', freq.threshold = 0.0625) {
  ## getting alleles with  min frequencies 
  filtered.data <- readRDS("./reads_vc_all_list_min_freq_filtered.RDS")[[locus]]
  filtered.data <- lapply(seq_along(filtered.data), function(i) {
    d <- filtered.data[[i]]
    d <- d[,c(1,3)]
    row.names(d) <- NULL
    d$Sample <- names(filtered.data[i])
    d$Freq <- sapply(d$Reads, function(r) {r/sum(d$Reads)})
    d} )
  melted.filtered.data <- do.call(rbind,filtered.data)
  
  #Selecting group
  if(group == "Blastocysts"){
    melted.filtered.data <- melted.filtered.data[-grep("M",melted.filtered.data$Sample),]
  }
  if(group == "Blastocysts_M"){
    melted.filtered.data <- melted.filtered.data[grep("BLM",melted.filtered.data$Sample),]
  }
  if(group == "Partenotes_M"){
    melted.filtered.data <- melted.filtered.data[grep("P",melted.filtered.data$Sample),]
  }
  # reordering the samples
  melted.filtered.data$Sample <- as.factor(melted.filtered.data$Sample)
  melted.filtered.data$Sample<-factor(melted.filtered.data$Sample, levels = rev(levels(melted.filtered.data$Sample)))
  # reordering variants
  melted.filtered.data$Variants<-factor(melted.filtered.data$Variants, levels = levels(melted.filtered.data$Variants)[order(levels(melted.filtered.data$Variants), decreasing = FALSE)])
  # Put frequencies as percentages
  melted.filtered.data$Percent <- signif(melted.filtered.data$Freq,2)*100
  # plot
  p =  ggplot(melted.filtered.data, aes(x= Sample, y = Percent, fill = Variants)) + 
    geom_bar(position = "fill", stat =  "identity")+
    scale_y_continuous(labels = scales::percent_format()) + geom_text(data = melted.filtered.data, aes(x= Sample, y = Freq, label = paste0(Percent,"%")), position = position_stack(vjust = 0.5 )) + theme_pubr() + coord_flip()
  return(p)
}
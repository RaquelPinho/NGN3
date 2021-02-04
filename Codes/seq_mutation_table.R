

getSeqTablePerSample <- function(crisprset, locus = "NGN3", sample, references) {
  # Getting the alleles filtered by sample
  reads_vc_all_list <- readRDS("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R/reads_vc_all_list_min_freq_filtered.RDS")
  allele_sample <-reads_vc_all_list[[locus]][[sample]]
  allele_sample <- allele_sample %>%  arrange(desc(Reads)) 
  allele_sample$Sample <- sample 
  allele_sample$Locus <- locus
  allele_sample <- allele_sample[,c(5,1:4,6)]
  rownames(allele_sample) = NULL
  
  # Geting the sequence for the reads
  
  sample_alns <- crisprset[["crispr_runs"]][[sample]][["alns"]]
  data_al_sample <- data.frame(Allele = sample_alns@elementMetadata$allele,Seq = sample_alns@elementMetadata$seq) 
  data <- data_al_sample %>% group_by(Seq,Allele) %>% count(Seq) %>% arrange(desc(n))
  # Combining the data and the seq
  allele_sample$Seq <- data$Seq[match(allele_sample$Variants, data$Allele)]
  allele_sample$Reads_sample <- data$n[match(allele_sample$Variants, data$Allele)]
  if(length(allele_sample$Variants == "Other") !=0){
    allele_sample <- allele_sample[allele_sample$Variants != "Other",]
  }  
    # Checking the differences 
  #ref <- unique(data[which(data$Allele == "no variant"),])
  #ref <- ref$Seq[which(ref$n == max(ref$n))]
  ridx <- grep(locus, names(references))
  if(ridx == 1) {
    ref <- as.character(reverseComplement(references[ridx]))
  } else {
    ref <- as.character(references[ridx])
  }
  Mutation <- lapply(seq_along(allele_sample$Seq), function(i){
    bio_algm <- Biostrings::pairwiseAlignment(allele_sample$Seq[i],ref)
    ref_s <- as.character(bio_algm@subject)
    allele_s <- as.character(bio_algm@pattern)
    Alignement = data.frame(Sample = sample,Locus = locus, type = c("reference","allele"),Sequence = c(ref_s,allele_s), Variant = rep(as.character(allele_sample$Variants[i]),2))
    
  })
   Mutation <- do.call(rbind,Mutation)
   return(list(allele_info = allele_sample, Alignement = Mutation))
}
names_samples <- names(crisprset$crispr_runs)
tables <- lapply(names_samples,function(sample) { getSeqTablePerSample(crisprset = crisprset, locus = "NGN3",references = references, sample = sample)})
names(tables) <- names_samples
Allele_info <- lapply(tables, function(list){
  Allele_info <- list[['allele_info']]
})
 Allele_info <- do.call(rbind,Allele_info)

 Alignments <- lapply(tables, function(list){
  Alg <- list[['Alignement']]
})
 Alignements <- do.call(rbind,Alignments)
 
write.table(Allele_info,file = "./Allele_info_target_NGN3.csv", sep = "\t", row.names = F, quote = F)
write.table(Alignements,file = "./Alignments_NGN3.csv", sep = "\t", row.names = F, quote = F)
 #####################################################################################################################
##################### AFTER THIS IS BAD ##########################################################################
####################################################################################################################

Mutation <- do.call(rbind,Mutation)
Mutation <- sapply(allele_sample$Seq, function(s) {
  diffc=diag(attr(adist(ref,s4, counts = TRUE), "trafos"))
  Mutation =regmatches(ref,regexpr("[^M]",diffc))
  Mutation})
getting_mutation <- function(allele,reference) {
  ref<- reference
  mut_dia<- diag(attr(adist(ref,allele, counts = TRUE), "trafos"))
  idx <- gregexpr("[^M]",mut_dia)[[1]]
  if(idx ==-1){
    print("no mutation")
  } else {
    dt_index <-data.frame(type = unlist(strsplit(mut_dia,""))[idx], index = idx)
    ref_split <-unlist(strsplit(ref,""))
    allele_split <-unlist(strsplit(allele,""))
    mut_allele <- ref_split 
       if(length(dt_index$type == "D") > 0) {
          dels<-dt_index$index[dt_index$type == "D"]
          mut_allele[dels] <- "-"
       } 
      if(length(dt_index$type == "I") > 0) {
         ins <- dt_index$index[dt_index$type == "I"]
         ins_base <- allele_split[ins]
         names(ins_base) <-ins
       }
  }
  # Getting the sequential mutation
  is.sequential <- function(x){
         abs(diff(x)) == 1
     }  
  seq_ins <- list()
  seq_i<- is.sequential(dt_index$index)
    inc <- which(seq_i == FALSE)
    if(length(inc) > 0) {
    inc <- inc + 1
    seq_ins <- lapply(seq_along(inc), function(i){
      seq_ins <- c(inc[i],inc[i+1]-1)
      if(i == length(inc)){
        seq_ins <- c(inc[i],length(seq_i)+1)
      }
      seq_ins
      })
    if(inc[1]>=2) {
      seq_1 <- c(1,inc[1]-1)
      seq_ins <-list(seq_1, unlist(seq_ins))
       }
    }
    seq_ins_indx <- lapply(seq_ins, function(seqi){
      seq_ins_indx <- dt_index$index[seqi]
      Insert <- ifelse(identical(as.character(dt_index$type[seqi]),rep("I",length(seqi))), ins_base[as.character(seq_ins_indx)],NA)
      Insert
    }) 
    }
    
  
    
  
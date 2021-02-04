### This function create a target file for the reagion spannig the same size of the target region but at 300 bp upstream or downstream of the cut site.
## @param references DNAStringSet containing the target sequences
## @param amp_seq path to .fasta containing the amplicon sequences
## @param gd GRange object containing the target sequences
## @param target_loc list containing guide binding site, target position (cut) in relation to references regions and the PAM position
## @return error_control DNAStringSet containing the regions and 
## @return gd_error_control gd containing the range in the amplicon
Error_evaluation_region <- function(references, amp_seq, gd ){
  ## Geting sequence from fasta file
  amplicons <- Biostrings::readDNAStringSet(amp_seq)
  amp_seq <- lapply(seq_along(amplicons),FUN= function(i) {
    toString(amplicons[i])
  })
  names(amp_seq) <- names(amplicons)
  
  ## matching the order of gd and amp_seq
  locus <- sapply(strsplit(names(amp_seq),":"), "[", 1)
  order_names <- sapply(seq_along(locus), function(i){grep(locus[i],gd$names)})
  amp_seq <- amp_seq[order_names]
  locus <- locus[order_names]
  ## Getting position of the cut site for each locus
  cut <- sapply(locus, function(locus) {
    cut <- target_loc[[locus]]$target_loc
  })
  
  ## Getting the start and the end position as well as the sequence of the control error region upstream the cut
  start_up <-start(gd) + cut - 300
  end_up <- start_up + width(gd) - 1
  name <- paste0(names(amp_seq),":",start_up,"-",end_up)
  name <- gsub(" ","",name)
  seq <- substr(amp_seq,start_up,end_up)
  names(seq) <- name
  control_up < seq
  ## Creating the genomic range objects containing the informaion on gd for those control sequences
  
  gd_up <- gd
  start(gd_up) <- start_up
  end(gd_up) <- end_up
  
  ## Getting the start and the end position as well as the sequence of the control error region downstream the cut
  start_down <-start(gd) + cut + 300
  end_down <- start_down + width(gd) - 1
  name <- paste0(names(amp_seq),":",start_down,"-",end_down)
  name <- gsub(" ","",name)
  seq <- substr(amp_seq,start_down,end_down)
  names(seq) <- name
  control_down <- seq
  
  ## Creating the genomic range objects containing the informaion on gd for those control sequences
  gd_down <- gd
  end(gd_down) <- end_down 
  start(gd_down) <- start_down 
  
  ## return
  return(list(control_up =control_up,gd_up=gd_up,control_down=control_down,gd_down=gd_down))
  
}
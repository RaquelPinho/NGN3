# Get to know the working directory
getwd()
setwd("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R/")
# loading CrisprVariants GenomicAlignemts
library(CrispRVariants)
library(GenomicAlignments)
library(rBLAST)
library(ggmsa)
# loading the crisprset object
load("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R/Data/CrisprSets_log_del.RData")
# getting the crispr_set for NGN3
crispr_set <- All_sets_long_del[[1]] 

# use the get chimeras function - it is pro sample 
chimeras_bl04 <- getChimeras(crispr_set, sample = "BL04")
chimeras_bl04

# To isolate the sequence we can tranform the GAlignment object into GRange object
gr_chimeras_bl04 <- as(chimeras_bl04, "GRanges")
## then get the DNAStringSet 
dnas_chimeras_bl04 <- gr_chimeras_bl04$seq
## then I can use the 
sequences <- sapply(dnas_chimeras_bl04, as.character)
head(sequences)

# But first we can see how many unique chimeras theres was in this sample and what are the alignemet from GenomicAlignment
# For that I use the cigar descriptive.
#. M, X, = Sequence match or mismatch
#. I Insertion to the reference
#. D Deletion from the reference
#. N Skipped region from the reference
#. S Soft clip on the read - maybe real chimeras
#. H Hard clip on the read - maybe duplication (is matchiing the reference sequence, in this case the amplicon sequence)
#. P Silent deletion from the padded reference

cigar_bl04 <- chimeras_bl04@cigar
cigar_bl04 <- sort(table(cigar_bl04),decreasing = T)
most_frequent_chimera <- grep(names(cigar_bl04[1]), chimeras_bl04@cigar)
seq_mf_chimera <- sequences[most_frequent_chimera][1]

# getting the DNA string of the most frequent chimera

DNASt_mf_chimera <- 

# I will try to use the ggmsa to align the read of the most frequent quimera to the target and off target 

# Fist I will load the .fasta files of the reference genome for the target / off-target region
ref_seq_long <- "./Data/ref_seq.fasta"


# we can then use blast to align the chimeras sequence to the reference sequence and figure it out if they are real duplicate or if we can align it to other parte of the genome

# Get to know the working directory
getwd()
setwd("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R/")
# loading CrisprVariants GenomicAlignemts
library(CrispRVariants)
library(GenomicAlignments)
library(rBLAST)
library(ggmsa)
library(annotate)
library(Biostrings)
# loading the crisprset object
load("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R/Data/CrisprSets_log_del.RData")
# getting the crispr_set for NGN3
crispr_set <- All_sets_long_del[[1]] 

# use the get chimeras function - it is pro sample 
chimeras_bl04 <- getChimeras(crispr_set, sample = "BL04")
chimeras_bl04



################################################################################################################################
# I got the function that created the crispr run to check how the chimeras are selected
# to get the function I used : 'getMethod("readsToTarget",signature = c(reads = "GAlignments", target= "GRanges"))
# because when you use 'showMethods(readsToTarget)' there are three different signatures, depending on the class of the parameters
# https://rdrr.io/github/markrobinsonuzh/CrispRVariants/src/R/initialisers.R
# https://rdrr.io/github/markrobinsonuzh/CrispRVariants/src/R/findChimeras.R
##################################################################################

#  Find chimeric reads, assuming that the GAlignments object
#' does not contain multimapping reads. That is, read names that appear
#' more than ones in the file are considered chimeras.  Chimeric reads
#' are reads that cannot be mapped as a single, linear alignment.  Reads
#' from structual rearrangements such as inversions can be mapped as chimeras.

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

# So the reads that have duplicates are chimeras. we can combine the duplicate reads and try to blast 
# To check the most frequent combination of cigar labels among the duplicated reads:

View(crispr_set[[".->crispr_runs"]][["BL04"]][[".->chimera_combs"]])
mf_comb_chimeras <- head(crispr_set[[".->crispr_runs"]][["BL04"]][[".->chimera_combs"]])

## merging the two alignments of the same read for the most frequent combination


# getting the DNA string of the most frequent chimera

DNASt_mf_chimera <- dnas_chimeras_bl04[most_frequent_chimera[1]]

# I will try to use the ggmsa to align the read of the most frequent quimera to the target and off target 

# Fist I will load the .fasta files of the reference genome for the target / off-target region
# saving the path
ref_seq_long <- "./Data/ref_seq.fasta"
# import the sequences as DNAstring set
ref_seq_long <- readDNAStringSet(ref_seq_long, format = "fasta" )
# transform strings into character vectors
ref_seq_long_cha <- sapply(ref_seq_long, as.character)
# we can then use blast to align the chimeras sequence to the reference sequence and figure it out if they are real duplicate or if we can align it to other parte of the genome

### Using blast with the annotate package
# http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html
# this bioconductor package has a function `blastSequences` that sends a query for ncbi site an retrn the results.
# I will try to use the sequence of most frequent quimera to test it:

blast_mf <- blastSequences(seq_mf_chimera, database = "nr")
head(blast_mf)
blast_mf_2<- blastSequences(seq_mf_chimera, database = "nr", hitListSize = 2)
blast_mf
blast_df <- blastSequences(seq_mf_chimera, database = "nr", as = "data.frame")

bl_df <- blastSequences(seq_mf_chimera, database = "nr", as = "data.frame")
### Using the pairwiseAlignment function from BioStrings

# To use this function we need to use the DNAString/ DNAStringSet objects
# geting the complementary reverse of the  reference sequence because the pairwise alignment only works if both sequences are on the same direction.

ref2 <- readDNAStringSet("./Data/ref_seq_cr.fasta", "fasta")
as.character(ref2[[1]])

# running paorwise alignment

Biostrings::pairwiseAlignment(DNASt_mf_chimera,ref_seq_long[[1]])
Biostrings::pairwiseAlignment(DNASt_mf_chimera,reverseComplement(ref2[[1]]))
Biostrings::pairwiseAlignment(reverseComplement(ref2[[1]]),ref_seq_long[[1]])

# To do the pairwise alignemet with the sequences of the off-target sites

AlignTest <- lapply(ref_seq_long, FUN = function(ref) {x = pairwiseAlignment(DNASt_mf_chimera, ref)})

# This is an object PairwiseAlignmentsSingleSubject and there are function to used in this type of the object

as.matrix(AlignTest)
nchar(AlignTest)


#DECIPHER





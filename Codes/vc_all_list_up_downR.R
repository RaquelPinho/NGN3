# Loading reference, GRange objects, trancript database and metadatas
## I saved the fasta files with the reference sequences and the GRange objects of the target / off target regions in a .RData file
load("./NGN3_cr_Robj.RData")
## where gd is the GRange of the target regions,
## amp_ranges is the GRanges of  the amplicon regions with seqlevels as the ncbi id in the transcript database
## the references is a DNAStringSet, with the references sequences of the loci of target and off-target sites, 
## locimetadata is a df with the information of the loci,
## metadata is the metadata of the samples with the information of file paths and group information
# Loading the data for the up/downstream regions that were created with the D:/Raquel/Desktop/Post_doc/Codes/Error_evaluation_region.r script function
load("./Up_downstream_objects.RData")
## where control_up/control_down is a list of sequences for the regions 300bp up/downstream of the cut site in the target and off target regions
## where gd_up/gd_down is the GRange objects of the regions in control_up/control_down
# target locus and sequence of guide binding region
## Getting crispr_sets for the taget region and for the 300 bp upstream and downstream region of the cut
load("./crispr_sets_up_down.RData") 
## where crispr_set_target is the crispr set for the target region in all samples
## where crispr_set_up is the crispr set for the upstream region in all samples
## where crispr_set_down is the crispr set for the downstream region in all samples
target_loc <- Target_position(gd,locimetadata, references)
# creating matrix of read counts per alleles 
# for the target region
vc_all_target <- variantCounts(crispr_set_target)
# for the uspstream region
vc_all_up <- variantCounts(crispr_set_up) 
# for the downstream region
vc_all_down <- variantCounts(crispr_set_down)

## making list of the variant counts 
vc_all_list_up_down <- list(downstream = vc_all_down, target = vc_all_target, upstream = vc_all_up)
# save vc_all_list and target_loc object
save(vc_all_list_up_down, target_loc, file = "./vc_all_list_up_down.RData")

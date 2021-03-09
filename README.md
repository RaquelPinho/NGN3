# Efficiency of CRISPR-Cas9 targeting neurogenin-3 gene and off-target effects in porcine model.

## Investigators

Pablo Ross, PhD (?):  Supervision, Project adminstration, Funding acquisition.

Insung Park (UC Davis): Investigation, data collection, writing, review.

Sergio Navarro-Serna (University of Murcia, Department of Physiology) Investigation, data collection, writing, review.

Raquel Pinho, PhD (UC Davis), rdpinho@ucdavis.edu: Data analysis and visualization (PacBio)


## Project Summary

In this project we use the CRISPR/Cas9 system targeting the porcine neurogenin-3 gene in blastocysts and parthenotes to analyze the mutation and cleavage efficiency at the target region, as well as the identification of multiple variant formation (mosaicism) through high-throughput long read amplicon sequencing using PacBio technology. 
This project will also analyze seven off-target sites using the PacBio technology in long (more than 1 kb) amplicons. 

### Samples: 

#### For analysis of mosaicism:

31 IVF blastocysts (group = BLM), 28 parthenotes (group = PTM) , one control WT (parthenote)(group = PTCL).

#### For off-target analysis:

15 blastocysts (BL) were collected after IVF and zygote electroporation. 
One blastocyst collected after IVF and not electroporated was used as a control (group = BLCL).


### Methods:

#### Mosaicism analysis:

Blastocysts/parthenotes (76) were lysed and nested PCR for the target region, (primers for a shorter amplicon was used for sanger seq (730 ish?) and a long one for pacbio, samples for mosaicism were only tested for the target region). 
These samples were sequenced using the chain termination method (Sanger sequencing), processed at the (Genewiz). The long-read single-molecule real time sequencing (SMRT, PacBio, Sequel II sequencing) methodology for the target region was also perfomed at the (Genome center, UC Davis).

#### For off-target analysis:

Blastocysts (n = 16) were lysed using lysis buffer (Epicentre, cat# QE09050), and nested PCR reactions for target and 7 off-target sites (BreakingCas). The amplicons were sent for sequencing using the SMRT technology. 

Samples were divided in 4 groups (BL - blastocysts, BLM - blatocysts mosaic, PTM - partenotes mosaic, CL - wildtype).

Primers, adapters and detailed protocol described in the file: ./Data/Metadata/SN_primers_barcodes_NGN3.pdf

#### Data Analysis 

##### Mosaicism

74 samples divided in 3 groups (BL - blastocysts, BLM - blatocysts mosaic, PTM - partenotes mosaic) including two controls not injected (PTCL, BLCL). The target region will be analysed in all samples and the off target only in the BL group. The sequences were obtained through PacBio technology and the results will be compared to the sanger sequence of the same samples (BLM, PTM).


# Efficiency of CRISPR-Cas9 targeting neurogenin-3 gene and off-target effects in porcine model.

## Investigators

Pablo Ross, PhD (company?):  Supervision, Project adminstration, Funding acquisition.
Insung Park (UC Davis): Investigation, data collection, writing, review.
Sergio Navarro (Institution?) Investigation, data collection, writing, review.
Raquel Pinho, PhD (UC Davis), rdpinho@ucdavis.edu: Data analysis and visualization (PACBIO)

## Project Summary

In this project we use the CRISPR/Cas9 system targeting the porcine neurogenin-3 gene in blastocysts and parthenotes to analyze the mutation and cleavage efficiency at the target region, as well as the identification of multiple variant formation (mosaicism) through high-throughput long read amplicon sequencing using PACBIO technology. (ask Insung which one ?). This project will also analyze seven off-target sites using the PACBIO technology in long (more than 1 kb) amplicons. 

### Samples: 

#### For analysis of mosaicism:

31 IVF blastocysts (group = BLM), 28 parthenotes (group = PTM) , one control wt (IVF or parthenote?) (group = PTCL)

#### For off-target analysis:

15 blastocysts (BL) were collected after IVF and zygote electroporation. One blastocyst collected after IVF and not elecrtroporated was used as a control (group = BLCL).?


### Methods:

#### Mosaicism analysis:

Blastocysts/parthenotes (76) were lysed and nested PCR for the target region (which primer ? only target? which amplicon seq?). These samples were sequenced using the chain termination method (Sanger sequencing), processed at the (UC Davis sequencing facility?). The long-read single-molecule real time sequencing (SMRT, PacBio) methodology for the target region was also perfomed at the (Genome center, UC Davis?).

#### For off-target analysis:

Blastocysts (n = 16) were lysed using lysis buffer (Epicentre, cat# QE09050), and nested PCR reactions for target and 7 off-target sites (How was the off-target sites decided?).The amplicons were sent for sequencing using the SMRT technology. 

Samples were divided in 3 groups (BL - blastocysts, BLM - blatocysts mosaic, PTM - partenotes mosaic).

Primers, adapters and detailed protocol described in the file: ./Data/Metadata/SN_primers_barcodes_NGN3.pdf

#### Data Analysis 

##### Mosaicism

74 samples divided in 3 groups (BL - blastocysts, BLM - blatocysts mosaic, PTM - partenotes mosaic) including two controls not injected (PTCL, BLCL). The target region will be analysed in all samples and the off target only in the BL group. The sequences were obtained through PacBio technology and  the results will be compared to the sanger sequence of the same samples (BLM, PTM),


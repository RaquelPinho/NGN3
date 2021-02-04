---
title: "Off-target evaluation of CRISPRs/Cas9 system targeting gene NGN3 in porcine embryos"
author: "Raquel"
date: "4/29/2020"
output: 
    ioslides_presentation:
    css: temp1.css
runtime: shiny
---

```{r setup, include = FALSE }
knitr::opts_knit$set(root.dir = "D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R")
knitr::opts_chunk$set(echo = FALSE)

```

## NGN3 variants analysis

The analysis of the variantes of off target loci from NGN3 guides in porcines.The guide used was CCCGCGCAGCGCAUCCAACG.

### Libraries used  

I will use mainly the CrispRVariants package, for plotting ggplot2 and for importing the data rtracklayer.

```{r libraries, echo = FALSE, include=FALSE}
library("ggplot2")
library("CrispRVariants")
library("Rsamtools")
library("gdata")
library("rtracklayer")
library("gridExtra")
library("shiny")
library("plotly")

```

### Samples and metadata 

Samples from blastocysts derivides from IVF (BL) or parthenogenesis (PT) were analyzed with targeted sequencing for the target region (NGN3) and 7 predicted off-target sites (OF) 


```{r set_samples, include=F, echo=FALSE}
# Setting up the working directory
setwd("D:/Raquel/Desktop/Post_doc/NGN3/NGN3_R")

# Read samples files, each .bam file are the reads for all (target/offtarget) sites for one sample.
readfiles<-list.files("./Bam")

# Creating metadata table for samples:
Samples <- gsub(".sort.bam", "", readfiles)
Group <- vector(mode = "character", length = length(readfiles))
Files <- file.path("./Bam",readfiles)
Group[grep("PTM",readfiles)] <- "Parthenotes_M" 
Group[grep("BL",readfiles)] <- "Blastocysts" 
Group[grep("BLM",readfiles)] <- "Blastocysts_M" 
Group[grep("CL",readfiles)] <- "Control" 
metadata <- data.frame (Samples, Files = Files, Group)
metadata
# Loading references and bed file
## I saved the fasta files with the reference sequences and the GRange objects of the target / off target regions in a .RData file
load("./NGN3.RData")
# where gd2 is the GRange of the target regions,
# the references2 is a DNAStringSet, with the sequences, 
# crispr_set3_all, a crispr set with all the samples but targeting only the NGN3 locus
# locimetadata is a df with the information of the loci
# txlocigff is a TxDb of the regions containing the target and off target regions in the pig genome GCF_000003025.6_Sscrofa11.1_genomic


```
## Loci information

```{r locimetadata}
ui <- fluidPage(
    # define title of the panel
    titlePanel("Loci metadata"),
    # define input
    sidebarLayout( 
        sidebarPanel( 
            selectInput( inputId = "Loci", label = "Choose target / off-target (OF) region", choices = locimetadata$ID)),
        # define output 
        mainPanel(
            tableOutput("table")
        )))
# define server logic options
server <- function(input, output) {
    tableInput <-  reactive({
        switch(input$`Loci metadata`, locimetadata$ID)
    })
    output$table <-renderTable({
        table <- locimetadata[locimetadata$ID == input$Loci, -11]
        table <- t(as.matrix(table))
        table <- cbind(colnames(locimetadata)[-11],table)
        colnames(table) <- c("Feature", "Loci")
        table
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

```

## Variant frequency at the target site by group 

- The allele and mutation frequency was calculated using the CrispRVariants 
- The GCF_000003025.6_Sscrofa11.1_genomic database downloaded from the NCBI. 
- 1 control (CL) sample for the blastocysts and one control for the parthenotes 
- Displaying alleles with more the 10% frequency
- Samples were diveded in 4 groups: Blastocysts, Blastocysts_M (potencially mosaic), Partenotes_M ( potencially mosaic), Controls. The controls samples are also going to be presented when a group that is not control is chosen.

## 

```{r freq_group, out.height = "110%", out.width = "110%"}
vc <- variantCounts(crispr_set3_all)
ui2 <- fluidPage(
    # define title of the panel
    titlePanel("Group"),
    # define input
    sidebarLayout( 
        sidebarPanel( 
            selectInput( inputId = "Group", label = "Choose group", choices = unique(metadata$Group))),
        # define output 
        mainPanel(
            plotOutput("heat_s_plot")
        )))
# define server logic options
server2 <- function(input, output) {
    Input <-  reactive({
        switch(input$Group, unique(metadata$Group))
    })
    InputS <- reactive({
        order <- unique(c(which(metadata$Group %in% input$Group), which(metadata$Group == "Control")))
    })
    
    output$heat_s_plot <-renderPlot({
        plotFreqHeatmap(crispr_set3_all, min.freq = 0.1, order = InputS(), plot.text.size = 4, legend.text.size = 12, x.hjust = 1 )
    }, height = 750, width = 700)
}

# Run the application 
shinyApp(ui = ui2, server = server2)

```

## Variant frequency at the target site by sample

- The allele and mutation frequency was calculated using the CrispRVariants 
- The GCF_000003025.6_Sscrofa11.1_genomic database downloaded from the NCBI. 
- 1 control (CL) sample for the blastocysts and one control for the parthenotes 
- Displaying alleles with more the 10% frequency

## 

```{r freq_sample, out.height = "110%", out.width = "110%"}
vc <- variantCounts(crispr_set3_all)
ui2 <- fluidPage(
    # define title of the panel
    titlePanel("Samples"),
    # define input
    sidebarLayout( 
        sidebarPanel( 
            selectInput( inputId = "Samples", label = "Choose sample(s)", choices = metadata$Samples, multiple = TRUE)),
        # define output 
        mainPanel(
            plotOutput("heat_s_plot")
        )))
# define server logic options
server2 <- function(input, output) {
    Input <-  reactive({
        switch(input$Samples, metadata$Samples)
    })
    InputS <- reactive({
        order <- which(gsub('.sort.bam', '', colnames(vc)) %in% input$Samples)
    })
    
    output$heat_s_plot <- renderPlot({
        plotFreqHeatmap(crispr_set3_all, min.freq = 0.1, order = InputS(), plot.text.size = 4, legend.text.size = 12, x.hjust = 1 )
    }, height = 750, width = 700)
}

# Run the application 
shinyApp(ui = ui2, server = server2)

```

# Alignment plot 

- the alignemst plots were criated using the function PlotAlignrments and the annotation provided by the .GFF file from the loci downloaded at NCBI.



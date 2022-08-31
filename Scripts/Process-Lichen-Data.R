## Pre-Analysis Data processing for Lichen Surface Sterilization Data ##
# DADA2 has already been preformed #

# Input: 
# Outputs:



#### Load Data ####


# if (!requireNamespace("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("ampvis2", version = "3.14")

install.packages("remotes")
remotes::install_github("MadsAlbertsen/ampvis2")

library(phyloseq)
library(ggplot2)
library(rmarkdown)
library(dada2)
library(decontam)
library(phangorn)
library(DECIPHER)
library(vegan)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(reshape)
library(coin)
library(FSA)
library(ampvis2)
library(ape)
library(ShortRead)
library(tidyr)
library(gridExtra)


seqtab <- read.csv("Data/seqtab_final - seqtab_final(1).csv")
taxatab <- read.csv("Data/taxa - taxa.csv")
metadata <- read.csv("Data/lichen_metadata - metadata.csv")

colnames(seqtab)[1]
row.names(seqtab) <- seqtab[,1]
seqtab <- seqtab[,-1]

row.names(taxatab) <- taxatab[,1]
taxatab <- taxatab[,-1]

colnames(metadata)
row.names(metadata) <- metadata[,"full.2"]

seqtab <- as.matrix(seqtab)
taxatab <- as.matrix(taxatab)

### Create phyloseq object ####

# put what I assign the names to the data frames into the parentheses
ps <- phyloseq(otu_table(seqtab,taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxatab))

# extracts metadata: sample_data()
# extracts sequence counts: otu_table()
# extracts taxonomy table: taxa_table()

#### Decontaminate data ####

# example: ps_h20 <- subset_samples(ps_onlyfungi, SampleSubType!="Control_Kit")
# "ps_control" <- subset_samples("ps_data", SampleSubType!= )
?
## Using decontam to deal with our NC
#We do not have input DNA concentration, so we are using decontom's
#“prevalence” method. In this method, the prevalence 
#(presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence 
#in negative controls to identify contaminants.

#identify contaminants using the prevelence method
# Using the default threshold for a contaminant is 
#that it reaches a probability of 0.1 in the statistical test being performed.


#I tried to make a subset of the data with just the samples from the fukami lab??
ps_fukami <- subset_samples(ps, lab=="fukami")

#define the neg explicitly using Marina's code
#sample_data(ps)$is.neg <- sample_data(ps)$SampleType == "Kit"
sample_data(ps_fukami)$is.neg <- sample_data(ps_fukami)$treatment == "control"

# example: contamdf.prev.h20 <- isContaminant(ps_h20 , method="prevalence", neg="IsControl")
contamdf.prev <- isContaminant(ps_fukami , method="prevalence", neg="is.neg", threshold = 0.5)


#how many contaminants

# ex: table(contamdf.prev.h20$contaminant)
table(contamdf.prev$contaminant)

#which seq is contaminant

# ex: which(contamdf.prev.h20$contaminant) 
which(contamdf.prev$contaminant)

#Using a more strict threshold
#threshold=0.5, which will identify as contaminants 
#all sequences there are more prevalent in negative controls 
#than in positive samples


#Looking at the number of times these taxa were observed in 
#negative controls and positive samples:

# ex: ps.neg <- prune_samples(sample_data(ps_onlyfungi)$IsControl == "TRUE", ps_onlyfungi)
ps.neg <- prune_samples(sample_data(ps_fukami)$is.neg == "TRUE", ps_fukami)

#Make phyloseq object of presence-absence in true positive samples

# ex: ps.pos <- prune_samples(sample_data(ps_onlyfungi)$IsControl == "FALSE", ps_onlyfungi)
ps.pos <- prune_samples(sample_data(ps_fukami)$is.neg == "FALSE", ps_fukami)

#Make data.frame of prevalence in positive and negative samples 

#using prev threshold = 0.5

# ex: df.pres <- data.frame(prevalence.pos=taxa_sums(ps.pos.presence), prevalence.neg=taxa_sums(ps.neg.presence), contam.prev=contamdf.prev.h20.05$contaminant)
# ex: ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos, color=contam.prev)) + geom_point()
df.pres <- data.frame(prevalence.pos=taxa_sums(ps.pos), prevalence.neg=taxa_sums(ps.neg), contam.prev=contamdf.prev$contaminant)
ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos, color=contam.prev)) + geom_point()


#Get sequences that are contaminants

# ex: contaminants.h20 <- subset(contamdf.prev.h20.05, contaminant == "TRUE")
contaminants <- subset(contamdf.prev, contaminant == "TRUE")

#save these files as csv to investigate

# ex: write.csv(contaminants.h20, "contaminants.h20.05.csv")
write.csv(contaminants, "contaminants.csv")

#Combined in Excel while investigating and removed duplicates

# ex: contaminants <- read.csv("contaminants.combined.noduplicates.csv")
contaminants <- read.csv("contaminants.csv")

#We will want to remove these ASVs from further analyses

# ex: allTaxa = taxa_names(ps_onlyfungi)
# ex: newTaxa <- allTaxa[!(allTaxa %in% contaminants.list)]
# ex: ps_onlyfungi_NC = prune_taxa(newTaxa, ps_onlyfungi)\
allTaxa = taxa_names(ps_fukami)
newTaxa <- allTaxa[!(allTaxa %in% contaminants)]
ps_fukami_NC = prune_taxa(newTaxa, ps_fukami)

## goal: remove the contaminant taxa from ps (both labs) ##

allTaxa = taxa_names(ps)
newTaxa <- allTaxa[!(allTaxa %in% contaminants)]
ps_NC = prune_taxa(newTaxa, ps)

# ? how to see the names of the taxa that are contaminants ?

#### Rarefaction ####

## Using phyloseq loaded & controls removed ps object
#ex: ps.nocontrols

#defining function that is the wrapper -- Does not need to be changed
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
 #ex: x <- as(otu_table(physeq), "matrix")
 #ex: if (taxa_are_rows(physeq)) { x <- t(x) }


## This script is adapted from vegan `rarecurve` function
#ex: tot <- rowSums(x)
#ex: S <- rowSums(x > 0)
#ex: nr <- nrow(x)


# # #ex: 
rarefun <- function(i) {
  cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
  n <- seq(1, tot[i], by = step)
  if (n[length(n)] != tot[i]) {
    n <- c(n, tot[i])
  }
  y <- rarefy(x[i, ,drop = FALSE], n, se = se)
  if (nrow(y) != 1) {
    rownames(y) <- c(".S", ".se")
    return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
  } else {
    return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
  }
}
if (parallel) {
  out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
} else {
  out <- lapply(seq_len(nr), rarefun)
}
df <- do.call(rbind, out)
# # #

## Get sample data 
if (!is.null(sample_data(physeq, FALSE))) {
  sdf <- as(sample_data(physeq), "data.frame")
  sdf$Sample <- rownames(sdf)
  data <- merge(df, sdf, by = "Sample")
  labels <- data.frame(x = tot, y = S, Sample = rownames(x))
  labels <- merge(labels, sdf, by = "Sample")
}
# # #

## Add, any custom-supplied plot-mapped variables
if( length(color) > 1 ){
  data$color <- color
  names(data)[names(data)=="color"] <- deparse(substitute(color))
  color <- deparse(substitute(color))
}
if( length(label) > 1 ){
  labels$label <- label
  names(labels)[names(labels)=="label"] <- deparse(substitute(label))
  label <- deparse(substitute(label))
}
# # #

# # #
p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
p <- p + labs(x = "Sample Size", y = "Species Richness")
if (!is.null(label)) {
  p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                     size = 4, hjust = 0)
}
p <- p + geom_line()
if (se) { ## add standard error if available
  p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
}
if (plot) {
  plot(p)
}
invisible(p)


# # #

#plotting rarefaction curves, using step = 1000 reads
#ex: p = ggrare(ps.nocontrols, step = 1000, label = "SampleID", color = "SampleSubType") 

#drawing cut-off line at 20000 reads
#ex: p + geom_vline(xintercept = 20000, linetype="dashed")

#drawing cut-off line at 30000 reads
#ex: p + geom_vline(xintercept = 30000, linetype="dashed")

#drawing cut-off line at 10000 reads
#ex: p + geom_vline(xintercept = 10000, linetype="dashed")


#drawing cut-off line at 5000 reads
#ex: p + geom_vline(xintercept = 5000, linetype="dashed")


#ex: p
#Rarefaction:



## Let's look at library size to see distribution of reads across samples

## Look at library size
#ex: df <- as.data.frame(sample_data(ps.nocontrols)) # Put sample_data into a ggplot-friendly data.frame
#ex: df$LibrarySize <- sample_sums(ps.nocontrols)
#ex: df <- df[order(df$LibrarySize),]
#ex: df$Index <- seq(nrow(df))
#ex: ggplot(data=df, aes(x=Index, y=LibrarySize, color=Competion)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")
#ex: ggplot(data=df, aes(x=Index, y=LibrarySize, color=Treatment)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")



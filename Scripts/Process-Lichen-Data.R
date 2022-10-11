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
library(Rmisc)

set.seed(20123)

seqtab <- readRDS("Data/seqtab.RDS")
taxatab <- readRDS("Data/taxa.RDS")
metadata <- readRDS("Data/metadata.RDS")


### Create phyloseq object ####

# put what I assign the names to the data frames into the parentheses
ps <- phyloseq(otu_table(seqtab,taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(taxatab))

# extracts metadata: sample_data()
# extracts sequence counts: otu_table()
# extracts taxonomy table: taxa_table()

#### Remove Evernia ####

# remove specific genera from ps
taxatab_df <- as.data.frame(taxatab)
hosttaxa <- row.names(taxatab_df[which(taxatab_df$Genus == "g__Evernia"), ])
allTaxa = taxa_names(ps)
newTaxa <- allTaxa[!(allTaxa %in% hosttaxa)]
ps_final <- prune_taxa(newTaxa, ps)


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
x <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) { x <- t(x) }


## This script is adapted from vegan `rarecurve` function
tot <- rowSums(x)
S <- rowSums(x > 0)
nr <- nrow(x)


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
}


# # #

#plotting rarefaction curves, using step = 1000 reads
p = ggrare(ps_final, step = 1000, label = "location", color = "mean.tissue.nitrogen") 
# the label is a column name from the metadata of a variable I want to look at, can be blank

#drawing cut-off line at 20000 reads
#p + geom_vline(xintercept = 20000, linetype="dashed")

#drawing cut-off line at 30000 reads
#p + geom_vline(xintercept = 30000, linetype="dashed")

#drawing cut-off line at 10000 reads
#p + geom_vline(xintercept = 10000, linetype="dashed")


#drawing cut-off line at 5000 reads
#p + geom_vline(xintercept = 5000, linetype="dashed")



#Rarefaction:


## Maybe don't do anything with this
# Let's look at library size to see distribution of reads across samples

## Look at library size
#ex: df <- as.data.frame(sample_data(ps.nocontrols)) # Put sample_data into a ggplot-friendly data.frame
#ex: df$LibrarySize <- sample_sums(ps.nocontrols)
#ex: df <- df[order(df$LibrarySize),]
#ex: df$Index <- seq(nrow(df))
#ex: ggplot(data=df, aes(x=Index, y=LibrarySize, color=Competion)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")
#ex: ggplot(data=df, aes(x=Index, y=LibrarySize, color=Treatment)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")


## next steps: rarefy to different size
# rarefy even depth fucntion
# take output and run an ordination -- look at marinas big code -- ITS ordinations actual ordination and figures ordinations



#Even depth function:
#ex: ps.nocontrols.rare9434 <- rarefy_even_depth(ps.nocontrols, sample.size = 9434, replace = FALSE, rngseed = 5311) 
set.seed(20123)
ps_final_rare <- rarefy_even_depth(ps_final, sample.size = 2095, rngseed = FALSE, replace = FALSE) 
#5 samples and 165 OTUs were removed at samples size 2095

p = ggrare(ps_final_rare, step = 1000, label = "location", color = "Landscape.position") 

# df <- as.data.frame(sample_data(ps_final)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(ps_final)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Competion)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Treatment)) + geom_point() + geom_hline(yintercept = 10000, linetype="dashed")


### Statistical data ####

DistBC <- phyloseq::distance(ps_final_rare, method = "bray", type="samples")

## mean.tissue.nitrogen -- not significant
set.seed(20123)
adonis2(DistBC ~ mean.tissue.nitrogen, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Mean_tissue_nitrogen <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$mean.tissue.nitrogen)
set.seed(20123)
permutest(Mean_tissue_nitrogen, permutations = 9999) 
plot(Mean_tissue_nitrogen, label = F)

## mean.annual.precip -- significant USE!!!
set.seed(20123)
adonis2(DistBC ~ mean.annual.precip, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Mean_annual_precip <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$mean.annual.precip)
set.seed(20123)
permutest(Mean_annual_precip, permutations = 9999) 
plot(Mean_annual_precip, label = F)

## infected -- not significant
set.seed(20123)
adonis2(DistBC ~ infected, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Infected <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$infected)
set.seed(20123)
permutest(Infected, permutations = 9999) 
plot(Infected, label = F)

## location -- significant
set.seed(20123)
adonis2(DistBC ~ location, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Location <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$location)
set.seed(20123)
permutest(Location, permutations = 9999) 
plot(Location, label = F)

#Landscape.position -- significant USE!!!
set.seed(20123)
adonis2(DistBC ~ Landscape.position, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Landscape_position <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$Landscape.position)
set.seed(20123)
permutest(Landscape_position, permutations = 9999) 
plot(Landscape_position, label = F)

#winter.light.trans -- not significant
set.seed(20123)
adonis2(DistBC ~ winter.light.trans, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Winter_light_trans <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$winter.light.trans)
set.seed(20123)
permutest(Winter_light_trans, permutations = 9999) 
plot(Winter_light_trans, label = F)

#canopy.cover.measured -- not significant
set.seed(20123)
adonis2(DistBC ~ canopy.cover.measured, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Canopy_cover_measured <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$canopy.cover.measured)
set.seed(20123)
permutest(Canopy_cover_measured, permutations = 9999) 
plot(Canopy_cover_measured, label = F)

#proportion.infected -- not significant
set.seed(20123)
adonis2(DistBC ~ proportion.infected, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999) 

Proportion.infected <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$proportion.infected)
set.seed(20123)
permutest(Proportion.infected, permutations = 9999) 
plot(Proportion.infected, label = F)

#max.etr - photosynthesis capacity -- not significant
set.seed(20123)
adonis2(DistBC ~ max.etr, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Max.etr <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$max.etr)
set.seed(20123)
permutest(Max.etr, permutations = 9999) 
plot(Max.etr, label = F)

#mean.n15 -- significant!! maybe Use a little
set.seed(20123)
adonis2(DistBC ~ mean.n15, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Mean.n15 <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$mean.n15)
set.seed(20123)
permutest(Mean.n15, permutations = 9999) 
plot(Mean.n15, label = F)

#evernia.abundance -- not significant
set.seed(20123)
adonis2(DistBC ~ evernia.abundance, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Evernia.abundance <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$evernia.abundance)
set.seed(20123)
permutest(Evernia.abundance, permutations = 9999) 
plot(Evernia.abundance, label = F)

#Pb.e -- significant p=.04 maybe use a little
set.seed(20123)
adonis2(DistBC ~ Pb.e, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Lead <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$Pb.e)
set.seed(20123)
permutest(Lead, permutations = 9999) 
plot(Lead, label = F)

#Cd.e -- not significant p=.06
set.seed(20123)
adonis2(DistBC ~ Cd.e, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Cadmium <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$Cd.e)
set.seed(20123)
permutest(Cadmium, permutations = 9999) 
plot(Cadmium, label = F)

#aspect (direction it faces) -- significant USE!!!
set.seed(20123)
adonis2(DistBC ~ aspect, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Aspect <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$aspect)
set.seed(20123)
permutest(Aspect, permutations = 9999) 
plot(Aspect, label = F)

#num.asco (number of ascomata) -- significant
set.seed(20123)
adonis2(DistBC ~ num.asco, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 

Num.asco <- betadisper(DistBC, as(sample_data(ps_final_rare), "data.frame")$num.asco)
set.seed(20123)
permutest(Num.asco, permutations = 9999) 
plot(Num.asco, label = F)

#c.per.e -- significant
set.seed(20123)
adonis2(DistBC ~ c.per.e, as(sample_data(ps_final_rare), "data.frame"), permutations = 9999, na.action = na.omit) 



### ordination ####

ps_final_rare_ord <- ordinate(
  physeq = ps_final_rare, 
  method = "PCoA", 
  distance = "bray")

#ordination plot
plot_ordination(ps_final_rare, ps_final_rare_ord, color = "mean.annual.precip") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  labs(x="PCoA 22.1%" , y="PCoA 16.6%") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) 

#maybe use, maybe don't use, maybe use aspect instead
plot_ordination(ps_final_rare, ps_final_rare_ord, color = "Landscape.position", shape = "Landscape.position") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Landscape.position)) +
  labs(x = "PCoA1 (22.1%)", y = "PCoA2 (16.6%)") +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.margin = margin(c(.1,.1,.1,.1))
  ) 

#probable convergence failure


### estimate richness ####
richness <- estimate_richness(ps_final_rare, measures = "Observed")

richness <- cbind(richness, sample_data(ps_final_rare))


richness$location <- as.factor(richness$location)

set.seed(20123)
kruskal_test(Observed ~ location, distribution = approximate(nresample = 9999), data = richness)

dunnTest(Observed ~ location, data = richness, method = "bh") 

#effect size, average species count
summarySE(richness, groupvars = "location", measurevar = "Observed")



plot_richness(ps_final_rare, measures = "Observed", x = "location", color = "location") + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank()
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "richness")

### Genera distribution ####
ps_final_rare.RA <- transform_sample_counts(ps_final_rare, function(x) x / sum(x))

#define standard error function
se <- function(x) sqrt(var(x)/length(x))

AvgRA_o.g <- tax_glom(ps_final_rare.RA, taxrank = "Genus", NArm = FALSE)

#variation (originally .0002)
AvgRA99_O_variation <- filter_taxa(AvgRA_o.g, function(x) var(x) > .0005, TRUE)
#mean (originally .005)
#AvgRA99_O_mean <- filter_taxa(AvgRA_o.g, function(x) mean(x) > .008, TRUE)

df_o_v <- psmelt(AvgRA99_O_variation)
#df_o_m <- psmelt(AvgRA99_O_mean)

#group and calculate mean, sd and se for different taxonomic levels
avgs_o_v <- df_o_v %>% group_by(location, Family, Genus) %>%
  dplyr::summarise(mean = 100*mean(Abundance), sd = 100*sd(Abundance), se = 100*se(Abundance))

#avgs_o_m <- df_o_m %>% group_by(location, Family, Genus) %>%
#  dplyr::summarise(mean = 100*mean(Abundance), sd = 100*sd(Abundance), se = 100*se(Abundance))


ggplot(avgs_o_v, aes(x = Family, y = mean, fill = Genus)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
                width = .4, position = position_dodge(.9)) +
  facet_grid(~location) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 20)) + 
  ylab("Mean Relative Abundance")

#ggplot(avgs_o_m, aes(x = Family, y = mean, fill = Genus)) + 
 # geom_bar(stat = "identity", position = position_dodge()) +
  #geom_errorbar(aes(ymin = (mean - se), ymax = (mean + se)), 
   #             width = .4, position = position_dodge(.9)) +
  #facet_grid(~location) + 
  #theme_classic() +
  #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
   #     text = element_text(size = 20)) + 
  #ylab("Mean Relative Abundance")


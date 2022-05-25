# Phyloseq tutorial
library(tidyverse)
library(phyloseq)

# Some nice tutorials available on the website: https://joey711.github.io/phyloseq/
# Nice overview in phyloseq paper: https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0061217.g003

# Do the examples below with some real COPSAC 16S data.
# Does not work outside COPSAC, but can be substituted with any other phyloseq object for these examples.
load("/Volumes/Research/Microbiology/Shared/DADA2 Golden Russel Jul 2020/dada2_gtdb_2020_07_13.RData")

# 5 "slots" of connected data. It's super handy that phyloseq takes care of all the organization when subsetting, transforming, pruning etc.
phy

sample_names(phy)
taxa_names(phy)
nsamples(phy)
ntaxa(phy)

####
# otu_table is a special num. matrix that contains abundances/counts
otu_table(phy)[1:5, 1:5]
class(otu_table(phy))
otu_table(phy) %>% as.matrix %>% class  # does not hard convert
otu_table(phy) %>% as("matrix") %>% class  # does hard convert
otu_table(phy) %>% as.data.frame %>% class  # does hard convert, takes long time
otu_table(phy) %>% as("matrix") %>% as.data.frame %>% class  # does hard convert, goes faster!
# Get sample sums/totals
sample_sums(phy)
# Get totals for each taxon
taxa_sums(phy)

####
# tax_table is a special char. matrix that contains taxonomy information on features/ASVs. Can also put other feature metadata in there, but preferably to the right of the taxonomy as the hierarchy is required for tax_glom, see below
tax_table(phy)[1:10,]
tax_table(phy) %>% as.matrix %>% class  # does not hard convert
tax_table(phy) %>% as("matrix") %>% class  # does hard convert
# tax_table(phy) %>% as.data.frame %>% class  # does hard convert, takes long time
tax_table(phy) %>% as("matrix") %>% as.data.frame %>% class  # does hard convert, goes faster!
# Inspect the names of each rank
rank_names(phy)

####
# sample_data is a special data.frame which contains sample metadata
sample_data(phy) %>% head
sample_data(phy) %>% class
sample_data(phy) %>% as.data.frame %>% class  # does not hard convert
sample_data(phy) %>% data.frame %>% class # does hard convert
sample_data(phy) %>% as("data.frame") %>% class # does hard convert
# another handy function to inspect/extract sample metadata.
get_variable(phy) %>% head
get_variable(phy, c("SampleID", "ABCNO")) %>% head
# Inspect variable names
sample_variables(phy)

####
# phy_tree is a relational tree structure of taxa
phy_tree(phy)

###
# refseq are the reference sequences. Can be extracted and saved, blasted, compared to other studies, or whatever's needed
refseq(phy)
refseq(phy) %>% as.character %>% head
# 

######
#####
#### Utilities

# Subsetting
phy_small <- subset_samples(phy, Tray == "Tray01")
# Prune on tax abundances, good idea to always do after subsetting
phy_small <- phy_small %>% prune_taxa(taxa_sums(.) > 0, .)

# If you want eg. only 1 month Tracheal samples
phy_t1m <- subset_samples(phy, Sampletype == "Trachea" & Time == "1m")
# remember to prune taxa after!

# Subset on taxonomy
phy_small %>% tax_table %>% head
phy_small_fuso <- phy_small %>% subset_taxa(Genus == "Fusobacterium")
phy_small_fuso %>% tax_table

# prune samples
important_samples <- sample_names(phy) %>% sample(5)
phy_important <- prune_samples(important_samples, phy)

# Convert counts to relative abundances
phy_relabu <- phy_small %>% transform_sample_counts(function(x) x/sum(x))
sample_sums(phy_relabu)

# Extract counts and convert to z-scaled log relative abundances with pseudocount. Useful for downstream analysis, eg PLS, RF.
zrelabucounts <- phy_relabu %>% 
  otu_table %>% 
  as("matrix") %>% 
  apply(1, function(x) log(x + min(x[x != 0])) %>% scale) 
# now taxa are columns!
zrelabucounts[1:10,1:10]
rownames(zrelabucounts) <- sample_names(phy_relabu)

########
########
# Agglomerating

phy_small_family <- tax_glom(phy_small, "Family")
# Chooses most abundante OTU as archetype - only has consequences for taxa name, tree and refseq
phy_small
phy_small_family
taxa_names(phy_small_family)

phy_small_treecut <- tip_glom(phy_small, h = 0.1)
# Chooses archetype in the same way. So archetype taxonomy is carried forward.
phy_small_treecut


#######
#######  Assignment, adding variables, overwriting elements
#######

# Manually replace OTU data by extraction, conversion, replacement
phy_small_empty <- phy_small
empty_otus <- otu_table(phy_small)
empty_otus[,] <- 0
otu_table(phy_small_empty) <- empty_otus
sample_sums(phy_small_empty)

# Manually replace OTU data from scratch, this could be eg. DESeq2 normalized counts
cool_matrix <- matrix(0, nrow = ntaxa(phy_small), ncol = nsamples(phy_small))
dim(cool_matrix)
manual_otus <- otu_table(cool_matrix, taxa_are_rows = TRUE) # otu_table function is context dependent
otu_table(phy_small_empty) <- manual_otus  # taxa names do not match!
rownames(manual_otus) <- taxa_names(phy_small)
otu_table(phy_small_empty) <- manual_otus  # sample names do not match!
colnames(manual_otus) <- sample_names(phy_small)
otu_table(phy_small_empty) <- manual_otus  # Succes!

# Add variables
phy_small_added <- phy_small
# one at a time
sample_data(phy_small_added)$newvar <- rnorm(nsamples(phy_small_added))

# a whole chunk
stuff_to_add <- data.frame(ABCNO = get_variable(phy_small_added, "ABCNO"),
                           newvar2 = rnorm(nsamples(phy_small_added)),
                           newvar3 = rnorm(nsamples(phy_small_added)))
sample_data(phy_small_added) <- phy_small_added %>% 
  sample_data %>% 
  as("data.frame") %>% 
  rownames_to_column %>% 
  left_join(stuff_to_add, by = "ABCNO") %>% 
  column_to_rownames %>% 
  sample_data


#######
#######  Alpha diversity
#######

alpha <- estimate_richness(phy_small)
alpha_sub <- estimate_richness(phy_small, measures = c("Observed", "Shannon"))
head(alpha)
head(alpha_sub)
# this can be put back in sample_data

#######
#######  Beta diversity
#######
library(vegan)

phy_small
phy_small_bray_d <- distance(phy_small, "bray")
phy_small_bray_o <- ordinate(phy_small, "PCoA", phy_small_bray_d)
plot_ordination(phy_small, phy_small_bray_o, color = "Lot.number") + 
  stat_ellipse(level = 0.68) +
  scale_color_brewer(palette = "Set1")
adonis(phy_small_bray_d ~ Lot.number, data = get_variable(phy_small))

phy_small_unifrac_d <- distance(phy_small, "unifrac")
phy_small_unifrac_o <- ordinate(phy_small, "PCoA", phy_small_unifrac_d)
plot_ordination(phy_small, phy_small_unifrac_o, color = "Lot.number") + 
  stat_ellipse(level = 0.68) +
  scale_color_brewer(palette = "Set1")
adonis(phy_small_unifrac_d ~ Lot.number, data = get_variable(phy_small))

###
### MRAdat - a useful taxa summary dataframe
###

make_mradat <- function(phy){
  data.frame(tax = taxa_names(phy),
             mra = rowMeans(phy %>% transform_sample_counts(function(x) x/sum(x)) %>% 
                              otu_table),
             prevalence = rowMeans(otu_table(phy) > 0),
             tax_table(phy) %>% as("matrix") %>% as.data.frame)
}


mradat <- phy_small %>% make_mradat
head(mradat)

# Choose taxa above the line in the plot below
chosen_taxa <- mradat %>% 
  filter(prevalence >= 0.25*(1e-4/mra)^.5) %>% 
  pull(tax) %>% 
  as.character

## MRA vs. prevalence of entire dataset, useful for making cutoffs for choosing common taxa for DA testing or machine learning
## Needs fullscreen
ggplot(mradat, aes(mra, prevalence)) +
  geom_point(aes(color = Phylum)) +
  ggrepel::geom_label_repel(aes(label = Species), data = . %>% filter(tax %in% chosen_taxa), size = 3, label.padding = .1) +
  scale_x_log10() +
  stat_function(fun = function(x) 0.25*(1e-4/x)^.5, linetype = "dashed", inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 1))

# Can also do mradats on agglomerated versions
mradat_f <- phy_small_family %>% make_mradat
mradat_f %>% arrange(desc(prevalence)) %>% head

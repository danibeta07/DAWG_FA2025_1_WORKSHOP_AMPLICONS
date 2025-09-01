amplicons_workshop
================

DADA2 Analysis Pipeline for 16S and ITS Amplicon Data project
PRJNA880162
URL:<https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA880162&o=acc_s%3Aa>
Bacterial 16S and fungal ITS amplicons were generated from DNA isolated
from human breastmilk and human infant fecal samples. Fecal samples were
collected from healthy infants at 1 and 6 months of age. Milk samples
were collected from infant mothers at 1 month postnatal clinic visits.
The dataset is unique in that the majority of infants were not exposed
to antibiotics and were exclusively breastfed, allowing better control
of two major potential gut microbiome and mycobiome modulators.

## Load libraries

``` r
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(vegan)
library(knitr)
setwd("~/Downloads/1_workshop_amplicons")
#setwd("/storage/group/efc5279/default/DaniB/DAWG/FA2025/1_workshop_amplicons")
```

# 16S 
## Define file paths

``` r
base_path <- "~/Downloads/1_workshop_amplicons"
#base_path <- "/storage/group/efc5279/default/DaniB/DAWG/FA2025/1_workshop_amplicons"
path_16S <- file.path(base_path, "16S_FASTQ")
filt_path_16S <- file.path(base_path, "16S_filtered")
taxa_db_16S <- file.path(base_path, "silva_nr99_v138.1_train_set.fa")
```

## 16S dataset

``` r
fnFs_16S <- sort(list.files(path_16S, pattern = "_1.fastq", full.names = TRUE))
fnRs_16S <- sort(list.files(path_16S, pattern = "_2.fastq", full.names = TRUE))
#extract sample names
sample.names_16S <- sapply(strsplit(basename(fnFs_16S), "_"), `[`, 1)
```

## Inspect reads quality

``` r
plotQualityProfile(fnFs_16S[10:11])
plotQualityProfile(fnRs_16S[10:11])
```

## Filter and trimming

``` r
filtFs_16S <- file.path(filt_path_16S, paste0(sample.names_16S, "_F_filt.fastq.gz"))
filtRs_16S <- file.path(filt_path_16S, paste0(sample.names_16S, "_R_filt.fastq.gz"))
names(filtFs_16S) <- sample.names_16S
names(filtRs_16S) <- sample.names_16S

# Set filtering parameters based on V4 region and quality plots
out_16S <- filterAndTrim(fnFs_16S, filtFs_16S, fnRs_16S, filtRs_16S,
                       truncLen = c(240, 240),
                       maxN = 0, #discards reads with ambigous bases
                       maxEE = c(2, 2), #maximum number of expected errors
                       truncQ = 2, #truncates read with a quality score less or equals to 2
                       rm.phix = TRUE, #removes reads that match phix genome
                       compress = TRUE,
                       multithread = TRUE)
print("Reads surviving 16S filtering:")
print(out_16S)
```

## Core DADA2 pipeline

``` r
#Build specific error model from your dataset. Here we figure out error rates by looking at the composition and quality scores of the reads
errF_16S <- learnErrors(filtFs_16S, multithread = TRUE)
errR_16S <- learnErrors(filtRs_16S, multithread = TRUE)

#Find all the identical read sequences and collapses them into unique sequences
derepFs_16S <- derepFastq(filtFs_16S, verbose = TRUE)
derepRs_16S <- derepFastq(filtRs_16S, verbose = TRUE)
names(derepFs_16S) <- sample.names_16S
names(derepRs_16S) <- sample.names_16S

#Take the unique sequences from dereplication and applies the custom error model to decide whether a rare unique sequence is a true biological variant or just a sequencing error of a more abundant sequence. It then corrects these errors to infer the original, error-free biological sequences (ASVs).
dadaFs_16S <- dada(derepFs_16S, err = errF_16S, multithread = TRUE)
dadaRs_16S <- dada(derepRs_16S, err = errR_16S, multithread = TRUE)

#Take the denoised forward and reverse reads and merges them into a single, full-length sequence. It aligns them using their overlapping region and ensures the overlap is nearly identical
mergers_16S <- mergePairs(dadaFs_16S, derepFs_16S, dadaRs_16S, derepRs_16S, verbose = TRUE)

#Organize all the merged sequences into a feature table
seqtab_16S <- makeSequenceTable(mergers_16S)

#Remove chimeric ASVs from the table
seqtab.nochim_16S <- removeBimeraDenovo(seqtab_16S, method = "consensus", multithread = TRUE, verbose = TRUE)
```

## Taxononmy classification

``` r
taxa_16S <- assignTaxonomy(seqtab.nochim_16S, taxa_db_16S, multithread = TRUE)
```

## Phyloseq object

``` r
metadata_path <- file.path(base_path, "filtered_metadata.csv")
samdf <- read.csv(metadata_path)
rownames(samdf) <- samdf$Run
samdf_16S <- samdf[sample.names_16S, ]

ps_16S <- phyloseq(otu_table(seqtab.nochim_16S, taxa_are_rows = FALSE),
               sample_data(samdf_16S),
               tax_table(taxa_16S))

saveRDS(ps_16S, file = file.path(base_path, "16S_phyloseq_object.rds"))

print(ps_16S)
```

# ITS 
## Define file paths

``` r
path_ITS <- file.path(base_path, "ITS_FASTQ")
filt_path_ITS <- file.path(base_path, "ITS_filtered")
taxa_db_ITS <- file.path(base_path, "sh_general_release_dynamic_s_19.02.2025_dev.fasta")
```

## ITS dataset

``` r
fnFs_ITS <- sort(list.files(path_ITS, pattern = "_1.fastq", full.names = TRUE))
fnRs_ITS <- sort(list.files(path_ITS, pattern = "_2.fastq", full.names = TRUE))
sample.names_ITS <- sapply(strsplit(basename(fnFs_ITS), "_"), `[`, 1)
```

## Inspect reads quality

``` r
plotQualityProfile(fnFs_ITS[10:11])
plotQualityProfile(fnRs_ITS[10:11])
```

## Filter and trimming

``` r
filtFs_ITS <- file.path(filt_path_ITS, paste0(sample.names_ITS, "_F_filt.fastq.gz"))
filtRs_ITS <- file.path(filt_path_ITS, paste0(sample.names_ITS, "_R_filt.fastq.gz"))
names(filtFs_ITS) <- sample.names_ITS
names(filtRs_ITS) <- sample.names_ITS

out_ITS <- filterAndTrim(fnFs_ITS, filtFs_ITS, fnRs_ITS, filtRs_ITS,
                       trimLeft = c(22, 20),
                       maxN = 0,
                       maxEE = c(2, 5),
                       truncQ = 2,
                       rm.phix = TRUE,
                       compress = TRUE,
                       multithread = TRUE)
print("Reads surviving ITS filtering:")
print(out_ITS)
```

## Core DADA2 pipeline

``` r
errF_ITS <- learnErrors(filtFs_ITS, multithread = TRUE)
errR_ITS <- learnErrors(filtRs_ITS, multithread = TRUE)

derepFs_ITS <- derepFastq(filtFs_ITS, verbose = TRUE)
derepRs_ITS <- derepFastq(filtRs_ITS, verbose = TRUE)
names(derepFs_ITS) <- sample.names_ITS
names(derepRs_ITS) <- sample.names_ITS

dadaFs_ITS <- dada(derepFs_ITS, err = errF_ITS, multithread = TRUE)
dadaRs_ITS <- dada(derepRs_ITS, err = errR_ITS, multithread = TRUE)

mergers_ITS <- mergePairs(dadaFs_ITS, derepFs_ITS, dadaRs_ITS, derepRs_ITS, verbose = TRUE)
seqtab_ITS <- makeSequenceTable(mergers_ITS)
seqtab.nochim_ITS <- removeBimeraDenovo(seqtab_ITS, method = "consensus", multithread = TRUE, verbose = TRUE)
```

## Taxononmy classification

``` r
taxa_ITS <- assignTaxonomy(seqtab.nochim_ITS, taxa_db_ITS, multithread = TRUE)
```

## Phyloseq object

``` r
metadata_path <- file.path(base_path, "filtered_metadata.csv")
samdf <- read.csv(metadata_path)
rownames(samdf) <- samdf$Run

samdf_ITS <- samdf[sample.names_ITS, ]

ps_ITS <- phyloseq(otu_table(seqtab.nochim_ITS, taxa_are_rows = FALSE),
               sample_data(samdf_ITS),
               tax_table(taxa_ITS))

saveRDS(ps_ITS, file = file.path(base_path, "ITS_phyloseq_object.rds"))
print(ps_ITS)
```

## Data Analysis 
### Alpha diversity 16S

``` r
#Load phyloseq objects
ps_16S <- readRDS(file.path(base_path, "16S_phyloseq_object.rds"))
ps_ITS <- readRDS(file.path(base_path, "ITS_phyloseq_object.rds"))

# Alpha diversity is the diversity within a single sample.
# We will use plot_richness to create boxplots of different diversity metrics.
p_alpha_16S <- plot_richness(ps_16S, x = "Disease", measures = c("Observed", "Shannon")) +
  geom_boxplot(aes(fill = Disease), alpha = 0.8) +
  labs(title = "16S Alpha Diversity", x = "Sample Group") +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(size = 12))
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
print(p_alpha_16S)
```

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->
### Beta diversity 16S

``` r
# Beta diversity measures the differences between samples. We will use PCoA
# ordination on a Bray-Curtis dissimilarity matrix.

# Remove samples with zero reads. Normalize counts and calculate Bray-Curtis distance
ps_16S_pruned <- prune_samples(sample_sums(ps_16S) > 0, ps_16S)
ps_rel_16S <- transform_sample_counts(ps_16S_pruned, function(x) x / sum(x))
bray_dist_16S <- phyloseq::distance(ps_rel_16S, method = "bray")

# Perform PCoA
pcoa_16S <- ordinate(ps_rel_16S, method = "PCoA", distance = bray_dist_16S)

# Create the plot
# Replace "Your_Metadata_Variable" with the same column name as before.
p_pcoa_16S <- plot_ordination(ps_rel_16S, pcoa_16S, color = "Disease") +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "16S PCoA (Bray-Curtis)", color = "Sample Group") +
  stat_ellipse(aes(group = Disease), level = 0.95) +
  theme_classic()

print(p_pcoa_16S)
```

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->
### Relative abundance 16S

``` r
# Transform counts to relative abundance and agglomerate at the Phylum level
ps_phylum_16S <- tax_glom(ps_16S, taxrank = "Phylum")
ps_phylum_rel_16S <- transform_sample_counts(ps_phylum_16S, function(x) x / sum(x))

# Create the plot
p_bar_16S <- plot_bar(ps_phylum_rel_16S, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  labs(title = "16S Relative Abundance by Phylum", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

    ## Warning in psmelt(physeq): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

``` r
print(p_bar_16S)
```

    ## Warning: Removed 34 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

    ## Warning: Removed 34 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Alpha diversity ITS

``` r
p_alpha_ITS <- plot_richness(ps_ITS, x = "Disease", measures = c("Observed", "Shannon")) +
  geom_boxplot(aes(fill = Disease), alpha = 0.8) +
  labs(title = "ITS Alpha Diversity", x = "Sample Group") +
  theme_classic() +
  theme(legend.position = "none", strip.text = element_text(size = 12))

print(p_alpha_ITS)
```

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Beta diversity ITS

``` r
# Normalize counts and calculate Bray-Curtis distance
ps_ITS_pruned <- prune_samples(sample_sums(ps_ITS) > 0, ps_ITS)
ps_rel_ITS <- transform_sample_counts(ps_ITS_pruned, function(x) x / sum(x))
bray_dist_ITS <- phyloseq::distance(ps_rel_ITS, method = "bray")

# Perform PCoA
pcoa_ITS <- ordinate(ps_rel_ITS, method = "PCoA", distance = bray_dist_ITS)

# Create the plot
# Replace "Your_Metadata_Variable" with the same column name as before.
p_pcoa_ITS <- plot_ordination(ps_rel_ITS, pcoa_ITS, color = "Disease") +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "ITS PCoA (Bray-Curtis)", color = "Sample Group") +
  stat_ellipse(aes(group = Disease), level = 0.95) +
  theme_classic()

print(p_pcoa_ITS)
```

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Relative abundance ITS

``` r
ps_phylum_ITS <- tax_glom(ps_ITS_pruned, taxrank = "Phylum")
ps_phylum_rel_ITS <- transform_sample_counts(ps_ITS_pruned, function(x) x / sum(x))

# Create the plot
p_bar_ITS <- plot_bar(ps_phylum_rel_ITS, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  labs(title = "ITS Relative Abundance by Phylum", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

    ## Warning in psmelt(physeq): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

``` r
print(p_bar_ITS)
```

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

### FUNGuild

``` r
otu_table_its <- as.data.frame(otu_table(ps_ITS))

# Transpose the OTU table so taxa become rows and samples become columns
# This is the crucial step to match the taxonomy table's structure.
otu_table_transposed <- as.data.frame(t(otu_table_its))

#Extract and format the taxonomy table
tax_table_its <- as.data.frame(tax_table(ps_ITS))
tax_table_its$taxonomy <- paste(
  tax_table_its$Kingdom, tax_table_its$Phylum, tax_table_its$Class,
  tax_table_its$Order, tax_table_its$Family, tax_table_its$Genus,
  tax_table_its$Species, sep = ";")

# Combine the transposed OTU table with the taxonomy column.
otu_table_with_taxonomy <- cbind(otu_table_transposed, taxonomy = tax_table_its$taxonomy)

#Add the '#OTUID' column 
otu_table_with_taxonomy <- cbind("#OTUID" = rownames(otu_table_with_taxonomy), otu_table_with_taxonomy)

#Define the output file path.
corrected_output_path <- file.path(base_path, "FUNGuild_corrected_input.txt")

# Save the table
write.table(otu_table_with_taxonomy, file = corrected_output_path, sep = "\t", row.names = FALSE, quote = FALSE)
```

### Run FUNGuild

``` python
cd ~/Downloads/1_workshop_amplicons/FUNGuild
python Guilds_v1.1.py -otu ~/Downloads/1_workshop_amplicons/FUNGuild_corrected_input.txt -db fungi

#Result saved to '/Users/daniela/Downloads/1_workshop_amplicons/FUNGuild_input_taxa.guilds.txt'
```

### Plot FUNGuild results

``` r
#Define paths
funguild_input_file <- file.path(base_path, "FUNGuild_corrected_input.txt")
funguild_output_file <- file.path(base_path, "FUNGuild_corrected_input.guilds.txt")

funguild_input <- read.delim(funguild_input_file, sep = "\t", header = TRUE)
funguild_results <- read.delim(funguild_output_file, sep = "\t", header = TRUE)

funguild_results$Taxon <- funguild_input$"X.OTUID" 

ps_long <- psmelt(ps_ITS)
```

    ## Warning in psmelt(ps_ITS): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

``` r
ps_long_guilds <- dplyr::left_join(ps_long, funguild_results, by = c("OTU" = "Taxon"))

#Filter out unassigned guilds
ps_guilds_filtered <- ps_long_guilds[!is.na(ps_long_guilds$Guild) & ps_long_guilds$Guild != "", ]

# Calculate relative abundance
guild_rel_abundance <- ps_guilds_filtered %>%
  group_by(Sample, Disease, Guild) %>%
  summarise(GuildAbundance = sum(Abundance), .groups = 'drop') %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = GuildAbundance / sum(GuildAbundance)) %>%
  ungroup()

#Create the plot
p_guilds <- ggplot(guild_rel_abundance, aes(x = Sample, y = Relative_Abundance, fill = Guild)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Disease, scales = "free_x", nrow = 1) +
  labs(
    title = "Fungal Guild Composition by Disease State",
    x = "Samples",
    y = "Relative Abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"), # Adjusts the size of the colored boxes
    legend.text = element_text(size = 9),     # Adjusts the font size of legend text
    legend.title = element_text(size = 9)     # Adjusts the font size of the legend title
  )

print(p_guilds)
```

    ## Warning: Removed 100 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
ggsave("Guilds.pdf", p_guilds, height = 10, width = 30)
```

    ## Warning: Removed 100 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

### Human associated

``` r
# Create a list of genera known to be associated with humans
# This list is based on common human pathogens and commensals from literature and CDC lists.
human_fungal_genera <- c(
  "Candida", "Aspergillus", "Cryptococcus", "Pneumocystis", "Mucor",
  "Rhizopus", "Fusarium", "Scedosporium", "Malassezia", "Trichophyton",
  "Histoplasma", "Coccidioides", "Blastomyces", "Penicillium", "Saccharomyces"
)

# Extract the taxonomy table from your phyloseq object
tax_table_its <- as.data.frame(tax_table(ps_ITS))

# Check which of your ASVs belong to these human-associated genera
matches <- grepl(paste(human_fungal_genera, collapse = "|"), tax_table_its$Genus, ignore.case = TRUE)

# Create a new table showing only the human-associated fungi found in your data
human_associated_taxa <- tax_table_its[matches, ]

# Get the list of ASV IDs for the human-associated taxa we found
human_asv_ids <- rownames(human_associated_taxa)

# Prune the phyloseq object to keep only these human-associated ASVs
ps_its_human <- prune_taxa(human_asv_ids, ps_ITS)

# Transform the counts to relative abundance for accurate comparison
ps_its_human_rel <- transform_sample_counts(ps_its_human, function(OTU) OTU / sum(OTU))

# Melt the data into a long format suitable for ggplot2
plot_data <- psmelt(ps_its_human_rel)
```

    ## Warning in psmelt(ps_its_human_rel): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

``` r
# Create the bar plot
p_human_fungi <- ggplot(plot_data, aes(x = Genus, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "summary", fun = "mean") + # Bar height is the mean abundance across all samples
  labs(
    title = "Abundance of Human-Associated Fungal Genera",
    x = "Fungal Genus",
    y = "Mean Relative Abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right"
  )

print(p_human_fungi)
```

    ## Warning: Removed 228 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
ggsave("HumanFungi.pdf", p_human_fungi, height = 10, width = 30)
```

    ## Warning: Removed 228 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

``` r
## By Disease
p_human_fungi_by_disease <- ggplot(plot_data, aes(x = Genus, y = Abundance, fill = Genus)) +
  geom_bar(stat = "summary", fun = "mean") +
  facet_wrap(~Disease, scales = "free_y") + # Creates separate panels for Y and N
  labs(
    title = "Abundance of Human-Associated Fungi by Disease Status",
    x = "Fungal Genus",
    y = "Mean Relative Abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none", # Hiding legend as colors match x-axis labels
    strip.text = element_text(size = 12, face = "bold")
  )

# Print the plot
print(p_human_fungi_by_disease)
```

    ## Warning: Removed 228 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

![](amplicons_workshop_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
ggsave("HumanFungiByDisease.pdf", p_human_fungi_by_disease, height = 10, width = 30)
```

    ## Warning: Removed 228 rows containing non-finite outside the scale range
    ## (`stat_summary()`).

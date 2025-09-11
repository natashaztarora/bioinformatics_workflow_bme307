# R studio

### 1. Package Installation and Library Loading

#### Install BiocManager {: data-toc-label="" }

```r
install.packages("BiocManager")
```

#### Install Bioconductor packages {: data-toc-label="" }

```r
BiocManager::install(c("phyloseq", "Biostrings", "S4Vectors", "IRanges", "XVector"))
```

#### Install CRAN packages {: data-toc-label="" }

```r
install.packages(c("vegan", "ape", "data.table", "Rcpp", "forcats", "tidyverse"))
```

Quick check:

```r
library(phyloseq)
packageVersion("phyloseq")
```

#### Install **qiime2R**  {: data-toc-label="" }

```r
install.packages("remotes")
remotes::install_url(
  "https://github.com/jbisanz/qiime2R/archive/refs/heads/master.zip",
  dependencies = TRUE
)
```

#### Install **phyloseq-extended**  {: data-toc-label="" }

```r
install.packages("remotes")   # only needed once, can skip if already installed
library(remotes)
remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
```

#### Alternative installation (only if `install_github()` fails with HTTP error 401)  {: data-toc-label="" }

```r
# install.packages("gitcreds")   # if not already installed
# library(gitcreds)
# gitcreds::gitcreds_get()       # check credentials
# gitcreds::gitcreds_delete()    # delete stored GitHub token if invalid
# remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
```

Quick check:

```r
library(phyloseq.extended)
packageVersion("phyloseq.extended")
```

#### Load libraries {: data-toc-label="" }

```r
library(qiime2R)
library(tidyverse)
library(vegan)
library(ape)
library(Rcpp)
library(data.table)
library(phyloseq.extended)
```

#### Version check {: data-toc-label="" }

```r
pkgs <- c("phyloseq", "qiime2R", "microbiome", "vegan", "ape", "tidyverse")
sapply(pkgs, \(p) as.character(packageVersion(p)))
```

---

### 2. Import QIIME2 files to phyloseq 

```r
# Define path to QIIME2 files
files_path <- "/Path/to/Practical_materials_uploads/QIIME2_files"

# (Optional) Read tree file to check contents
read_qza(file.path(files_path, "rooted-tree.qza"))

# Convert QIIME2 artefacts to a phyloseq object
pseq <- qza_to_phyloseq(
  features  = file.path(files_path, "table.qza"),
  taxonomy  = file.path(files_path, "taxonomy.qza"),
  metadata  = file.path(files_path, "metadata.tsv"),
  tree      = file.path(files_path, "rooted-tree.qza")
)

# Display phyloseq object
pseq
```

???+ info "Expected Output"
    ```
    phyloseq-class experiment-level object 
    otu_table()   OTU Table:         [ 502 taxa and 18 samples ] 
    sample_data() Sample Data:       [ 18 samples by 2 sample variables ]
    tax_table()   Taxonomy Table:    [ 502 taxa by 7 taxonomic ranks ]
    phy_tree()    Phylogenetic Tree: [ 502 tips and 500 internal nodes ]
    ```

```r
# Save and reload phyloseq object
save(pseq, file = "pseq.RData")
load("pseq.RData")
```

---

### 3. Explore the `phyloseq` Object

In this section we explore different `phyloseq` accessors to get familiar with the contents of our object.

#### Inspect the three core tables {: data-toc-label="" }

```r
head(otu_table(pseq))     # abundance matrix
head(tax_table(pseq))     # taxonomy
head(sample_data(pseq))   # metadata
```

#### Basic information  {: data-toc-label="" }

```r
ntaxa(pseq)               # number of taxa
nsamples(pseq)            # number of samples
sample_names(pseq)        # sample IDs
head(taxa_names(pseq))    # taxon (ASV) IDs
```

#### Relabel the ASVs with numbers {: data-toc-label="" }

```r
original_ids <- taxa_names(pseq)
asv_labels <- paste0("ASV", seq_len(ntaxa(pseq)))
taxa_names(pseq) <- asv_labels
tax_table(pseq) <- cbind(tax_table(pseq), OriginalID = original_ids)

head(taxa_names(pseq))
```

#### More accessors {: data-toc-label="" }

```r
sample_sums(pseq)                     # total reads per sample
head(taxa_sums(pseq))                 # total reads per taxon
sample_sums(pseq)[sample_names(pseq)=="IL10-2"]  # depth of one sample
rank_names(pseq)                      # available taxonomic ranks
sample_variables(pseq)                # metadata variables
get_variable(pseq, "type")            # values of a metadata variable
```

#### Extract sample-specific information  {: data-toc-label="" }

```r
# Counts for one sample (if taxa are rows, samples are columns)
otu_table(pseq)[, "WT4"]

# Metadata for one sample
sample_data(pseq)["WT4", ]

# Prune the object to only WT4
prune_samples("WT4", pseq)
```

#### Taxonomy queries {: data-toc-label="" }

```r
get_taxa_unique(pseq, "Phylum")
get_taxa_unique(pseq, "Family")
get_taxa_unique(pseq, "Class")
get_taxa_unique(pseq, "Genus")
```

---


### 4. Read depth

#### Show distribution of read counts {: data-toc-label="" }
```r
# total reads per sample
ss <- sample_sums(pseq)

hist(
  ss,
  breaks = 30,
  xaxt = "n",
  main = "Distribution of total reads per sample",
  xlab = "Reads per sample"
)
axis(1, at = pretty(ss))
```

```r
# Convert sample_data to a plain data.frame
meta <- as.data.frame(sample_data(pseq))

# Assuming your metadata has a column named 'type'
pseq_sum <- data.frame(
  sampleID    = sample_names(pseq),
  read_depth  = ss,
  type        = meta$type,
  row.names   = sample_names(pseq),
  check.names = FALSE
)

# Peek at the first few rows
head(pseq_sum)
```


```r
# Order samples by depth for a clearer plot
pseq_sum$sampleID <- factor(
  pseq_sum$sampleID,
  levels = pseq_sum$sampleID[order(pseq_sum$read_depth, decreasing = TRUE)]
)

read_depth_sum <- ggplot(pseq_sum, aes(x = sampleID, y = read_depth, fill = type)) +
  geom_col() +
  labs(
    title = "Total number of reads",
    x = "Samples (ordered by depth)",
    y = "Number of reads"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

read_depth_sum
```

### 5. Rarefaction curve and alpha diversity

#### Rarefaction curves with `vegan`

```r
library(phyloseq)
library(vegan)

# 1) Get counts as a plain matrix with SAMPLES IN ROWS (what vegan expects)
mat <- as(otu_table(pseq), "matrix")
if (taxa_are_rows(pseq)) mat <- t(mat)

# Optional: drop all-zero samples (rare but safer)
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
```

```r
# 2) Quick rarefaction curves
hist(rowSums(mat), main = "Library sizes", xlab = "Reads per sample")

raremax <- min(rowSums(mat))

rarecurve(
  mat,
  step = 1000,
  label = TRUE,
  cex = 0.5
)
```

```r
# 3) use phyloseq.extended to plot the rarefaction curves
library(phyloseq.extended)

phyloseq.extended::ggrare(
  pseq,
  step = 1000,
  color = "type",
  se = TRUE
)
```

```r
# install once (uncomment if needed)
# remotes::install_github("jbisanz/phyloseq-extended")

library(phyloseq.extended)

ggrare(
  pseq,
  step = 1000,
  color = "type",
  se = TRUE
)
```


```r
# 4) Plot Observed Richness

ObsR_plot <- plot_richness(
  pseq,
  x = "type",
  color = "type",
  measures = c("Observed", "Shannon")
) + 
  geom_boxplot()

ObsR_plot   # just print it
```

### 6. Normalize the data

#### Normalize read counts to **relative abundance** so samples are comparable.

```r
library(phyloseq)

# Normalize to proportions within each sample
pseq_normal <- transform_sample_counts(pseq_nonzero, function(x) x / sum(x))

# Inspect the normalized object
head(otu_table(pseq_normal))

# Sanity check: each sample should sum to ~1
summary(sample_sums(pseq_normal))
```

### 7. Taxonomic composition

#### Check unassigned/missing taxonomy

```r
library(phyloseq)
library(dplyr)

stopifnot(!is.null(tax_table(pseq_normal)))  # ensure taxonomy is present

check_taxonomy_quality <- function(ps) {
  tt <- as(tax_table(ps), "matrix")
  # ensure character
  tt <- apply(tt, 2, function(x) as.character(x))
  # label vector to detect common “unknowns”
  is_unassigned <- function(x) {
    y <- tolower(trimws(x))
    y %in% c("unassigned", "unknown", "uncultured", "incertae sedis", "incertae_sedis")
  }
  out <- apply(tt, 2, function(col) {
    n_total <- length(col)
    n_na    <- sum(is.na(col))
    n_blank <- sum(trimws(col) == "", na.rm = TRUE)
    n_unas  <- sum(is_unassigned(col), na.rm = TRUE)
    c(Total = n_total, Missing = n_na, Blank = n_blank, Unassigned = n_unas)
  })
  as.data.frame(t(out))
}

check_taxonomy_quality(pseq_normal)
```

```r
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
display.brewer.all()

# Parameters
rank  <- "Genus"   # e.g., "Phylum", "Class", "Order", "Family", "Genus"
topN  <- 15        # keep top-N taxa globally; rest collapsed into "Other"

# 1) Agglomerate to chosen rank
stopifnot(rank %in% rank_names(pseq_normal))
ps_glom <- tax_glom(pseq_normal, taxrank = rank, NArm = TRUE)

# 2) Melt to tidy long format
df <- psmelt(ps_glom) %>%
  mutate(
    Taxon = dplyr::coalesce(!!sym(rank), "Unassigned")
  )

# 3) Identify global top-N taxa by mean relative abundance
top_taxa <- df %>%
  group_by(Taxon) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = topN) %>%
  pull(Taxon)

# 4) Collapse others to "Other" and re-sum within sample
df_bar <- df %>%
  mutate(Taxon_collapsed = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
  group_by(Sample, Taxon_collapsed) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  left_join(
    as(sample_data(ps_glom), "data.frame") %>% tibble::rownames_to_column("Sample"),
    by = "Sample"
  )

# 5) Plot
TopTaxa <- ggplot(df_bar, aes(x = Sample, y = Abundance, fill = Taxon_collapsed)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Sample",
    y = "Relative abundance",
    fill = rank,
    title = paste0("Taxonomic composition (", rank, ", top ", topN, ")")
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

TopTaxa  + scale_fill_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(16))

```
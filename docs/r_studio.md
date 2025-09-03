# R studio

### 1. Package Installation and Library Loading

#### Install BiocManager
```r
install.packages("BiocManager")
````

#### Install Bioconductor packages

```r
BiocManager::install(c("phyloseq", "Biostrings", "S4Vectors", "IRanges", "XVector"))
```

#### Install CRAN packages

```r
install.packages(c("vegan", "ape", "data.table", "Rcpp", "forcats", "tidyverse"))
```

Quick check:

```r
library(phyloseq)
packageVersion("phyloseq")
```

#### Install **qiime2R**

```r
install.packages("remotes")
remotes::install_url(
  "https://github.com/jbisanz/qiime2R/archive/refs/heads/master.zip",
  dependencies = TRUE
)
```

#### Install **phyloseq-extended**

```r
install.packages("remotes")   # only needed once, can skip if already installed
library(remotes)
remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
```

#### Alternative installation (only if `install_github()` fails with HTTP error 401)

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

#### Load libraries

```r
library(qiime2R)
library(tidyverse)
library(vegan)
library(ape)
library(Rcpp)
library(data.table)
library(phyloseq.extended)
```

#### Version check

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

Absolutely 👍 I’ll draft a clean **Markdown/Quarto section** for your website that mirrors your script but presents it as an instructional tutorial with runnable code chunks and nicely formatted outputs.

Here’s a template:

---

### 3. Explore the `phyloseq` Object

In this section we explore different `phyloseq` accessors to get familiar with the contents of our object.

#### Inspect the three core tables

```{r}
head(otu_table(pseq))     # abundance matrix
head(tax_table(pseq))     # taxonomy
head(sample_data(pseq))   # metadata
````

#### Basic information

```{r}
ntaxa(pseq)               # number of taxa
nsamples(pseq)            # number of samples
sample_names(pseq)        # sample IDs
head(taxa_names(pseq))    # taxon (ASV) IDs
```

#### Relabel the ASVs with numbers

```{r}
original_ids <- taxa_names(pseq)
asv_labels <- paste0("ASV", seq_len(ntaxa(pseq)))
taxa_names(pseq) <- asv_labels
tax_table(pseq) <- cbind(tax_table(pseq), OriginalID = original_ids)

head(taxa_names(pseq))
```

#### More accessors

```{r}
sample_sums(pseq)                     # total reads per sample
head(taxa_sums(pseq))                 # total reads per taxon
sample_sums(pseq)[sample_names(pseq)=="IL10-2"]  # depth of one sample
rank_names(pseq)                      # available taxonomic ranks
sample_variables(pseq)                # metadata variables
get_variable(pseq, "type")            # values of a metadata variable
```

#### Extract sample-specific information

```{r}
# Counts for one sample (if taxa are rows, samples are columns)
otu_table(pseq)[, "WT4"]

# Metadata for one sample
sample_data(pseq)["WT4", ]

# Prune the object to only WT4
prune_samples("WT4", pseq)
```

#### Taxonomy queries

```{r}
get_taxa_unique(pseq, "Phylum")
get_taxa_unique(pseq, "Family")
get_taxa_unique(pseq, "Class")
get_taxa_unique(pseq, "Genus")
```

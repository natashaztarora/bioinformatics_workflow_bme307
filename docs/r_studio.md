# R studio

### 1. Package Installation and Library Loading

#### 1. Install BiocManager
```r
install.packages("BiocManager")
````

#### 2. Install Bioconductor packages

```r
BiocManager::install(c("phyloseq", "Biostrings", "S4Vectors", "IRanges", "XVector"))
```

#### 3. Install CRAN packages

```r
install.packages(c("vegan", "ape", "data.table", "Rcpp", "forcats", "tidyverse"))
```

Quick check:

```r
library(phyloseq)
packageVersion("phyloseq")
```

#### 4. Install **qiime2R**

```r
install.packages("remotes")
remotes::install_url(
  "https://github.com/jbisanz/qiime2R/archive/refs/heads/master.zip",
  dependencies = TRUE
)
```

#### 5. Install **phyloseq-extended**

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

#### 6. Load libraries

```r
library(qiime2R)
library(tidyverse)
library(vegan)
library(ape)
library(Rcpp)
library(data.table)
library(phyloseq.extended)
```

#### 7. Version check

```r
pkgs <- c("phyloseq", "qiime2R", "microbiome", "vegan", "ape", "tidyverse")
sapply(pkgs, \(p) as.character(packageVersion(p)))
```


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
<style>
.md-typeset p {
    text-align: justify;
}
</style>

# 2. Exploring the Mouse Gut Microbiome in R

In the previous part of the course, we processed the raw sequencing reads using **QIIME 2**. After quality control, primer trimming, denoising, taxonomic classification and removal of non-target sequences, we obtained a set of **amplicon sequence variants (ASVs)** representing the bacterial communities detected in our samples.

We will now use **R** and the [`phyloseq`](https://joey711.github.io/phyloseq/) package to explore and compare these microbial communities.

Our dataset contains gut microbiome samples from **18 mice**, divided into three groups:

- 6 wild-type (**WT**) mice
- 6 **IL-10-deficient** mice
- 6 **MUC2-deficient** mice

Our main biological question is:

> **Do we observe differences in gut microbial composition and diversity across the three groups of mice?**

To investigate this question, we will explore the data from several complementary perspectives:

1. **Sequencing depth and rarefaction curves** — How deeply was each mouse sequenced, and is the sequencing depth sufficient to compare the samples?
2. **Alpha diversity** — How diverse is the microbial community within each mouse? How many ASVs are detected in each mouse? How evenly are reads distributed among those ASVs?
3. **Normalization** — Convert read counts to relative abundance so that we can compare the relative composition of microbial communities across samples.

4. **Taxonomic composition** — Which bacterial groups are present, and what does their relative composition look like across the mice?

5. **Beta diversity** — How different are the microbial communities between mice and between experimental groups of mice?

Throughout the analysis, remember that there is an important distinction between:

- **Descriptive observations** — patterns that we can see in our data and figures.
- **Statistical evidence** — tests of whether observed differences are larger than expected under a null hypothesis.
- **Biological interpretation** — what these results may mean biologically and what conclusions the experiment allows us to make.

!!! question "Exercise 1 — Before looking at the data"

    In pairs, discuss what gut microbiome differences you might expect among **wild-type, IL-10-deficient and MUC2-deficient mice**.

    Write down **one hypothesis** before examining the data. We will return to your hypothesis at the end of the analysis.


---

## 1. Preparing R

### 1.1 Set the working directory

Open **RStudio** and make sure your working directory is the main `BME307_course_material` folder.

You can check your current working directory with:

```r
getwd()
```

And see the files and folders available there with:

```r
list.files()
```

You should see folders including:

```text
raw_data
metadata
quality_control
reference
qiime2
r_analysis
```

If necessary, set the working directory using:

```r
setwd("PATH/TO/BME307_course_material")
```

!!! tip

    In RStudio, you can also use **Session → Set Working Directory → Choose Directory...**

    Remember that file paths tell R where to find your data, just as they did when we worked in the command line.


### 1.2 Install the required packages

Before starting the analysis, we need to make sure that the required R packages are installed. Most of the packages used in this practical can be installed directly from CRAN. The `phyloseq` package - for microbiome data manipulation - is distributed through *Bioconductor*. First install `BiocManager` if necessary, and then install `phyloseq`:

```r
install.packages(c("tidyverse", "vegan", "ape", "picante"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```

!!! info "Installing versus loading packages"

    A package only needs to be *installed once* on your computer. However, it must be *loaded again every time you start a new R session*.

### 1.3 Load the packages

At the beginning of the analysis, load the packages that we will use:

```r
library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)
library(picante)
```

If these commands run without an error, the packages are installed and ready to use.

---

### 1.4 Load the microbiome dataset

During the QIIME 2 workflow we generated the main components needed for microbiome analysis:

- an **ASV read count table** - referred to as the otu_table
- a **taxonomic classification** table - referred to as the tax_table
- a **sample metadata** table - referred to as the sample_data
- and a **rooted phylogenetic tree** - referred to as the phy_tree

These components can be combined into a single R object using the `phyloseq` package.

For this practical, we have already created this object for you.

A phyloseq object allows us to keep the abundance data, taxonomic information, sample information and phylogenetic relationships together while we analyse the microbiome.


Load the object:

```r
ps <- readRDS("r_analysis/bme307_phyloseq_raw_counts.rds")
```

Now simply type:

```r
ps
```

You should see a summary of the object.

??? info "Expected output"

    The phyloseq object contains:

    - **18 samples**
    - **497 bacterial ASVs**
    - an ASV abundance table
    - sample metadata
    - taxonomic classifications
    - a rooted phylogenetic tree




---

## 2. Exploratory analysis

### 2.1 Explore the phyloseq object

We can access the different components individually:

```r
otu_table(ps)
tax_table(ps)
sample_data(ps)
phy_tree(ps)
```

Because these objects can be large, it is often more useful to inspect only their first few rows:

```r
head(otu_table(ps))
head(tax_table(ps))
head(sample_data(ps))
```

We can also ask simple questions about the dataset:

```r
nsamples(ps)
ntaxa(ps)
sample_names(ps)
sample_sums(ps)
rank_names(ps)
```

??? info "Other useful checks"

    These additional accessors can be useful when you want to inspect a specific part of the object. You do not need to memorize them.

    ```r
    # Look at ASV names and total reads per ASV
    head(taxa_names(ps))
    head(taxa_sums(ps))

    # List the variables stored in the sample metadata
    sample_variables(ps)

    # Extract one metadata variable
    get_variable(ps, varName = "type")

    # Inspect the read depth of one named sample
    sample_sums(ps)["WT4"]

    # Inspect the unique assignments at selected taxonomic ranks
    get_taxa_unique(ps, "Phylum")
    get_taxa_unique(ps, "Genus")
    ```

    Commands such as `head(map)` are not included because `map` is not a component of our current workflow. The metadata are already stored inside `sample_data(ps)`.

!!! question "Exercise 2 — Get to know your dataset"

    Use the commands above to answer the following questions:

    1. How many mice are included?
    2. How many ASVs are present?
    3. What are the four main components stored in the phyloseq object?
    4. What information is contained in the sample metadata?
    5. Which taxonomic ranks are available?
    6. Which mouse has the highest number of reads? Which has the lowest?

    Compare these numbers with the final output of the QIIME 2 workflow. Do they agree?


---

### 2.2 Explore read depth

The number of reads obtained from each sample after quality control and read processing is called its **sequencing depth**.

Different samples rarely contain exactly the same number of reads. Before comparing microbial diversity between samples, it is therefore useful to examine how sequencing depth varies across our dataset.

Remember that a mouse with more sequencing reads has had more opportunities for its microbial DNA sequences to be detected. **It does not necessarily mean that the original sample contained more bacteria.**

For example, a mouse with 80,000 reads should not be interpreted as containing twice as many bacteria as a mouse with 40,000 reads.

Calculate the total number of reads per mouse:

```r
sample_sums(ps)
```

We can turn this into a small data frame:

```r
depth_df <- data.frame(
    sample = sample_names(ps),
    reads = sample_sums(ps),
    group = sample_data(ps)$type
)
```

And visualize it:

```r
ggplot(depth_df, aes(x = sample, y = reads, fill = group)) +
    geom_col() +
    labs(
        x = "Mouse",
        y = "Number of reads",
        fill = "Group"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

!!! question "Exercise 3 — Sequencing depth"

    Look at the sequencing depth across the 18 mice.

    - Is sequencing depth identical across all samples?
    - Which sample has the highest and lowest sequencing depth?
    - Why could differences in sequencing depth become a problem when comparing diversity between samples?
    - Does having more sequencing reads necessarily mean that the original sample contained more bacteria?


!!! example "Homework — Is sequencing depth different among the three mouse groups?"

    We visualized sequencing depth across the 18 mice.

    Now consider a statistical question:

    > **Is there evidence that sequencing depth differs among WT, IL-10-deficient and MUC2-deficient mice?**

    Before choosing a statistical test, think about the structure of the data:

    - The response variable, **sequencing depth**, is numerical.
    - We have **three groups**.
    - The mice are **independent observations**.
    - There are only **six mice per group**.

    First visualize the sequencing depths by group:

    ```r
    ggplot(depth_df, aes(x = group, y = reads)) +
        geom_boxplot() +
        geom_jitter(width = 0.1) +
        theme_minimal()
    ```

    **Questions**

    1. What is the null hypothesis?
    2. Which statistical test could you use?
    3. What is the p-value?
    4. Is there evidence that sequencing depth differs among the three groups?
    5. Why might sequencing depth be important to check before comparing microbiome diversity?

    ??? info "Homework answer — Parametric or non-parametric?"

        A common parametric test for comparing three independent groups is a **one-way ANOVA**. ANOVA tests for differences among group means and assumes that the model residuals are approximately normally distributed and, in its classical form, that the groups have similar variances.

        Fit the model and inspect its residuals:

        ```r
        depth_aov <- aov(reads ~ group, data = depth_df)
        summary(depth_aov)

        qqnorm(residuals(depth_aov))
        qqline(residuals(depth_aov))
        shapiro.test(residuals(depth_aov))
        ```

        A common non-parametric alternative is the **Kruskal-Wallis test**. It is based on ranks and does not assume normally distributed residuals:

        ```r
        kruskal.test(reads ~ group, data = depth_df)
        ```

        With only six mice per group, it is important not to decide whether the data are "normal" from a single test alone. We should also examine the raw observations, the boxplot and the model diagnostics.

        In this dataset, **both tests give a non-significant result and lead to the same practical conclusion**: we do not have sufficient evidence that sequencing depth differs among the three mouse groups.

        The ANOVA null hypothesis is that the three group means are equal. The Kruskal-Wallis null hypothesis is that the groups have the same underlying distribution of ranks.

        A **p-value below the chosen significance threshold**, commonly 0.05, provides evidence against the corresponding null hypothesis. A non-significant result does **not** prove that sequencing depths are identical.


---

### 2.3 Rarefaction curves

We have seen that sequencing depth differs between mice.

Why does this matter for diversity?

A sample with more sequencing reads has had more opportunities to detect **rare ASVs**. Therefore, a deeply sequenced sample could appear richer simply because we looked at it more deeply.

A **rarefaction curve** shows the relationship between sequencing depth and the expected number of ASVs detected:

> **As we sample more and more reads from a microbiome, how many ASVs do we expect to detect?**

At this stage, we are using rarefaction only as a **diagnostic visualization**. Generating a rarefaction curve does not modify or rarefy the phyloseq object.

`vegan` expects samples in rows, so we first prepare the count table:

```r
otu_mat <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
}
```

Generate the curves:

```r
rare <- vegan::rarecurve(
    otu_mat,
    step = 1000,
    tidy = TRUE
)
```

Add the experimental-group information:

```r
rare <- rare %>%
    left_join(
        data.frame(
            Site = sample_names(ps),
            group = sample_data(ps)$type
        ),
        by = "Site"
    )
```

Now visualize the curves:

```r
ggplot(
    rare,
    aes(
        x = Sample,
        y = Species,
        group = Site,
        color = group
    )
) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_vline(
        xintercept = 25000,
        linetype = "dashed"
    ) +
    labs(
        title = "Rarefaction curves of the 18 mouse microbiomes",
        subtitle = "The dashed line marks 25,000 reads",
        x = "Number of sampled reads",
        y = "Expected number of observed ASVs",
        color = "Mouse group"
    ) +
    theme_classic() +
    theme(legend.position = "top")
```

Each curve represents **one mouse**. The curves do not necessarily have the same length because the original samples do not all contain the same number of reads.

!!! question "Exercise 4 — Understand the rarefaction curves"

    - Why do some curves extend farther along the x-axis than others?
    - What happens to the number of detected ASVs as more reads are sampled?
    - What would it mean if a curve begins to flatten?
    - Have the curves begun to flatten around **25,000 reads**?
    - What would happen to a mouse with fewer than 25,000 reads if we later chose to rarefy the data to that depth?


!!! question "Checkpoint — Is our dataset suitable for diversity analysis?"

    Based on the sequencing-depth plot, the statistical comparison of read depth and the rarefaction curves:

    - Are all samples sequenced equally deeply?
    - Have the curves begun to flatten?
    - Do any samples appear insufficiently sequenced?
    - Is sequencing depth systematically associated with mouse group?
    - What limitations should we remember when comparing diversity?

---

## Microbial diversity analysis

We will now examine diversity from several complementary perspectives:

- **Alpha diversity** — diversity within each individual mouse;
- **Taxonomic composition** — the relative composition of bacterial groups;
- **Beta diversity** — differences in complete community composition between mice.


## 3. Alpha diversity

**Alpha diversity** describes diversity **within a single microbial community**.

### 3.1 Before we look at alpha diversity...

As we saw from the rarefaction curves, samples with more sequencing reads have more opportunities to detect rare ASVs. This can influence estimates of alpha diversity. To compare diversity across mice at the same sampling depth, we will therefore **rarefy the dataset to the lowest read depth in our dataset**.

Rarefaction randomly subsamples the same number of reads from each sample. Note that this leads to comparable sampling effort, but it also discards reads and may remove rare ASVs.

```r
set.seed(123)

ps_rare <- rarefy_even_depth(
    ps,
    sample.size = LOWEST_READ_DEPTH,
    rngseed = 123,
    replace = FALSE,
    verbose = FALSE
)

This creates a new phyloseq object, ps_rare, in which all samples have the same number of reads. The original ps object remains unchanged.

To measure diversity, there are different metrics available; each of these captures different properties of a microbial community. 

| Metric | What does it consider? | Main question |
| --- | --- | --- |
| **Observed ASVs** | Richness | How many different ASVs were detected? |
| **Pielou's evenness** | Evenness | How evenly are reads distributed among the detected ASVs? |
| **Shannon diversity** | Richness + evenness | How diverse is the community considering both the number and distribution of ASVs? |
| **Faith's PD** | Phylogenetic relationships | How much phylogenetic diversity is represented? |

First calculate Observed richness and Shannon diversity from the rarefied object (ps_rare), then add the sample metadata:

```r
alpha <- estimate_richness(
    ps_rare,
    measures = c("Observed", "Shannon")
)

alpha$sample <- sample_names(ps_rare)

metadata <- data.frame(sample_data(ps_rare))
alpha$type <- metadata[alpha$sample, "type"]

table(alpha$type, useNA = "ifany")
```

Pielou's evenness is derived from Shannon diversity and Observed richness:

```r
alpha$Pielou <- alpha$Shannon / log(alpha$Observed)
```

Although Pielou is presented before Shannon below to match the teaching sequence, its calculation necessarily uses the Shannon value.

### 3.2 Observed ASVs

Observed ASVs is a measure of **richness**: the number of different ASVs detected in each mouse. It does not consider how abundant or evolutionarily distinct those ASVs are.

```r
ggplot(alpha, aes(x = type, y = Observed, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1) +
    labs(
        x = "Mouse group",
        y = "Observed ASVs"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
```

### 3.3 Pielou's evenness

Pielou's evenness focuses on how evenly reads are distributed among the detected ASVs. Values closer to 1 indicate a more even community.

```r
ggplot(alpha, aes(x = type, y = Pielou, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1) +
    labs(
        x = "Mouse group",
        y = "Pielou's evenness"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
```

### 3.4 Shannon diversity

Shannon diversity combines richness and evenness. It increases when more ASVs are present and when their abundances are more evenly distributed.

```r
ggplot(alpha, aes(x = type, y = Shannon, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1) +
    labs(
        x = "Mouse group",
        y = "Shannon diversity"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
```

### 3.5 Faith's phylogenetic diversity

Observed richness treats every ASV as a different feature. Faith's phylogenetic diversity additionally uses the **phylogenetic tree** and measures the total branch length represented in each community.

```r
otu_alpha <- as(otu_table(ps_rare), "matrix")

if (taxa_are_rows(ps_rare)) {
    otu_alpha <- t(otu_alpha)
}

faith <- picante::pd(
    otu_alpha,
    phy_tree(ps_rare),
    include.root = TRUE
)

alpha$Faith_PD <- faith[alpha$sample, "PD"]
```

```r
ggplot(alpha, aes(x = type, y = Faith_PD, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1) +
    labs(
        x = "Mouse group",
        y = "Faith's phylogenetic diversity"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
```

!!! question "Exercise 5 — Explore alpha diversity"

    Compare the three mouse groups using the four metrics.

    - What patterns do you observe for **Observed ASVs**?
    - What patterns do you observe for **Pielou's evenness**?
    - What patterns do you observe for **Shannon diversity**?
    - What patterns do you observe for **Faith's PD**?
    - Do all four metrics show the same pattern?
    - Why might different metrics give different views of the same microbial community?


### 3.6 Statistical analysis for alpha diversity

For alpha diversity, we are comparing diversity values among **three independent groups of mice**. We will use the **Kruskal-Wallis test**.

```r
kruskal.test(Observed ~ type, data = alpha)
kruskal.test(Pielou ~ type, data = alpha)
kruskal.test(Shannon ~ type, data = alpha)
kruskal.test(Faith_PD ~ type, data = alpha)
```

The Kruskal-Wallis test asks whether there is evidence that the distributions differ among the three groups. If an overall test provides evidence of a difference, we can explore pairwise comparisons. For example:

```r
pairwise.wilcox.test(
    alpha$Shannon,
    alpha$type,
    p.adjust.method = "BH"
)
```

!!! question "Exercise 6 — Alpha-diversity statistics"

    - Which visual patterns are supported by the statistical tests?
    - Are there patterns that looked different visually but are not supported statistically?
    - If an overall test suggests a difference, which pairwise comparisons would you investigate?

---
## 4. Normalize the data

So far, our ASV table contains the original read counts for each sample.

For analyses where we want to compare the **relative abundance of ASVs or bacterial taxa**, we convert the read counts in each sample to relative abundance.

This creates a new phyloseq object, `ps_norm`, which we will use for the next analyses:

```r
ps_norm <- transform_sample_counts(
    ps,
    function(x) x / sum(x)
)
```

Each value now represents the **proportion of reads within that sample** assigned to an ASV. The values in each sample sum to 1 rather than to the original sequencing depth.

Check the result:

```r
sample_sums(ps_norm)
```
!!! warning "Relative abundance ≠ absolute abundance"

    If a bacterial group represents 20% of the sequencing reads in a sample, this does **not** tell us the absolute number of bacterial cells that were present.

    Relative-abundance normalization also does not remove the compositional nature of microbiome data: when the relative abundance of one taxon increases, the proportions of other taxa must collectively decrease.



## 5. Taxonomic composition

We will now examine **which bacterial groups are present** in our samples. We can now group ASVs according to their taxonomy. We will start at the **phylum level**.

### 5.1 Group ASVs at the phylum level

```r
ps_phylum <- tax_glom(
    ps_norm,
    taxrank = "Phylum",
    NArm = FALSE
)
```

`tax_glom()` combines ASVs that have the same taxonomic lineage up to the selected rank. Their abundances are summed, so the object now represents phylum-level groups rather than individual ASVs. Two ASVs carrying the same lower-rank label are only merged when their classifications at the preceding ranks also agree.

We use `NArm = FALSE` so that sequences without a Phylum assignment are not silently discarded. Remember that placeholder labels such as `"uncultured"` are not automatically well-defined biological taxa. This is why the complete taxonomic lineages must always be inspected.


### 5.2 Visualize the bacterial composition

Because `ps_phylum` was created from `ps_norm`, the y-axis represents **relative abundance**:

```r
plot_bar(ps_phylum, x = "sample_id", fill = "Phylum") +
    labs(
        x = "Mouse",
        y = "Relative abundance",
        fill = "Phylum"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Each bar represents **one mouse**.

!!! question "Exercise 7 — Explore taxonomic composition"

    Look at the figure before trying to draw conclusions.

    - Which bacterial phyla appear to dominate the samples?
    - Do mice from the same experimental group appear similar?
    - Do you notice any patterns among WT, IL-10-deficient and MUC2-deficient mice?
    - Is there variation between individual mice?
    - Would looking at a different taxonomic level change what you can see?

    At this stage, simply describe what you **observe**.


### 5.3 Explore other taxonomic levels


You are not restricted to phylum.

For example, you could replace:

```r
taxrank = "Phylum"
```

with another available taxonomic rank.

Check the available ranks with:

```r
rank_names(ps)
```

!!! tip "Make the plot your own"

    The plots in this practical are deliberately simple.

    For your poster, you can explore different labels and themes using [`ggplot2`](https://ggplot2.tidyverse.org/).

    For example:

    ```r
    + theme_classic()
    ```

    or:

    ```r
    + labs(title = "Gut microbiome composition")
    ```

    First make sure you understand the analysis. **Then** improve the appearance of the figure.


!!! question "Try it — Explore different taxonomic ranks"

    Repeat the taxonomic composition analysis using other available taxonomic ranks.

    For example, try:

    - Domain
    - Class
    - Order
    - Family
    - Genus

    As you move towards more specific taxonomic ranks, look carefully at the legend.

    **What happens to the number of different taxonomic groups?**


!!! tip "Too many taxa to visualize?"

    At lower taxonomic ranks, especially **Family** and **Genus**, there may be too many taxa to display clearly.

    One option is to show only the **15 most abundant taxa** and combine everything else into a category called **Other**.

    ??? example "Show code: Keep the top 15 genera and group the rest as Other"

        First, agglomerate ASVs that share the same Genus assignment:

        ```r
        ps_genus <- tax_glom(
            ps_norm,
            taxrank = "Genus",
            NArm = FALSE
        )
        ```

        !!! note "What does `NArm = FALSE` mean?"

            `NArm` controls what happens to ASVs that **do not have an assignment at the selected taxonomic rank**.

            - `NArm = TRUE` removes taxa for which the selected rank is `NA`.
            - `NArm = FALSE` keeps them.

            We use `NArm = FALSE` here because we do not want to silently discard sequences simply because they could not be classified to Genus level.

            This is especially important at lower taxonomic ranks, where some ASVs may not have sufficiently confident assignments.

        The object is already in relative abundance because it was created from `ps_norm`:

        ```r
        ps_genus_rel <- ps_genus
        ```

        Identify the 15 most abundant genus-level groups across the complete dataset:

        ```r
        top15 <- names(
            sort(taxa_sums(ps_genus_rel), decreasing = TRUE)
        )[1:15]
        ```

        Save their complete taxonomic classification **before** changing any labels:

        ```r
        top15_taxonomy <- as.data.frame(
            tax_table(ps_genus_rel)[top15, ]
        )

        top15_taxonomy
        ```

        This table allows you to inspect the complete lineage of each of the 15 most abundant groups:

        **Domain → Phylum → Class → Order → Family → Genus**

        Now rename everything outside the top 15 as `"Other"`:

        ```r
        tax_table(ps_genus_rel)[
            !(taxa_names(ps_genus_rel) %in% top15),
            "Genus"
        ] <- "Other"
        ```

        Agglomerate the relabelled object again:

        ```r
        ps_genus_top15 <- tax_glom(
            ps_genus_rel,
            taxrank = "Genus",
            NArm = FALSE
        )
        ```

        `tax_glom()` respects the complete lineage up to Genus, so rows from different higher-level lineages may remain separate internally even though they are all labelled `"Other"`. In the stacked bar plot, however, all rows carrying that label are displayed using the same `"Other"` fill category.

        Finally, plot the result:

        ```r
        plot_bar(
            ps_genus_top15,
            x = "sample_id",
            fill = "Genus"
        ) +
            labs(
                x = "Mouse",
                y = "Relative abundance",
                fill = "Genus"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1)
            )
        ```

        The legend now shows the **15 most abundant genus-level groups plus "Other"**.

        To inspect their full taxonomic lineages again:

        ```r
        top15_taxonomy
        ```

    Try changing `"Genus"` to another taxonomic rank. What happens as you move from Phylum → Class → Order → Family → Genus?



!!! example "Homework — Taxonomic resolution"

    **1. Look at your taxonomic composition plots at different ranks.**

    Examine the legends from **Phylum → Class → Order → Family → Genus**.

    Do all taxonomic assignments reach the same level of resolution?

    **2. How many of our 497 ASVs are classified at each taxonomic rank?**

    Determine how many ASVs have an assignment at:

    - Phylum
    - Class
    - Order
    - Family
    - Genus

    What happens as you move towards more specific taxonomic ranks?

??? success "Homework answer — How well were our ASVs classified?"

    There are two related but different questions we can ask:

    1. **How many ASVs received an assignment at each taxonomic rank?**
    2. **How many distinct taxonomic groups are represented at each rank?**

    These are not the same thing. Many different ASVs can receive the same taxonomic assignment.

    ```r
    tax_df <- as.data.frame(tax_table(ps))

    resolution <- data.frame(
        rank = colnames(tax_df),

        assigned_ASVs = sapply(
            tax_df,
            function(x) sum(!is.na(x) & x != "")
        ),

        distinct_groups = sapply(
            tax_df,
            function(x) length(unique(x[!is.na(x) & x != ""]))
        )
    )

    resolution$total_ASVs <- ntaxa(ps)

    resolution$percent_assigned <- round(
        100 * resolution$assigned_ASVs / resolution$total_ASVs,
        1
    )

    resolution
    ```

    For this dataset, the number of ASVs assigned at each rank is:

    | Rank | Assigned ASVs | Total ASVs | % assigned |
    | --- | ---: | ---: | ---: |
    | Domain | 497 | 497 | 100.0 |
    | Phylum | 496 | 497 | 99.8 |
    | Class | 496 | 497 | 99.8 |
    | Order | 494 | 497 | 99.4 |
    | Family | 484 | 497 | 97.4 |
    | Genus | 441 | 497 | 88.7 |

    The `distinct_groups` column tells us something different: **how many different non-missing taxonomic labels occur at each rank**.

    For example, although **441 ASVs have a Genus-level assignment**, they correspond to only **64 distinct non-missing Genus labels**.

    Therefore:

    **441 assigned ASVs ≠ 441 genera**

    Multiple ASVs can have the same taxonomic classification. For example, several different ASVs may all be classified as *Blautia*.

    This is also why `tax_glom(ps, taxrank = "Genus")` changes the structure of the dataset: ASVs sharing the same taxonomic lineage through Genus are combined into the same taxonomic group.

    As we move towards more specific taxonomic ranks, fewer ASVs are usually classified.

    This reflects the limited taxonomic resolution of short-read 16S sequencing and the information available in the reference database.

    Note that a non-missing taxonomic label does not necessarily mean that the organism has been identified precisely. Labels such as `uncultured` may still occur.


!!! warning "Always inspect your taxonomic assignments"

    A non-missing taxonomic assignment does **not necessarily mean that we have a clean, informative biological name**.

    After summarizing classification success, always inspect the actual taxonomic labels.

    ??? example "Inspect the taxonomic assignments"

        For example, look at every unique Genus-level assignment:

        ```r
        unique(tax_df$Genus)
        ```

        Count the distinct non-missing Genus-level labels:

        ```r
        length(
            unique(
                tax_df$Genus[!is.na(tax_df$Genus) & tax_df$Genus != ""]
            )
        )
        ```

        You can repeat this for other ranks:

        ```r
        unique(tax_df$Phylum)
        unique(tax_df$Class)
        unique(tax_df$Order)
        unique(tax_df$Family)
        unique(tax_df$Genus)
        ```

        Or count the number of distinct groups at every rank at once:

        ```r
        sapply(
            tax_df,
            function(x) length(unique(x[!is.na(x) & x != ""]))
        )
        ```

    When inspecting the results, you may encounter labels such as:

    - `NA`
    - `"uncultured"`
    - names ending in `_group`
    - names based on uncultured or incompletely characterized lineages
    - family-level names appearing in the Genus column
    - database-specific placeholder names

    These labels are important to notice.

    For example, in our Genus column we find labels such as `"uncultured"`, `"Lachnospiraceae_NK4A136_group"`, `"Muribaculaceae"` and `"Rikenellaceae_RC9_gut_group"`.

    Therefore, **"assigned at Genus level" does not automatically mean "identified as a well-characterized named genus."**

    Taxonomic classification depends on the reference database, the sequenced 16S region, sequence similarity and the confidence of the classifier. Always inspect the actual assignments before interpreting the biological resolution of your dataset.

---

## 6. Beta diversity

Alpha diversity asks: **How diverse is each individual mouse?**

Beta diversity asks: **How different are the microbial communities between mice?** 

For our analyses, we will use the normalized object, `ps_norm`.


| Metric | Uses abundance? | Uses phylogeny? | Main emphasis |
| --- | --- | --- | --- |
| **Bray-Curtis** | Yes | No | Differences in ASV abundance |
| **Unweighted UniFrac** | No | Yes | Presence/absence of phylogenetic lineages |
| **Weighted UniFrac** | Yes | Yes | Abundance-weighted phylogenetic differences |

### 6.1 Calculate the distance matrices

**Bray-Curtis** considers differences in abundance but not phylogenetic relationships:

```r
bray <- phyloseq::distance(
    ps_norm,
    method = "bray"
)
```

**Unweighted UniFrac** considers presence/absence together with phylogenetic relationships:

```r
unifrac_unweighted <- UniFrac(
    ps_norm,
    weighted = FALSE
)
```

**Weighted UniFrac** combines abundance with phylogenetic relationships:

```r
unifrac_weighted <- UniFrac(
    ps_norm,
    weighted = TRUE
)
```

### 6.2 Visualize Bray-Curtis with PCoA

A distance matrix contains pairwise differences between all samples. **Principal Coordinates Analysis (PCoA)** represents these relationships in fewer dimensions.

```r
ord_bray <- ordinate(
    ps_norm,
    method = "PCoA",
    distance = bray
)

plot_ordination(
    ps_norm,
    ord_bray,
    color = "type"
) +
    geom_point(size = 4) +
    theme_minimal()
```

Each point represents **one mouse**.

### 6.3 Visualize the other beta-diversity metrics


#### Unweighted UniFrac

```r
ord_unweighted <- ordinate(
    ps_norm,
    method = "PCoA",
    distance = unifrac_unweighted
)

plot_ordination(ps_norm, ord_unweighted, color = "type") +
    geom_point(size = 4) +
    theme_minimal()
```

#### Weighted UniFrac

```r
ord_weighted <- ordinate(
    ps_norm,
    method = "PCoA",
    distance = unifrac_weighted
)

plot_ordination(ps_norm, ord_weighted, color = "type") +
    geom_point(size = 4) +
    theme_minimal()
```

!!! question "Exercise 8 — Compare the beta-diversity metrics"

    - Do all three metrics show the same pattern?
    - Does considering **abundance** change what you see?
    - Does considering **phylogenetic relationships** change what you see?
    - Is within-group variability similar for all three mouse groups?

!!! tip "There is no universally 'best' beta-diversity metric"

    These metrics use different information and answer slightly different ecological questions. Do not choose a metric simply because its PCoA produces the clearest-looking separation.


### 6.4 Statistical analysis for beta diversity

We can use **PERMANOVA** (*Permutational Multivariate Analysis of Variance*) to ask whether microbial community composition is associated with mouse group.

```r
meta_beta <- data.frame(sample_data(ps_norm))
meta_beta <- meta_beta[labels(bray), , drop = FALSE]

stopifnot(identical(rownames(meta_beta), labels(bray)))
```

```r
adonis2(bray ~ type, data = meta_beta, permutations = 999)
adonis2(unifrac_unweighted ~ type, data = meta_beta, permutations = 999)
adonis2(unifrac_weighted ~ type, data = meta_beta, permutations = 999)
```

Examine both the **p-value** and the **R² value**, which describes the proportion of variation in the distance matrix associated with mouse group.

!!! question "Exercise 9 — PERMANOVA"

    - Is there statistical evidence that community composition is associated with mouse group?
    - What proportion of variation is associated with group?
    - Do all three beta-diversity metrics give the same result?
    - How do the statistical results compare with the PCoA plots?

??? info "Optional — PERMDISP"

    PERMANOVA can be influenced by differences in **within-group dispersion**. We can investigate this using **PERMDISP**. For Bray-Curtis:

    ```r
    dispersion_bray <- betadisper(bray, meta_beta$type)
    plot(dispersion_bray)
    permutest(dispersion_bray, permutations = 999)
    ```

    You can repeat this for another distance matrix by replacing `bray`.

---

## 7. Putting everything together

You have now looked at the microbiome from several different perspectives:

**Sequencing depth and rarefaction curves**

> Did we sequence the samples deeply enough to make meaningful diversity comparisons?

**Alpha diversity**

> How diverse is the microbial community within each mouse?

**Taxonomic composition**

> Which bacterial groups were detected, and what patterns can we see across mice?

**Beta diversity**

> How different are the complete microbial communities between mice?


!!! question "Exercise 10 — Final interpretation"

    Return to the prediction or hypothesis you wrote at the beginning.

    Discuss your results in pairs.

    - What patterns did you observe in taxonomic composition?
    - What did the alpha-diversity metrics show?
    - What did the beta-diversity metrics show?
    - Did all metrics tell the same story?
    - Which observations were supported statistically?
    - How much variation was present between individual mice?
    - Does the evidence support your original prediction?
    - What additional information or experiments would you want before making stronger biological conclusions?

    There is not necessarily one single observation that summarizes the complete dataset. Try to combine the different pieces of evidence.


---

## 8. Important limitations

Before interpreting the experiment, consider some of its limitations.

!!! warning "Things to keep in mind"

    **Sample size**

    There are only six mice per group. How might this affect the conclusions that can be drawn?

    **Relative abundance**

    16S rRNA gene sequencing primarily provides information about the relative composition of the microbial community. What information about absolute bacterial abundance is missing?

    **Rarefaction**

    Rarefaction (subsampling) standardizes sequencing depth, but it also discards reads. What are the consequences of that trade-off?

    **Taxonomic resolution**

    Short-read 16S rRNA gene sequencing may not reliably distinguish microorganisms at species or strain level. How does this limit biological interpretation?

    **Taxonomy is not function**

    Knowing which bacterial taxa are present does not automatically tell us what those microorganisms are doing.

    **Association is not causation**

    Think carefully about which conclusions are supported by this experiment and which would require additional experiments.


---

## 9. Going further with your figures

For this practical, the plotting code has deliberately been kept simple so that the **biological analysis remains the main focus**.

For your final poster, however, you may want to improve the appearance of your figures.

Most of the figures in this practical use `ggplot2`, which means you can modify them by adding additional layers.

For example:

```r
+ theme_classic()
```

Change labels:

```r
+ labs(
    title = "Your title",
    x = "Your x-axis label",
    y = "Your y-axis label"
)
```

Modify text:

```r
+ theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
)
```

You can explore additional options in:

- [`ggplot2` documentation](https://ggplot2.tidyverse.org/)
- [`phyloseq` documentation](https://joey711.github.io/phyloseq/)
- [`vegan` documentation](https://vegandevs.github.io/vegan/)

!!! tip "A good plotting workflow"

    1. Make sure the analysis is scientifically appropriate.
    2. Create a simple figure.
    3. Understand what the figure is showing.
    4. Only then modify labels, themes, colours and layout for your poster.

    **Do not choose an analysis method simply because it produces the prettiest plot.**


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

1. **Sequencing depth** — How deeply was each mouse sequenced, and is the sequencing depth sufficient to compare the samples?
2. **Taxonomic composition** — What does the microbial composition look like across the mice?
3. **Alpha diversity** — How diverse is the microbial community within each mouse? How many ASVs are detected in each mouse? How evenly are reads distributed among those ASVs?
4. **Beta diversity** — How different are the microbial communities between mice and between experimental groups of mice?

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

Before starting the analysis, we need to make sure that the required R packages are installed. Most of the packages used in this practical can be installed directly from CRAN. The ⁠ phyloseq ⁠ package - for microbiome data manipulation - is distributed through *Bioconductor*. First install ⁠ BiocManager ⁠ if necessary, and then install ⁠ phylosec:

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

## 2. Loading the microbiome dataset

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


### 2.1 Exploring the phyloseq object

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

## 3. Sequencing depth

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


---

## 4. Taxonomic composition

We will now examine **which bacterial groups are present** in our samples.

The ASV table contains hundreds of individual ASVs. Looking at all of them simultaneously would make the overall composition difficult to interpret.

We can therefore group ASVs according to their taxonomy. Here, we will start at the **phylum level**.


### 4.1 Group ASVs at the phylum level

```r
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
```


### 4.2 Convert counts to relative abundance

For taxonomic composition plots, it is common to convert the counts into **relative abundance**:

```r
ps_phylum_rel <- transform_sample_counts(
    ps_phylum,
    function(x) x / sum(x)
)
```

Each taxon is now represented as a **proportion of the reads within that sample**.

!!! warning "Relative abundance ≠ absolute abundance"

    If a bacterial group represents 20% of the sequencing reads in a sample, this does **not** tell us the absolute number of bacterial cells that were present.


### 4.3 Visualize the bacterial composition

```r
plot_bar(ps_phylum_rel, x = "sample_id", fill = "Phylum") +
    labs(
        x = "Mouse",
        y = "Relative abundance",
        fill = "Phylum"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Each bar represents **one mouse**.

!!! question "Exercise 4 — Explore taxonomic composition"

    Look at the figure before trying to draw conclusions.

    - Which bacterial phyla appear to dominate the samples?
    - Do mice from the same experimental group appear similar?
    - Do you notice any patterns among WT, IL-10-deficient and MUC2-deficient mice?
    - Is there variation between individual mice?
    - Would looking at a different taxonomic level change what you can see?

    At this stage, simply describe what you **observe**.


### 4.4 Explore other taxonomic levels

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


---

## 5. Rarefaction

We have seen that sequencing depth differs between mice.

Why does this matter for diversity?

A sample with more sequencing reads has had more opportunities to detect **rare ASVs**. Therefore, a deeply sequenced sample could appear richer simply because we looked at it more deeply.

A **rarefaction curve** allows us to explore the relationship between sequencing depth and the number of ASVs detected.

The question is:

> **As we sample more and more reads from a microbiome, how many ASVs do we expect to detect?**


### 5.1 Prepare the ASV count matrix

`vegan` expects samples in rows, so we first prepare the count table:

```r
otu_mat <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
}
```


### 5.2 Generate the rarefaction curves

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

Each curve represents **one mouse**.

The curves do not necessarily have the same length because the original samples do not all contain the same number of sequencing reads.

!!! question "Exercise 5 — Understand the rarefaction curves"

    Look carefully at the curves.

    - Why do some curves extend farther along the x-axis than others?
    - What happens to the number of detected ASVs as more reads are sampled?
    - What would it mean if a curve begins to flatten?
    - Have the curves begun to flatten around **25,000 reads**?
    - Would every mouse be able to reach **50,000 reads**?
    - What would happen to a mouse with fewer reads than the rarefaction depth we choose?

    Based on the curves, what sequencing depth would **you** consider reasonable for comparing these samples?


### 5.3 Choosing a rarefaction depth

Choosing a rarefaction depth involves a trade-off.

At a **higher depth**, more sequencing information is retained from deeply sequenced samples, but low-depth samples may need to be excluded.

At a **lower depth**, more samples can be retained, but more sequencing reads are discarded from deeply sequenced samples.

Let's compare several possible depths:

```r
depths <- c(10000, 25000, 30000, 50000, 88500)

data.frame(
    depth = depths,
    mice_retained = sapply(
        depths,
        function(x) sum(sample_sums(ps) >= x)
    )
)
```

!!! question "Exercise 6 — Choose a rarefaction depth"

    Consider both the sequencing-depth plot and the rarefaction curves.

    - How many mice would remain at each depth?
    - Would all three experimental groups remain represented?
    - At what point do the curves begin to flatten?
    - What do we gain by choosing a higher depth?
    - What do we lose?

    **Which depth would you choose, and why?**

For the rest of this practical, we will use a rarefaction depth of **25,000 reads per mouse**.


### 5.4 Rarefy the dataset

```r
ps_rarefied <- rarefy_even_depth(
    ps,
    sample.size = 25000,
    rngseed = 42,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
)
```

Check the resulting object:

```r
ps_rarefied
sample_sums(ps_rarefied)
```

!!! question "Exercise 7 — What did rarefaction do?"

    Compare `ps` and `ps_rarefied`.

    - How many mice remain?
    - How many reads does each mouse now contain?
    - Did the number of ASVs change?
    - Why might some ASVs disappear during rarefaction?
    - What is an advantage of rarefaction?
    - What is a disadvantage?


---

## 6. Alpha diversity

**Alpha diversity** describes diversity **within a single microbial community**.

There is no single definition of diversity. Different metrics capture different properties of a microbial community.

We will focus on three complementary metrics:

| Metric | What does it consider? | Main question |
| --- | --- | --- |
| **Observed ASVs** | Richness | How many different ASVs were detected? |
| **Shannon diversity** | Richness + evenness | How diverse is the community considering both the number and distribution of ASVs? |
| **Faith's PD** | Phylogenetic relationships | How much phylogenetic diversity is represented? |


### 6.1 Observed ASVs and Shannon diversity

```r
alpha <- estimate_richness(
    ps_rarefied,
    measures = c("Observed", "Shannon")
)
```

Add sample identifiers and metadata:

```r
alpha$sample <- sample_names(ps_rarefied)

metadata <- data.frame(sample_data(ps_rarefied))

alpha$type <- metadata$type

table(alpha$type, useNA = "ifany")
```


### 6.2 Faith's phylogenetic diversity

Observed richness treats every ASV as a different feature.

Faith's phylogenetic diversity additionally uses the **phylogenetic tree** to consider the evolutionary relationships among the ASVs.

```r
otu_alpha <- as(otu_table(ps_rarefied), "matrix")

if (taxa_are_rows(ps_rarefied)) {
    otu_alpha <- t(otu_alpha)
}

faith <- picante::pd(
    otu_alpha,
    phy_tree(ps_rarefied),
    include.root = TRUE
)

alpha$Faith_PD <- faith$PD
```


### 6.3 Visualize alpha diversity

For example, Shannon diversity:

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

Now use the same structure to explore:

```text
Observed
Faith_PD
```

!!! question "Exercise 8 — Explore alpha diversity"

    Compare the three mouse groups using the different alpha-diversity metrics.

    - What patterns do you observe for **Observed ASVs**?
    - What patterns do you observe for **Shannon diversity**?
    - What patterns do you observe for **Faith's PD**?
    - Do all three metrics show the same pattern?
    - Why might different metrics give different views of the same microbial community?

    Do not worry about statistical significance yet. We will test the patterns later.


??? info "Optional — Pielou's evenness"

    Another alpha-diversity metric is **Pielou's evenness**, which focuses specifically on how evenly reads are distributed among the detected ASVs.

    If you want to explore it:

    ```r
    alpha$Pielou <- alpha$Shannon / log(alpha$Observed)

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

    How does the pattern compare with Shannon diversity?


---

## 7. Beta diversity

Alpha diversity asks:

> **How diverse is each individual mouse?**

Beta diversity asks:

> **How different are the microbial communities between mice?**

There are several ways to define the difference between two microbial communities.

We will compare four commonly used metrics:

| Metric | Uses abundance? | Uses phylogeny? |
| --- | --- | --- |
| **Jaccard** | No | No |
| **Bray-Curtis** | Yes | No |
| **Unweighted UniFrac** | No | Yes |
| **Weighted UniFrac** | Yes | Yes |

This gives us two useful questions when thinking about beta diversity:

> **Do we care only about whether an ASV was detected, or also about how abundant it was?**

> **Do we want to consider the evolutionary relationships among the ASVs?**


### 7.1 Calculate the distance matrices

**Jaccard** focuses on presence/absence:

```r
jaccard <- phyloseq::distance(
    ps_rarefied,
    method = "jaccard",
    binary = TRUE
)
```

**Bray-Curtis** also considers abundance:

```r
bray <- phyloseq::distance(
    ps_rarefied,
    method = "bray"
)
```

**Unweighted UniFrac** incorporates presence/absence and phylogenetic relationships:

```r
unifrac_unweighted <- UniFrac(
    ps_rarefied,
    weighted = FALSE
)
```

**Weighted UniFrac** additionally incorporates abundance:

```r
unifrac_weighted <- UniFrac(
    ps_rarefied,
    weighted = TRUE
)
```


### 7.2 Visualize beta diversity with PCoA

A distance matrix contains the pairwise differences between all samples. To visualize these relationships, we can use **Principal Coordinates Analysis (PCoA)**.

We will first work through **Bray-Curtis** together:

```r
ord_bray <- ordinate(
    ps_rarefied,
    method = "PCoA",
    distance = bray
)
```

Plot it:

```r
plot_ordination(
    ps_rarefied,
    ord_bray,
    color = "type"
) +
    geom_point(size = 4) +
    theme_minimal()
```

Each point represents **one mouse**.

!!! question "Exercise 9 — Explore the Bray-Curtis PCoA"

    Look at the distribution of the samples.

    - Do some mice appear closer together than others?
    - Do you observe any clustering by mouse group?
    - Is there variation among mice belonging to the same group?
    - Are there any samples that appear unusual?

    Again, describe what you **see** before testing it statistically.


### 7.3 Compare beta-diversity metrics

Now repeat the PCoA using the other distance metrics.

#### Jaccard

```r
ord_jaccard <- ordinate(
    ps_rarefied,
    method = "PCoA",
    distance = jaccard
)

plot_ordination(
    ps_rarefied,
    ord_jaccard,
    color = "type"
) +
    geom_point(size = 4) +
    theme_minimal()
```

#### Unweighted UniFrac

```r
ord_unweighted <- ordinate(
    ps_rarefied,
    method = "PCoA",
    distance = unifrac_unweighted
)

plot_ordination(
    ps_rarefied,
    ord_unweighted,
    color = "type"
) +
    geom_point(size = 4) +
    theme_minimal()
```

#### Weighted UniFrac

```r
ord_weighted <- ordinate(
    ps_rarefied,
    method = "PCoA",
    distance = unifrac_weighted
)

plot_ordination(
    ps_rarefied,
    ord_weighted,
    color = "type"
) +
    geom_point(size = 4) +
    theme_minimal()
```

!!! question "Exercise 10 — Compare the beta-diversity metrics"

    Compare the four PCoA plots.

    - Do all four metrics show the same pattern?
    - Does considering **abundance** change what you see?
    - Does considering **phylogenetic relationships** change what you see?
    - Is within-group variability similar for all three mouse groups?

    Think about **why** the metrics might produce different patterns.

!!! tip "There is no universally 'best' beta-diversity metric"

    These metrics use different information and therefore answer slightly different ecological questions.

    Do not choose a metric simply because its PCoA produces the clearest-looking separation.


---

# 8. Statistical analysis

So far, we have deliberately focused on **exploring and visualizing the data**.

Only now will we ask whether some of the patterns we observed are supported statistically.


## 8.1 Alpha-diversity statistics

For alpha diversity, we are comparing diversity values among **three independent groups of mice**.

We will use the **Kruskal-Wallis test**.

### Observed ASVs

```r
kruskal.test(Observed ~ type, data = alpha)
```

### Shannon diversity

```r
kruskal.test(Shannon ~ type, data = alpha)
```

### Faith's PD

```r
kruskal.test(Faith_PD ~ type, data = alpha)
```

The Kruskal-Wallis test asks whether there is evidence that the distributions differ among the three groups.

If an overall test provides evidence of a difference, we can explore the pairwise comparisons.

For example:

```r
pairwise.wilcox.test(
    alpha$Shannon,
    alpha$type,
    p.adjust.method = "BH"
)
```

!!! question "Exercise 11 — Alpha-diversity statistics"

    Compare your statistical results with the figures you examined earlier.

    - Which visual patterns are supported by the statistical tests?
    - Are there patterns that looked different visually but are not supported statistically?
    - If an overall test suggests a difference, which pairwise comparisons would you investigate?


## 8.2 Beta-diversity statistics: PERMANOVA

For beta diversity, our data are represented by **distance matrices** rather than one diversity value per mouse.

We can use **PERMANOVA** (*Permutational Multivariate Analysis of Variance*) to ask whether microbial community composition is associated with mouse group.

Create a metadata data frame:

```r
meta <- data.frame(sample_data(ps_rarefied))
```

### Bray-Curtis

```r
adonis2(
    bray ~ type,
    data = meta,
    permutations = 999
)
```

We can repeat the analysis for the other beta-diversity metrics:

```r
adonis2(
    jaccard ~ type,
    data = meta,
    permutations = 999
)

adonis2(
    unifrac_unweighted ~ type,
    data = meta,
    permutations = 999
)

adonis2(
    unifrac_weighted ~ type,
    data = meta,
    permutations = 999
)
```

Two useful values to examine are:

- the **p-value**;
- the **R² value**, which describes the proportion of variation in the distance matrix associated with the grouping variable.

!!! question "Exercise 12 — PERMANOVA"

    Compare the PERMANOVA results with your PCoA plots.

    - Is there statistical evidence that community composition is associated with mouse group?
    - What proportion of variation is associated with group?
    - Do all four beta-diversity metrics give the same result?
    - How do the statistical results compare with what you observed visually?


??? info "Optional — PERMDISP"

    PERMANOVA can be influenced by differences in **within-group dispersion**.

    In other words, one group may contain microbiomes that are much more variable than another group.

    We can investigate this using **PERMDISP**.

    For Bray-Curtis:

    ```r
    dispersion_bray <- betadisper(
        bray,
        meta$type
    )

    plot(dispersion_bray)

    permutest(
        dispersion_bray,
        permutations = 999
    )
    ```

    You can repeat this for the other distance matrices by replacing `bray`.

    Questions to consider:

    - Does one group appear more dispersed than another?
    - Is the pattern of dispersion the same for every beta-diversity metric?
    - How might differences in dispersion affect your interpretation of PERMANOVA?


---

# 9. Putting everything together

You have now looked at the microbiome from several different perspectives:

**Sequencing depth and rarefaction**

> Did we sequence the samples deeply enough to make meaningful diversity comparisons?

**Taxonomic composition**

> Which bacterial groups were detected, and what patterns can we see across mice?

**Alpha diversity**

> How diverse is the microbial community within each mouse?

**Beta diversity**

> How different are the complete microbial communities between mice?

**Statistical analysis**

> Which of the patterns we observed are supported statistically?


!!! question "Exercise 13 — Final interpretation"

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

# 10. Important limitations

Before interpreting the experiment, consider some of its limitations.

!!! warning "Things to keep in mind"

    **Sample size**

    There are only six mice per group. How might this affect the conclusions that can be drawn?

    **Relative abundance**

    16S sequencing primarily provides information about the relative composition of the microbial community. What information about absolute bacterial abundance is missing?

    **Rarefaction**

    Rarefaction standardizes sequencing depth, but it also discards reads. What are the consequences of that trade-off?

    **Taxonomic resolution**

    Short-read 16S sequencing may not reliably distinguish microorganisms at species or strain level. How does this limit biological interpretation?

    **Taxonomy is not function**

    Knowing which bacterial taxa are present does not automatically tell us what those microorganisms are doing.

    **Association is not causation**

    Think carefully about which conclusions are supported by this experiment and which would require additional experiments.


---

# 11. Going further with your figures

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
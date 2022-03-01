‘Downstream’ Analysis: from counts to Differential Expression
================
Nirupama Tamvada

**Reminder:** When answering the questions, it is not sufficient to
provide only code. Add explanatory text to interpret your results
throughout.

### Background

The dataset used for this assignment has been published by [Li et al. in
2016](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005511),
as well as included in a follow-up study by [Houston-Ludlam et al. also
in
2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0159197).
You should recognize the follow-up study, “Comparative Transcriptome
Profiling of Human Foreskin Fibroblasts Infected with the Sylvio and Y
Strains of Trypanosoma cruzi” as this is the paper you were asked to
review for the Paper Critique. The data used in this follow-up study was
a combination of newly generated data (for the Sylvio strain and one
replicate of Y strain), as well as the data generated in the original Li
et al. article (all for the Y strain). In this assignment, you’ll be
working with the Y strain data generated in the original article.

Although the processed data for this study are not available in GEO, the
*raw* RNA-Seq reads have been submitted to NCBI SRA under the accession
[SRP043008](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP043008).
Fortunately, they have also been preprocessed using a standardized
pipeline from raw reads to mapped gene and transcript level counts, and
formatted into `RangedSummarizedExperiment` objects (ready for analysis
in R!) as part of the [Recount2
project](https://jhubiostatistics.shinyapps.io/recount/). The
`RangedSummarizedExperiment` experiment class functions much the same as
the `SummarizedExperiment` class we learned about in lecture. We can
read the data into R in this format using the [recount
package](https://bioconductor.org/packages/release/bioc/html/recount.html).

We can also get metadata (mapping of sample IDs to time points,
infection status, batch, etc) from [Table
S1](https://journals.plos.org/plospathogens/article/file?type=supplementary&id=info:doi/10.1371/journal.ppat.1005511.s011)
in the original study.

Some helper code to read the data in is provided.

### Load libraries

Note: only a few libraries are provided for the initial steps; add
additional libraries as needed for the rest of the analyses.

``` r
library(recount)
library(openxlsx)
library(tidyverse)
library(edgeR)
library(pheatmap)
library(knitr)
```

### Question 1: Importing the data and getting familiar with it (3 POINTS)

-   First, we’ll download and read the gene-level counts
    `RangedSummarizedExperiment` object for project id “SRP043008” into
    R using the `recount` package. We’ll add a column called `sample` to
    the `colData` slot that contains the `run` variable (since this will
    serve as our sample ID and is easier). And we’ll also add the sample
    IDs to the colnames of the object. Note that there’s nothing you
    need to add here.

``` r
download_study(project="SRP043008")
```

    ## 2022-02-28 13:09:02 downloading file rse_gene.Rdata to SRP043008

``` r
load(file.path("SRP043008", "rse_gene.Rdata"))
colData(rse_gene)$sample <- colData(rse_gene)$run
colnames(rse_gene) <- colData(rse_gene)$sample
```

-   Now we’ll add sample IDs to the `colnames` of the object to the SRR
    sample IDs: Use the IDs in the `run` variable of the `colData` slot.
    Note that there’s nothing you need to add here.

``` r
colnames(rse_gene) <- colData(rse_gene)$sample
```

-   Investigate the object just loaded. How many genes are there?

``` r
length(unique(rowData(rse_gene)$gene_id))
```

    ## [1] 58037

There are 58037 genes

-   How many samples are there?

``` r
nrow(colData(rse_gene))
```

    ## [1] 27

There are 27 samples

-   Here, we convert the `RangedSummarizedExperiment` object into a
    [`DGEList`](https://rdrr.io/bioc/edgeR/man/DGEList.html) object (for
    use with `limma` and/or `edgeR`). As we discussed in lecture, this
    format is similar to `SummarizedExperiment`, but has some special
    slots, including one for normalization factors, that are important
    for methods like `limma` and `edgeR`. We’ll include the `colData` of
    the `RangedSummarizedExperiment` object as the `samples` data.frame
    in the `DGEList`. Note there’s nothing you need to add here.

``` r
dge <- DGEList(counts = assays(rse_gene)$counts,
               samples = colData(rse_gene))
```

-   Next, we’ll read in the metadata from [Table
    S1](https://journals.plos.org/plospathogens/article/file?type=supplementary&id=info:doi/10.1371/journal.ppat.1005511.s011)
    in the original study. Note that this is an Excel file with some
    nonstandard formatting, so this step needs to be done with care -
    the `read.xlsx` package is suggested since this can skip some rows,
    and fill in merged cells. Note there’s nothing you need to add here.

``` r
mdat <- read.xlsx("https://journals.plos.org/plospathogens/article/file?type=supplementary&id=info:doi/10.1371/journal.ppat.1005511.s011",
                  startRow = 3, fillMergedCells = TRUE) %>%
  mutate(sample=Accession.Number)
```

-   Subset the metadata to include only those samples that we have gene
    counts for (Table S1 also includes the non-host (non-human)
    samples), and add this information to the `samples` slot of the
    `DGEList` object (the `samples` slot of the `DGEList` is where the
    sample metadata is stored - analogous to `colData()` for a
    `SummarizedExperiment`). HINT: perform a join operation by sample id
    (beginning with “SRR”).

``` r
# Subsetting to include samples that we have gene counts for 
dge$samples <- dge$samples %>%
  left_join(mdat, by = "sample")
```

-   How many variables of interest are in our experimental design? Are
    these factors (categorical) or continuous? If factors, how many
    levels are there? List out the levels for each factor.

There are 4 variables of interest: developmental stage (hours post
infection or hpi), Infection status, and Batch (although this is is a
blocking factor). Although hpi can be considered as a factor (with 6
levels: 4, 6, 12, 24, 48, 72), we treat it as a continuous variable for
the analysis. Infection status is a categorical factor with two levels
and Batch is a factor with five levels.

``` r
# Levels of infection status
table(dge$samples$Infected)
```

    ## 
    ##  N  Y 
    ##  9 18

``` r
# Levels of batch 
table(dge$samples$Batch)
```

    ## 
    ##  A  B  C  D  E 
    ## 10  4  5  6  2

### Question 2: Remove lowly expressed genes (2 POINTS)

-   Remove lowly expressed genes by retaining genes that have CPM &gt; 1
    in at least 25% of samples.

``` r
# Create copy of dge
dge2 <- dge

# Convert to counts-per-million 
dge2$counts <- cpm(dge2$counts, log = FALSE, normalized.lib.sizes = FALSE)

# Keeping genes that have CPM > 1 in at least 25% of samples
keep <- which((rowSums(dge2$counts > 1)/27) >= 0.25)

# Keeping gene subset
dge <- dge[keep, ]
dge2 <- dge2[keep, ]
```

-   How many genes are there after filtering?

``` r
# Number of genes after filtering
length(keep)
```

    ## [1] 16681

There are 16681 genes after filtering

### Question 3: Data manipulation (2 POINTS)

The different character values of `Developmental.stage` refer to time
points - these can be thought of as categorical or on a continous axis.
In order to make graphing easier, it will be helpful to convert this
variable to a numeric representation.

-   Create a new column in the samples metadata tibble. Call it “hpi”
    (which stands for hours post infection) and populate it with the
    appropriate numeric values.

``` r
# Creating hpi column
dge$samples$hpi <- gsub("[^0-9]", "", dge$samples$Developmental.stage)

# Converting to numeric 
dge$samples$hpi <- as.numeric(dge$samples$hpi)
```

### Question 4: Assessing overall distributions (4 POINTS)

-   The expression values are raw counts. Calculate TMM normalization
    factors (and add them to your `DGEList` object).

``` r
# Calculate TMM normalization counts
dge <- calcNormFactors(dge, method = "TMM")
```

-   Examine the distribution of gene expression on the scale of
    $\\sf{log\_{2}}$ CPM across all samples using box plots (with
    samples on the x-axis and expression on the y-axis). Hint 1: Add a
    small pseudo count (e.g. 1) before taking the log2 transformation to
    avoid taking log of zero. Hint 2: To get the data in a format
    amenable to plotting, some data manipulation steps are required.
    Take a look at the `pivot_longer` function in the tidyr package.

``` r
# Convert counts matrix to dataframe for plotting 
dge2$counts <-  dge2$counts %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id")

# Convert counts dataframe to long format 
dge2$counts <-
  as.data.frame(dge2$counts) %>%
  pivot_longer(cols =  2:ncol(dge2$counts), values_to = "cpm",
                 names_to = "sample")


# Computing log2(CPM +1) counts
dge2$counts$log2cpm <- log2(dge2$counts$cpm + 1)

# Boxplots of gene expression distribution across samples 
ggplot(dge2$counts, aes(sample, log2cpm)) + geom_boxplot(aes(color = sample)) +
  labs(x = "Sample", color = "Sample", y = "Log2(CPM+1) gene expression counts") +
  theme_bw()
```

<img src="analysis_assignment_files/figure-gfm/boxplot_expression-1.png" width="50%" style="display: block; margin: auto;" />

-   Examine the distribution of gene expression in units of
    $\\sf{log\_{2}}$ CPM across all samples using overlapping density
    plots (with expression on the x-axis and density on the y-axis; with
    one line per sample and lines coloured by sample).

``` r
# Density plot of gene expression coloured by samples 
 ggplot(dge2$counts, aes(sample)) + geom_density(aes(x = log2cpm, color = sample)) +
  labs(x = "Log2(CPM+1) gene expression counts", color = "Sample", y = "Density") +
  theme_bw()
```

<img src="analysis_assignment_files/figure-gfm/density_expression-1.png" width="50%" style="display: block; margin: auto;" />

-   Which sample stands out as different, in terms of the distribution
    of expression values, compared to the rest?

The distribution of SRR1346028 stands out as the most different with a
relatively lower density around the range of 3-8 Log2(CPM + 1) gene
expression counts compared to the other samples.

### Question 5: Single gene graphing (3 POINTS)

-   Find the expression profile for the gene *OAS* (Ensembl ID
    ENSG00000089127). Make a scatterplot with hours post-infection on
    the x-axis and expression value in log2(CPM + 1) on the y-axis.
    Color the data points by infection status, and add in a regression
    line for each one.

``` r
# Filter for ENSG00000089127 gene expression data
gene_match <- grepl("ENSG00000089127", dge2$counts$gene_id)

# Subset original dataset
oas_subset <- dge2$counts[gene_match, ]

# Add in hpi information 
oas_subset <- oas_subset %>%
  left_join(dge$samples[, c("sample", "hpi", "Infected")], by = "sample")

# Scatterplot
ggplot(oas_subset, aes(hpi, log2cpm, color = Infected)) + geom_point(aes(color = Infected)) + labs(x = "Hours post infection", y = "Log2(CPM + 1)", color = "Infection status") +
  theme_bw() + geom_smooth(method="lm")
```

    ## `geom_smooth()` using formula 'y ~ x'

<img src="analysis_assignment_files/figure-gfm/hpi_scplot-1.png" width="50%" style="display: block; margin: auto;" />

-   Is there sign of interaction between infection status and hours post
    infection for **OAS**? Explain using what you observed in your graph
    from the previous question.

Yes, there appears to be an interaction between infection status and
hours post infection. In infected cells, we see an increasing trend with
time post infection, while there is no changing trend in uninfected
cells.

### Question 6: How do the samples correlate with one another? (4 POINTS)

-   Examine the correlation **between samples** using one or more
    heatmaps (i.e. samples should be on the x axis and the y axis, and
    the values in the heatmap should be correlations). Again, use the
    log2 transformed CPM values. Display batch, hpi, and infection
    status for each sample in the heatmap. Hint: Consider using
    pheatmap() with annotations and cor to correlate gene expression
    between each pair of samples.

``` r
# Correlation matrix
cor_samples <- cor(dge$counts)

# Metadata for heatmap 
cor_mdat <- dge$samples[, c("Infected", "Batch", "hpi")]
row.names(cor_mdat) <- colnames(cor_samples)

# Heatmap
pheatmap(cor_samples, annotation_col = cor_mdat, cluster_rows = TRUE)
```

<img src="analysis_assignment_files/figure-gfm/heatmap_corr-1.png" width="50%" style="display: block; margin: auto;" />

-   Among the variables batch, hpi, and infection status, which one
    seems to be most strongly correlated with clusters in gene
    expression data? Hint: Consider using ‘cluster\_rows=TRUE’ in
    pheatmap().

Overall, experimental variables seem to have the strongest correlation
with the clusters. The variable infection status seems to be the most
strongly correlated with clusters, as all the clusters seem to be
composed largely of infected samples. Quite a few of the clusters seem
to contain samples from Batch A as well, and so batch is also likelu
correlated with clusters to a lesser extent.

-   There is a sample whose expression values do not correlate as highly
    with other samples of the sample hpi, and in general. Identify this
    sample by its ID.

This sample is SRR1346028 (the same sample identified in the density
distribution as well).

### Question 7: Construct linear model for Differential expression analysis (4 POINTS)

-   First set up a model matrix with hpi, infection status, and the
    interaction between hpi and infection status as covariates. Then
    calculate variance weights with voom, and generate the mean-variance
    trend plot. Hint: use the `DGEList` object as input to voom, since
    you want to input raw counts here (voom will automatically transform
    to log2 CPM internally).

``` r
# Setting up model matrix
modm <- model.matrix(~ hpi*Infected, data = dge$samples)

head(modm) %>% kable
```

| (Intercept) | hpi | InfectedY | hpi:InfectedY |
|------------:|----:|----------:|--------------:|
|           1 |   4 |         0 |             0 |
|           1 |   4 |         1 |             4 |
|           1 |   4 |         1 |             4 |
|           1 |   4 |         1 |             4 |
|           1 |   6 |         0 |             0 |
|           1 |   6 |         1 |             6 |

``` r
# Calculating voom variance weights 
vw <- voom(dge, design = modm, plot = TRUE, span = 0.5)  
```

<img src="analysis_assignment_files/figure-gfm/mv_plot-1.png" width="50%" style="display: block; margin: auto;" />

-   Use limma to fit the linear model with the model matrix you just
    created.

``` r
# Fitting weighted limma
lmvoom <- lmFit(vw, modm)
# Empirical Bayes moderation of standard errors
lmvoom <- eBayes(lmvoom)
```

-   Print the 10 top-ranked genes by adjusted p-value for the
    hpi:Infected coefficient

``` r
# Top 10 ranked genes
signif(topTable(lmvoom, number = 10, coef = "hpi:InfectedY", sort.by = "p"), 3) %>%
  kable()
```

|                    |  logFC | AveExpr |    t | P.Value | adj.P.Val |    B |
|:-------------------|-------:|--------:|-----:|--------:|----------:|-----:|
| ENSG00000137285.9  | 0.0512 |   6.030 | 7.30 | 1.0e-07 |   0.00131 | 7.14 |
| ENSG00000137267.5  | 0.0421 |   6.590 | 6.90 | 2.0e-07 |   0.00178 | 6.09 |
| ENSG00000139531.12 | 0.0180 |   4.090 | 6.70 | 4.0e-07 |   0.00199 | 5.71 |
| ENSG00000103966.10 | 0.0326 |   6.400 | 6.48 | 6.0e-07 |   0.00264 | 4.99 |
| ENSG00000281028.1  | 0.0322 |   4.600 | 6.22 | 1.2e-06 |   0.00409 | 4.45 |
| ENSG00000106785.14 | 0.0471 |   6.450 | 5.92 | 2.7e-06 |   0.00546 | 3.53 |
| ENSG00000002549.12 | 0.0534 |   6.610 | 5.92 | 2.7e-06 |   0.00546 | 3.52 |
| ENSG00000147509.13 | 0.0459 |   0.803 | 5.92 | 2.7e-06 |   0.00546 | 4.19 |
| ENSG00000272921.1  | 0.0180 |   4.750 | 5.84 | 3.3e-06 |   0.00546 | 3.42 |
| ENSG00000038210.12 | 0.0308 |   5.590 | 5.81 | 3.6e-06 |   0.00546 | 3.26 |

### Question 8: Interpret model (2 POINTS)

-   For the gene CDC20 (Ensembl id ENSG00000117399), what is the numeric
    value of the coeffcient of the hpi term? Interpret this value for
    infected and non-infected samples?

``` r
# Fixing genes column name
coef_table <- topTable(lmvoom, number = Inf)  %>%
  rownames_to_column(var = "gene_id")
  
# Extracting hpi coefficient for CDC20
coef_table  %>%
   filter(grepl("ENSG00000117399", gene_id)) %>%
  kable()
```

| gene\_id           |        hpi | InfectedY | hpi.InfectedY |  AveExpr |        F | P.Value | adj.P.Val |
|:-------------------|-----------:|----------:|--------------:|---------:|---------:|--------:|----------:|
| ENSG00000117399.13 | -0.0686189 | 0.0851386 |     0.0162618 | 5.624199 | 60.16853 |       0 |         0 |

The numeric value for the hpi term coefficient is -0.0686. For
non-infected samples, with a unit increase in hours per infection, the
expression of the CDC20 gene (in log2-fold-changes) changes by -0.0697.
For infected samples, with a unit increase in hours per infection, the
expression of the CDC20 gene (in log2-fold-changes) changes by
(0.0162618-0.0696) -0.0519545.

### Question 9: Quantify the number of genes differentially expressed (3 POINTS)

-   Using the linear model defined above, determine the number of genes
    differentially expressed by infection status *at any time point* at
    an FDR (use adjust.method = “fdr” in `topTable`) less than 0.05.
    Hint: in other words, test the null hypothesis that all coefficients
    involving infection status are equal to zero.

``` r
topTable(lmvoom, number = Inf,
coef = c("InfectedY", "hpi:InfectedY"), adjust.method = "fdr", p.value = 0.05) %>% nrow()
```

    ## [1] 2006

The number of genes differentially expressed by infection status are
2006.

### Question 10: Interpret the interaction term (3 POINTS)

-   Explain what you are modeling with the interaction term. For a
    particular gene, what does a significant interaction term mean?

For this model, the interaction term represents the partial dependency
of the infection status on the effect of hours per infection. This means
that the effect/slope of the hours post infection is different for the
two infection statuses. A significant interaction term for a particular
gene would thus mean that the effect of hours post infection on the
expression of that gene id different between an infected and uninfected
sample.

#### **Bonus Question** (2 POINTS)

-   Compare your DE results to those obtained by [Li et
    al. (2016)](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005511).
    Discuss any discrepancies. List at least three explanations for
    these discrepancies.

Many possible explanations for discrepancies, including: - The paper
reports DE genes for each time point separately, and here we consider
hpi as continuous.  
- The paper adjusts for uninfected (mock infection) by subtracting
counts (Infected-Non infected) instead of including infection status as
a covariate. - The paper removed a sample as an outlier.  
- The paper also included a filter for logFC when defining
significance.  
- The paper included batch as a covariate in the linear model.  
- The paper used q-values instead of fdr/BH.

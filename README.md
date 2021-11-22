DJExpress README
================
Lina Marcela Gallego-Paez
22/10/2021

-   [1 Installation](#installation)
-   [2 Details](#details)
-   [3 Example DJE analysis](#example-dje-analysis)
    -   [3.1 Import junction quantification files with
        DJEimport()](#import-junction-quantification-files-with-djeimport)
    -   [3.2 Annotate junctions with their respective gene of origin
        with
        DJEannotate()](#annotate-junctions-with-their-respective-gene-of-origin-with-djeannotate)
    -   [3.3 Filtering junctions for differential expression analysis
        using
        DJEprepare()](#filtering-junctions-for-differential-expression-analysis-using-djeprepare)
    -   [3.4 Test for Differential Junction Expression using
        DJEanalyze()](#test-for-differential-junction-expression-using-DJEanalyze)
    -   [3.5 Gene-wise Splice plots](#gene-wise-splice-plots)
    -   [3.6 Association between junction expression and external traits
        with
        DJEvsTrait()](#association-between-junction-expression-and-external-traits-with-djevstrait)
-   [4 Example Junction Co-expression Network Analysis
    (JCNA)](#example-junction-co-expression-network-analysis-jcna)
    -   [4.1 Prepare data for JCNA
        1-pass](#prepare-data-for-jcna-1-pass)
    -   [4.2 1-pass JCNA](#1-pass-jcna)
    -   [4.3 2-pass JCNA](#2-pass-jcna)

# 1 Installation

The devtools package provides install\_github() that allows you to
install the package from GitHub **(this option will be available once
the repository becomes public)**:

``` r
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("MauerLab/DJExpress")
```

DJExpress can be also installed from the source (DJExpress\_0.1.0.tar.gz
file available with granted access at
<https://gitlab.com/MauerLab/djexpress-source-file>) file using
install.packages():

``` r
install.packages(path_to_DJExpress_sourcefile, repos = NULL, type="source")
```

# 2 Details

The DJExpress package provides several convenient functions for
Differential Junction Expression analysis (DJE) as well as Junction
Co-expression Network Analysis (JCNA):

-   DJEimport()

-   DJEannotate()

-   DJEprepare()

-   DJEanalyze()

-   DJEvsTrait()

-   DJEplotSplice()

-   DJEspliceRadar()

-   JCNAprepare()

-   JCNA1pass()

-   JCNAgenePrepare()

-   JCNA2pass()

-   JCNAModTrait()

# 3 Example DJE analysis

DJExpress uses two types of files as the primary input for DJE analysis:

1.  STAR aligner-derived “SJ.out.tab” files containing splice junction
    counts per sample (generated from FASTQ or BAM files) in
    tab-delimited format, or outputs from any other junction
    quantification tool as long as they contain junction IDs as first
    columns, following the format chr: start: end: strand
    (e.g. chr1:123:456:1, where positive or negative strand are coded as
    1 and 2, respectively)

2.  The correspondent transcriptome annotation (gtf file) used for
    junction quantification.

Our recommended STAR alignment settings are the following:

**STAR –genomeDir GENOME –readFilesIn READ1 READ2 –runThreadN 4
–outFilterMultimapScoreRange 1 –outFilterMultimapNmax 20
–outFilterMismatchNmax 10 –alignIntronMax 500000 –alignMatesGapMax
1000000 –sjdbScore 2 –alignSJDBoverhangMin 1 –genomeLoad NoSharedMemory
–limitBAMsortRAM 70000000000 –readFilesCommand cat
–outFilterMatchNminOverLread 0.33 –outFilterScoreMinOverLread 0.33
–sjdbOverhang 100 –outSAMstrandField intronMotif –outSAMattributes NH HI
NM MD AS XS –sjdbGTFfile GENCODE\_ANNOTATION –limitSjdbInsertNsj 2000000
–outSAMunmapped None –outSAMtype BAM SortedByCoordinate –outSAMheaderHD
@HD VN:1.4 –outSAMattrRGline ID::<ID> –twopassMode Basic –outSAMmultNmax
1**

## 3.1 Import junction quantification files with DJEimport()

The first step in the DJE analysis is to merge “SJ.out.tab” files into a
single expression matrix (in case the primary input contains individual
files per sample) and to generate a matrix of junction coordinates for
junction-to-gene mapping. The input for DJEimport() is the path to the
folder where these junction quantification files are located.

As an example, we are going to use a few “SJ.out.tab” files produced by
STAR alignment using 3 TCGA colorectal (COADREAD) tumor samples and 3
matching normal colon samples from GTEx:

``` r
library(DJExpress)
in.file <- system.file("extdata", "junct.quant", package = "DJExpress")
print(in.file)
```

    ## [1] "/Users/paez/Library/R/4.0/library/DJExpress/extdata/junct.quant"

You can see that **in.file** contains the path to the folder
(“junct.quant”) where “SJ.out.tab” files from STAR alignment are found.

This is how individual “SJ.out.tab” files look like:

``` r
files <- list.files(in.file)
junctions <- readr::read_tsv(paste0(in.file,"/", files[1]), comment = "", col_names = FALSE)
head(junctions)
```

    ## # A tibble: 6 x 9
    ##   X1       X2    X3    X4    X5    X6    X7    X8    X9
    ##   <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 chr1  14830 14929     2     2     1     2     0    29
    ## 2 chr1  14830 14969     2     2     1    19   275    38
    ## 3 chr1  14830 15020     2     2     1     0    12    18
    ## 4 chr1  15013 25232     1     1     1     0     1    13
    ## 5 chr1  15039 15795     2     2     1     4   115    37
    ## 6 chr1  15039 16853     2     2     1     0    10     1

Now, let’s generate the junction quantification matrix and the
respective junction coordinates matrix. We are going to determine that 3
reads per junction is the minimal number of reads for a junction to be
considered expressed in the sample, which means that every read count
below 3 will be transformed to 0 (this is the default value for
**min.expressed** in **DJEimport()**):

``` r
out.file <- DJEimport(workDir = in.file, aligner="STAR", min.expressed = 3)
# summary of out.file:
summary(out.file)
```

    ##       Length Class      Mode
    ## quant 6      data.frame list
    ## coord 5      data.frame list

``` r
# Head of out.file$quant:
head(out.file$quant)
```

    ##                    GTEx1.aligned.bam.junctions GTEx2.aligned.bam.junctions
    ## chr1:14830:14929:2                           0                           0
    ## chr1:14830:14969:2                          19                          71
    ## chr1:14830:15020:2                           0                           0
    ## chr1:15013:25232:1                           0                           0
    ## chr1:15039:15795:2                           4                          30
    ## chr1:15039:16853:2                           0                           0
    ##                    GTEx3.aligned.bam.junctions TCGA1.aligned.junctions
    ## chr1:14830:14929:2                           0                       0
    ## chr1:14830:14969:2                           9                       0
    ## chr1:14830:15020:2                           0                       0
    ## chr1:15013:25232:1                           0                       0
    ## chr1:15039:15795:2                           0                       9
    ## chr1:15039:16853:2                           0                       0
    ##                    TCGA2.aligned.junctions TCGA3.aligned.junctions
    ## chr1:14830:14929:2                       0                       0
    ## chr1:14830:14969:2                      13                       0
    ## chr1:14830:15020:2                       0                       0
    ## chr1:15013:25232:1                       0                       0
    ## chr1:15039:15795:2                      15                       0
    ## chr1:15039:16853:2                       0                       0

``` r
# Head of out.file$coord:
head(out.file$coord)
```

    ##    chr start   end strand        junction_id
    ## 1 chr1 14830 14929      - chr1:14830:14929:2
    ## 2 chr1 14830 14969      - chr1:14830:14969:2
    ## 3 chr1 14830 15020      - chr1:14830:15020:2
    ## 4 chr1 15013 25232      + chr1:15013:25232:1
    ## 5 chr1 15039 15795      - chr1:15039:15795:2
    ## 6 chr1 15039 16853      - chr1:15039:16853:2

## 3.2 Annotate junctions with their respective gene of origin with DJEannotate()

We now need to associate junctions found in **out.file** to their
respective gene, so gene-level differential junction usage can be
calculated downstream in the analysis.

For this step we need the path to the transcriptome annotation file (or
gtf file) and the output object of DJEimport(). We are going to use an
example gtf file containing only chromosome 1 (chr1) for simplicity:

``` r
# Get path to gtf file:
gtf0 <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
print(gtf0)
```

    ## [1] "/Users/paez/Library/R/4.0/library/DJExpress/extdata/chr1.gtf.gz"

``` r
# Run DJEannotate():
ann.out <- DJEannotate(import.out = out.file, gtf = gtf0)

# Summary of ann.out:
summary(ann.out)
```

    ##                 Length Class      Mode     
    ## quant.annotated   6    data.frame list     
    ## featureID       420    -none-     character
    ## groupID         420    -none-     character

``` r
# head of ann.out$quant.annotated (the junction quantification matrix)
head(ann.out$quant.annotated)
```

    ##                    GTEx1.aligned.bam.junctions GTEx2.aligned.bam.junctions
    ## chr1:14830:14929:2                           0                           0
    ## chr1:14830:14969:2                          19                          71
    ## chr1:14830:15020:2                           0                           0
    ## chr1:15039:15795:2                           4                          30
    ## chr1:15039:16853:2                           0                           0
    ## chr1:15943:16606:2                           0                           3
    ##                    GTEx3.aligned.bam.junctions TCGA1.aligned.junctions
    ## chr1:14830:14929:2                           0                       0
    ## chr1:14830:14969:2                           9                       0
    ## chr1:14830:15020:2                           0                       0
    ## chr1:15039:15795:2                           0                       9
    ## chr1:15039:16853:2                           0                       0
    ## chr1:15943:16606:2                           0                       0
    ##                    TCGA2.aligned.junctions TCGA3.aligned.junctions
    ## chr1:14830:14929:2                       0                       0
    ## chr1:14830:14969:2                      13                       0
    ## chr1:14830:15020:2                       0                       0
    ## chr1:15039:15795:2                      15                       0
    ## chr1:15039:16853:2                       0                       0
    ## chr1:15943:16606:2                       6                       0

``` r
# head of ann.out$featureID (the junction ID list)
head(ann.out$featureID)
```

    ## [1] "chr1:14830:14929:2" "chr1:14830:14969:2" "chr1:14830:15020:2"
    ## [4] "chr1:15039:15795:2" "chr1:15039:16853:2" "chr1:15943:16606:2"

``` r
# head of ann.out$groupID (the junction-associated gene list)
head(ann.out$groupID)
```

    ## [1] "WASH7P" "WASH7P" "WASH7P" "WASH7P" "WASH7P" "WASH7P"

## 3.3 Filtering junctions for differential expression analysis using DJEprepare()

Once we have junctions annotated with their respective gene of origin
within **ann.out**, we can now filter them based on user-defined
expression cutoffs.

For this step, we need to provide the sample names (as defined by the
colnames in ann.out$quant.annotated) that correspond to the basal
experimental condition (e.g. WT/normal condition, non-treated condition,
etc). In the case of this example, the GTEx normal colon samples
correspond to our basal experimental condition.

Additionally, we are going to use the default expression cutoffs for
DJEprepare (minimum of read count mean per junction - **minMean** of 10
reads and minimum of read count variance per junction - **minVar** of
0):

``` r
Group1 <- colnames(ann.out$quant.annotated)[grep("GTEx", colnames(ann.out$quant.annotated))]
print(Group1) # GTEx sample names
```

    ## [1] "GTEx1.aligned.bam.junctions" "GTEx2.aligned.bam.junctions"
    ## [3] "GTEx3.aligned.bam.junctions"

``` r
# Run DJEprepare
prep.out <- DJEprepare(annotate.out = ann.out, Group1 = Group1,
                       minMean = 10,maxMean = Inf,minVar = 0,maxVar = Inf) # these are the default values for minMean, maxMean, minVar and MaxVar

# Summary of prep.out
summary(prep.out)
```

    ##               Length Class      Mode     
    ## JunctExprfilt  6     data.frame list     
    ## featureID     38     -none-     character
    ## groupID       38     -none-     character
    ## design        12     -none-     numeric

DJEprepare() contains expression-based filtered junction counts, gene
annotation and design matrix for differential expression analysis

## 3.4 Test for Differential Junction Expression using DJEanalyze()

At this point we have all the data we need for the differential junction
expression analysis.

An FDR cutoff of 0.05 and a minimal absolute log-fold change (logFC) of
2 are default paramenters for DJEanalyze(). Users are free to explore
and modify the multiple options for the function (type ?DJEanalyze for
more details):

``` r
# Run DJEanalyze
anlz.out <- DJEanalyze(prepare.out = prep.out, Group1 = Group1,
                       FDR = 0.05, logFC = 2)
```

    ## Total number of junctions:  38 
    ## Total number of genes:  4 
    ## Number of genes with 1 junction:  1 
    ## Mean number of junctions in a gene:  10 
    ## Max number of junctions in a gene:  19

With DJEanalyze, a tests for differential junction expression is
implemented using voom’s log2-counts per million (logCPM) and
observation-level weights based on mean-variance relationship (stored at
**anlz.out$v.norm**).

A linear model is then fit per junction using a provided experimental
design, and empirical Bayes moderated t-statistics are implemented to
assess the significance level of the observed expression changes.

The linear model framework of limma is also used in parallel to
calculate differential junction usage, where significant differences in
log-fold-changes in the fit model between junctions from the same gene
are tested (using the diffSplice function from limma and stored at
**anlz.out$ex.norm**).

DJExpress thereby identifies alternatively spliced regions in
transcripts based on two main features of splice junction expression: 1)
Quantitative changes in the abundance of individual junctions between
experimental groups, and 2) Differences in their expression levels
compared to the average expression of other junctions in the gene.

Following these criteria, splice junctions are classified based on their
absolute log-fold change (e.g. experimental condition A vs B) and their
relative log-fold change (target junction vs all other junctions in the
gene) in one of the following expression groups:

-   **Group 0**: Junctions without differential expression or
    differential usage.
-   **Group 1**: Junctions with equal levels of differential expression
    and differential usage, reflecting changes in splicing patterns
    between experimental conditions (in this case, both absolute and
    relative log-fold change values are similar, if not the same).
-   **Group 2**: Junctions with differential expression but no
    differential usage or vice-versa, implying the occurrence of
    generalized changes in expression across the gene, rather than the
    presence of a differentially spliced region (in this case, either
    the absolute or relative log-fold change value is not significant).
-   **Group 3**: Junctions with divergent levels of differential
    expression and differential usage, indicating concomitant changes in
    splicing and total gene expression (in this case, the absolute and
    relative log-fold change values can substantially vary from each
    other).

Let’s check the summary of **anlz.out** object:

``` r
summary(anlz.out)
```

    ##              Length Class      Mode
    ## v.norm        4     EList      list
    ## ex.norm      17     MArrayLM   list
    ## dje.out      25     data.frame list
    ## dje.sig       0     -none-     NULL
    ## logFC.plot    0     -none-     NULL
    ## volcano.plot  0     -none-     NULL
    ## model.fit     0     -none-     NULL
    ## group.par     0     -none-     NULL

-   **v.norm** correspond to the log-cpm values outputed by limma:voom
-   **ex.norm** contains the differential junction expression analysis
    output
-   **dje.out** correspond to the annotated **ex.norm** data set with
    additional information, including basic statistics (e.g. median,
    zero counts, etc) and DJE group for each junction.
-   **dje.sig** a subset of **dje.out** containing significant hits
    based on FDR and logFC cutoffs and the expression group of each
    junction.
-   **logFC.plot** Shows the regression plot of Absolute logFC \~
    Relative logFC for differentially expressed junctions (basal vs
    tested sample group)
-   **volcano.plot** Volcano plot of differential junction expression
-   **model.fit** Correspond to the confidence and prediction intervals
    for the linear regression (Absolute logFC \~ Relative logFC)
    required to define the expression group of each junction.
-   **group.par** Indicates the “Group” type for each junction (Group
    0,1,2 or 3) based on its position within the “relative vs absolute
    logFC” regression plot (**logFC.plot**) whose parameters are defined
    by the fit model (**model.fit**), and user-defined FDR and logFC
    cutoffs.

For most users, the top relevant element in DJEanalyze() output object
is **dje.sig** , a data frame with the differentially expressed
junctions considered significant based on the FDR and logFC cutoffs
*(for this particular example, since the gtf used only contained chr1
-and consequently, only genes in this chromosome- there were no
differentially expressed junctions detected and thus **dje.sig**,
**logFC.plot**, **volcano.plot**, **model.fit** and **group.par** are
not produced).*

The summary statistics and additional information about all junctions
analyzed by DJExpress() can be found in **dje.out**:

``` r
head(anlz.out$dje.out)
```

    ##             junctionID     GeneID     logFC          t    P.Value       FDR
    ## 1 chr1:569491:569920:2 AL669831.3 -7.399331 -2.4981874 0.01357649 0.4060932
    ## 2 chr1:569199:569362:2 AL669831.3 -6.920520 -2.3156498 0.02195099 0.4060932
    ## 3 chr1:568990:569195:2 AL669831.3  3.883564  1.4650514 0.14502801 0.9320353
    ## 4 chr1:569158:569317:2 AL669831.3 -3.658810 -1.1936361 0.23453015 0.9320353
    ## 5   chr1:15948:16606:2     WASH7P  3.528810  1.0169887 0.31081845 0.9320353
    ## 6 chr1:569297:569340:2 AL669831.3  2.318959  0.7760885 0.43893519 0.9320353
    ##   medianExp.group1   meanExp.group1 medianExp.group2   meanExp.group2
    ## 1                0 204.333333333333                0                0
    ## 2                0 150.333333333333                0                0
    ## 3                0                0               71 63.6666666666667
    ## 4                0               94                0                0
    ## 5                0               19                4 6.33333333333333
    ## 6                0                0                0 394.666666666667
    ##   medianNormExp.group1 meanNormExp.group1 medianNormExp.group2
    ## 1     10.8601062067676   13.3887641130092     7.44121906358489
    ## 2     10.8601062067676   13.2413183276152     7.44121906358489
    ## 3     10.2985733721812   9.96845760178822     14.6417753676264
    ## 4     10.8601062067676   13.0158266208891     7.44121906358489
    ## 5     10.2985733721812   12.2502876187697     10.6518290322903
    ## 6     10.2985733721812   9.96845760178822     7.44121906358489
    ##   meanNormExp.group2 group1WithCounts group2WithCounts zeroCounts.group1
    ## 1   7.44830555314337                1                0                 2
    ## 2   7.44830555314337                1                0                 2
    ## 3   12.4725591108128                0                2                 3
    ## 4   7.44830555314337                1                0                 2
    ## 5   10.1563459904198                1                2                 2
    ## 6   11.1849930475568                0                1                 3
    ##   zeroCounts.group2 neojunction logFC.ebayes AveExpr.ebayes   t.ebayes
    ## 1                 3        <NA>    -7.790190       10.41853 -2.4486829
    ## 2                 3        <NA>    -7.363267       10.34481 -2.2917382
    ## 3                 1 neojunction     2.447053       11.22051  0.8708416
    ## 4                 3        <NA>    -4.388091       10.23207 -1.3283461
    ## 5                 1        <NA>    -1.675349       11.20332 -0.4402782
    ## 6                 2 neojunction     1.110836       10.57673  0.3458094
    ##   P.Value.ebayes adj.P.Val.ebayes  B.ebayes
    ## 1     0.01547513        0.4425836 -4.592654
    ## 2     0.02329387        0.4425836 -4.593063
    ## 3     0.38521396        0.7921862 -4.595273
    ## 4     0.18605478        0.7921862 -4.594770
    ## 5     0.66036110        0.9934724 -4.595398
    ## 6     0.72996414        0.9934724 -4.595546

The column **neojunction** specifies the junctions that are uniquely
found in the tested condition (TCGA tumor samples in this example) and
never found in any sample from the basal condition.

When a junction annotation file is provided as an additional input to
DJEanalyze() (), **dje.out** contains an additional column called
**annotation**, where junctions not found in the gtf transcriptome
annotation file are indicated as ***unnanotated***.

***Note**: A junction annotation file can be produced using the JAn()
function in DJExpress. users need to provide the path to the gtf file
and JAn() generates a data frame with all exon-exon junctions found in
the gft and their correspondent gene ID (this can take some time
depending on the size of the gtf file).*

``` r
gtf <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
Jan.out <- JAn(gtf)

head(Jan.out)
```

    ##                gene           junction
    ## 1 ENSG00000223972.5 chr1:12228:12612:1
    ## 2 ENSG00000223972.5 chr1:12722:13220:1
    ## 3 ENSG00000223972.5 chr1:12058:12178:1
    ## 5 ENSG00000223972.5 chr1:12698:12974:1
    ## 6 ENSG00000223972.5 chr1:13053:13220:1
    ## 7 ENSG00000223972.5 chr1:13375:13452:1

## 3.5 Gene-wise Splice plots

One of the main features of DJExpress is the generation of interactive
gene-wise junction plots using **DJEplotSplice()** function. This type
of representation facilitates straight-forward visual inspection of
differential splicing across a gene of interest and exploration of
supplementary information about each junction’s expression, including
the above-mentioned classification based on absolute and relative
log-fold change patterns and basic statistics on expression levels
(e.g. mean and median expression in each experimental condition, number
of samples expressing the junction, etc.).

As example, we show here the gene-wise splice plot for the ENAH gene
indicating the presence of an exon inclusion event in tumor samples:

``` r
# load DJEanalyze ouptut object (30 TCGA COADREAD tumor vs 30 GTEx colon tissue samples):
data(DJEanlz)

# generate gene-wise splice plot:
iPlot.out <- DJEplotSplice(DJEanlz, geneID="ENAH", logFC = 0.5, FDR = 0.05)
iPlot.out$plot
```

<img src="ReadFig/tutorial_plotSplice.png" width="787" />

Up- and down-regulted Junctions (with both relative and absolute logFC
values above the specified threshold) are shown in red and blue,
respectively. The user can explore further statistical information about
each individual junction by hoovering over the junction dots in the
plot, which displays a box with summarized DJE information, including
relative and absolute logFC values, FDR values and expression group of
the selected junction.

When the path to a gtf file is provided, iPlot.out object contains an
additional gene model plot with exon-to-protein domain annotation and
the localization of a user-selected junction (e.g. the exon skipping
junction downregulated in ENAH):

``` r
# Indicate path to gtf file:
gtf0 <- "~/Downloads/gencode.v32lift37.annotation.gtf" # gtf downloaded from GENCODE website

# generate gene-wise splice plot:
plot <-DJEplotSplice(DJEanlz, geneID="ENAH",logFC = 0.5, FDR = 0.05,
                     gtf = gtf0,
                     target.junction = "chr1:225688773:225695652:2")

# Show gene model plot:
plot$JunctionToGene
```

<img src="ReadFig/DJEplotSplice_ENAH.png" width="1192" />

Colors within exonic regions in the gene model plot indicate the
presence of protein domains and/or post-translational modifications
(PTMs). The position of the selected junction within the gene model plot
is indicated by a dashed arc whose color correspond to the type of
differential expression (blue for downregulation and red for
upregulation).

This visualization strategy facilitates straight-forward visual
inspection of differential splicing across the entire gene (not
restricted to individual local splicing events) and exploration of
supplementary information about each junction’s expression, including
the relative position of target alternative splicing events to protein
domains and PTMs, which further facilitates the undestanding of the
functional consequences of such alterations in transcript and protein
structure.

## 3.6 Association between junction expression and external traits with DJEvsTrait()

In order to explore the physiologic significance of target alternative
splicing events, the user can make associations between junction
expression and external sample traits (e.g. clinical data, mutation
data, gene expression, etc) when available.

**DJEvsTrait()** function recieves as input:

1.  a **DJEanalyze** output object (analyze.out), and
2.  a numeric vector or a matrix of external sample traits (traitData).

User should also indicate the experimental condition that should be
excluded from the junction-trait association test (Group1, we want to
keep only tumor sample data).

**DJEvsTrait()** offers two type of association test: association based
on correlation, or association based on regression (test.type).

In the case of correlation analysis, **DJEvsTrait()** executes a big
matrix correlation test (“bicor”, “pearson”, “kendall” or “spearman”)
between the normalized junction expression contained within the
**DJEanalyze** output object and values in the trait data, for the
selected samples in **Group1**. The return object contains the matrices
with correlation coefficients and associated P-values for each
junction-trait pair.

When regression analysis is selected, **DJEvsTrait()** uses large matrix
operations adapted from the ***Matrix eQTL*** algorithms
(<http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/>) and
returns a data frame with significant junction-trait associations based
on the selected model in **useModel** (modelLINEAR, modelANOVA or
modelLINEAR\_CROSS).

Users have the additional option to supply a specific set of target
junctions (maximum 5 junction IDs) to generate SpliceRadar plots (using
**DJEspliceRadar** function) when test.type = “Correlation”.

For our example, we are going to use the expression of genes known to
code for splicing factors as external trait data from the TCGA COADREAD
samples. Our aim is to test potential associations between the
expression levels of splicing factors and the expression of splice
junctions, which could contribute to the characterization of regulatory
networks controlling alternative splicing events of interest.

Our example trait data looks like this:

``` r
SF <- system.file("extdata", "SF.expr.rds", package = "DJExpress")
SF.exp <- readRDS(SF)
print(SF.exp[,c(1:8)], row.names = FALSE)
```

    ##     ACIN1    AGGF1      AQR   ARGLU1     BAG2      BCAS1    BCAS2     BUB3
    ##  7.552864 5.587345 6.066815 7.076861 2.074548  4.3489730 4.494495 6.478300
    ##  7.053783 5.416257 4.330139 8.234791 2.230377  5.1375484 4.106623 7.122091
    ##  7.545562 5.097509 5.660424 7.551051 3.282097  6.4118113 4.770974 7.290347
    ##  7.381274 4.710511 5.348603 6.641707 4.510515  3.2940607 4.182552 7.234991
    ##  7.446583 4.644171 5.170452 6.859362 5.155595 -0.5479771 4.940280 6.876452
    ##  7.146317 5.147694 5.094728 7.849169 2.477678  2.4382041 4.782368 6.909769
    ##  6.538969 4.753975 4.666928 7.444170 4.952424  4.2047463 4.530543 7.103473
    ##  7.076956 4.982824 5.109519 6.120249 2.863316  6.3985415 4.808807 6.746150
    ##  7.189313 3.854566 5.113996 5.519869 3.421598  3.0286828 4.355646 7.424761
    ##  6.615288 4.500616 5.127311 5.957938 4.170996  4.6175814 5.412729 6.713595
    ##  6.069404 4.920051 4.398614 7.865325 3.211644  2.9211268 4.805560 7.008063
    ##  6.960877 4.398029 5.583253 6.787906 4.387801  4.1299934 4.850540 7.359703
    ##  6.909462 5.144576 4.477908 8.053648 2.483248  4.5048661 5.109640 6.885118
    ##  7.192344 5.130073 5.836868 6.075499 3.354691  2.9441545 4.937601 6.990735
    ##  7.058618 5.185318 5.511891 6.104409 2.781132  2.7857434 5.246682 6.954692
    ##  6.832429 5.094747 5.171595 7.211786 3.015575  4.8738161 5.095292 7.087915
    ##  7.221428 5.873249 6.148053 5.410874 4.351467  6.0306128 5.804583 7.679160
    ##  6.549775 5.818928 5.873622 7.986173 3.100379  3.1409853 4.492918 7.225574
    ##  7.405144 4.974626 5.657289 5.800262 4.283825  5.1377954 4.779288 6.719201
    ##  7.415361 5.615890 5.915267 6.943296 4.325100  6.3536787 5.757642 7.747063
    ##  6.848851 5.780174 6.119275 6.766806 3.688394  4.4958346 5.594975 7.479008
    ##  7.339277 5.358377 5.716494 5.815694 2.734738  4.9442019 5.291155 6.269628
    ##  7.234535 5.336751 5.890506 6.888185 1.336695  6.4180489 4.747839 6.467778
    ##  7.370953 5.070476 6.701108 6.839450 2.886882  6.3873475 5.161956 6.949431
    ##  7.045549 5.423561 5.702272 5.312090 4.628792  5.7806371 5.199481 7.175501
    ##  7.426605 5.137016 5.372259 7.260949 2.809443  7.5712496 4.694615 6.334887
    ##  7.062514 5.068249 4.827781 7.293300 3.837951  3.2426898 4.179654 6.812228
    ##  7.662508 4.393593 5.395414 7.319873 4.155953  4.3963637 4.778637 6.805694
    ##  6.910246 5.365189 5.115793 7.432051 1.900800  6.0665772 5.161969 6.603271
    ##  6.828736 5.364454 5.668844 6.356400 5.324406  3.7972869 5.874960 7.502759

Junction-trait association analysis based on biweight mid-correlation
(bicor) shows these results:

``` r
# define sample group for junction-trait association:
Group1 <- colnames(DJEanlz$v.norm$E)[grep("SRR", colnames(DJEanlz$v.norm$E))]

# Run DJEvsTrait:
DT.out <- DJEvsTrait(analyze.out = DJEanlz, Group1 = Group1,traitData = SF.exp,
                     coeff = 0.2,select.junctions = c("chr1:225692756:225695652:2",
                                                      "chr1:225688773:225692692:2",
                                                      "chr1:225688773:225695652:2"),
                     test.type = "Correlation", cor.method = "bicor")
```

    ## Allowing parallel execution with up to 2 working processes.
    ## [1] "subsetting associations to: chr1:225692756:225695652:2"
    ## [2] "subsetting associations to: chr1:225688773:225692692:2"
    ## [3] "subsetting associations to: chr1:225688773:225695652:2"

``` r
# Summary of DT.out:
summary(DT.out)
```

    ##             Length Class  Mode   
    ## TraitCor    385920 -none- numeric
    ## TraitPvalue 385920 -none- numeric
    ## sig.cor          3 -none- list

``` r
# Correlation coefficients:
head(DT.out$TraitCor[,c(1:10)])
```

    ##                                  ACIN1        AGGF1         AQR      ARGLU1
    ## chr16:67863992:67864292:2   0.33027264  0.161770864  0.14655085 -0.19782666
    ## chr17:73512698:73512826:1  -0.19339096  0.177985264  0.03910401  0.29665809
    ## chr3:150321287:150340171:1 -0.12109937 -0.015321447  0.11433838  0.15094432
    ## chr19:54970588:54970643:1  -0.03313364 -0.008815586 -0.36675411  0.22514066
    ## chr7:134853813:134855150:2  0.09934680  0.153334718  0.39719648  0.01736267
    ## chr19:54969701:54971945:1   0.13933309  0.401910104  0.18723459  0.15628675
    ##                                   BAG2        BCAS1       BCAS2        BUB3
    ## chr16:67863992:67864292:2   0.01183096 -0.013778808  0.08835670 -0.09855143
    ## chr17:73512698:73512826:1  -0.30689089 -0.056371682 -0.13733771  0.19555748
    ## chr3:150321287:150340171:1 -0.03725409 -0.008522488  0.03627020 -0.00957963
    ## chr19:54970588:54970643:1  -0.03049230 -0.101366952 -0.09463717  0.21727652
    ## chr7:134853813:134855150:2 -0.07880536 -0.056826845 -0.13440598  0.25771843
    ## chr19:54969701:54971945:1  -0.31797439  0.163953042 -0.15865312 -0.19529143
    ##                                   BUD13       BUD31
    ## chr16:67863992:67864292:2   0.310379904  0.14446102
    ## chr17:73512698:73512826:1   0.030813246  0.29721521
    ## chr3:150321287:150340171:1 -0.308716116 -0.11785015
    ## chr19:54970588:54970643:1  -0.398680758  0.12816278
    ## chr7:134853813:134855150:2 -0.009644815  0.08440710
    ## chr19:54969701:54971945:1   0.358594965  0.07900475

``` r
# Correlation p-values:
head(DT.out$TraitPvalue[,c(1:10)])
```

    ##                                 ACIN1      AGGF1        AQR    ARGLU1
    ## chr16:67863992:67864292:2  0.07466582 0.39307645 0.43965805 0.2946843
    ## chr17:73512698:73512826:1  0.30586162 0.34671266 0.83744804 0.1114087
    ## chr3:150321287:150340171:1 0.52382712 0.93595276 0.54742388 0.4259137
    ## chr19:54970588:54970643:1  0.86200967 0.96312371 0.04620474 0.2316263
    ## chr7:134853813:134855150:2 0.60144862 0.41853646 0.02975046 0.9274410
    ## chr19:54969701:54971945:1  0.46275102 0.02769423 0.32180880 0.4095248
    ##                                  BAG2     BCAS1     BCAS2      BUB3      BUD13
    ## chr16:67863992:67864292:2  0.95052351 0.9423902 0.6424359 0.6043772 0.09505788
    ## chr17:73512698:73512826:1  0.09903153 0.7673248 0.4692460 0.3003696 0.87159104
    ## chr3:150321287:150340171:1 0.84504360 0.9643490 0.8490889 0.9599301 0.09693744
    ## chr19:54970588:54970643:1  0.87291772 0.5940376 0.6188770 0.2487676 0.02909002
    ## chr7:134853813:134855150:2 0.67891899 0.7654980 0.4788746 0.1691444 0.95965766
    ## chr19:54969701:54971945:1  0.08682700 0.3866383 0.4023805 0.3010406 0.05166453
    ##                                BUD31
    ## chr16:67863992:67864292:2  0.4462792
    ## chr17:73512698:73512826:1  0.1107065
    ## chr3:150321287:150340171:1 0.5351055
    ## chr19:54970588:54970643:1  0.4997138
    ## chr7:134853813:134855150:2 0.6574294
    ## chr19:54969701:54971945:1  0.6781497

We have selected the 3 junction IDs involved in the exon inclusion event
in ENAH gene to generate **sig.cor** table with the subset of
significant associations to splicing factor expression. This table will
be use to generate the SpliceRadar plot using **DJEspliceRadar()**
function:

``` r
# sig.cor table:
summary(DT.out$sig.cor)

# define sample group for junction-trait association:
Sr.out <- DJEspliceRadar(DT.out, ordered.junction = "chr1:225688773:225692692:2")

Sr.out
```

<img src="ReadFig/tutorial_spliceradar.png" width="1123" />

In the SpliceRadar plot, the coefficient of top-ranked correlations
between the three ENAH junctions is used to map each junction-trait
association within a radar chart. Positive correlation coefficients are
located within the outer region and negative correlation coefficients
are found within the inner region of the radar chart.

In this example, we can see that both upregulated exon inclusion
junctions (in red and dark red) correlate in a consistent way with an
specific subset of splicing factors expression, and such association
pattern tends to be inversed for the downregulated exon exclusion
junction (in blue).

SpliceRadar plot concept thus allows for the simultaneous visual
inspection of relevant associations between the expression of selected
junctions and external traits (e.g. splicing factor expression), as well
as for the elucidatation of expression-trait patterns shared among
junctions of interest with potential biological relevance.

# 4 Example Junction Co-expression Network Analysis (JCNA)

The weighted junction co-expression network analysis module (JCNA) in
DJExpress provides an implementation of WGCNA algorithms (version
1.70.3, Langfelder & Horvath, 2008) in the context of splice junction
expression.

Before testing the JCNA module with custom data, we highly recommend
users to get familiarized to co-expression network analysis concepts and
pipelines by following WGCNA tutorials
(<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/>).

## 4.1 Prepare data for JCNA 1-pass

The initial inputs for JCNA module are:

1.  A **DJEanalyze** output object from where junction quantification
    and additional sample information is retrieved.

-   Alternatively, users can directly provide junction read counts
    table, such as the one produced by STAR alignment (the choice should
    be indicated in *input.type = c(“DJEanalyze.out”,
    “junction.counts”)* and the path to folder where individual junction
    quantification files are located should be indicated in the
    *workDir* argument. This folder should only contain junction
    quantification files).

-   As suggested by WGCNA guidelines
    (<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html>),
    sufficient sample size should be provided ( &gt;= 15 samples within
    single experimental conditions) in order to avoid the constructon of
    noisy expression networks that do not trully represent biologically
    meaningful associations.

2.  A vector or factor specifying sample names corresponding to the
    experimental condition that should be excluded from the JCNA
    analysis (Group1, we want to keep only tumor sample data). This
    should be indicated only when *input.type = “DJEanalyze.out”*.

3.  A numeric vector or a matrix of external sample traits. Samples
    should be as rows and traits as columns.

Additional options can be further explored by typing **?JCNAprepare**.

The first step in the JCNA module is to prepare data for a first round
of co-expression network construction. For this, outlier samples
(identified based on sample dendrogram clustering) as well as missing
entries, entries with weights below a threshold, and zero-variance
junctions are removed from the analysis.

After this, the soft-thresholding power for the network construction
should be chosen. For this, JCNAprepare contains a wrapper of the \*\_\_
pickSoftThreshold()\_\_\* function from WGCNA, which performs the
analysis of network topology, using some predifined set of
soft-thresholding testing powers as suggested in WGCNA guidelines. This
step generates two plots: 1) the scale-free fit index vs
soft-thresholding power plot, and 2) The mean connectivity vs
soft-thresholding power plot. This representation is helpful for the
selection of a proper soft-thresholding power that will be used later
during JCNA 1-pass step.

For this example, we are going to use a full version of the
**DJEanalyze()** output we used before. This object has to be loaded
from GitLab:

``` r
# Load DJEanalyze output for JCNA:
githubURL <- "https://gitlab.com/MauerLab/djexpress-DJEanalyze-output-file/raw/main/DJEanlz.total.RData"
load(url(githubURL))
DJEanlz <- analyze.coadread

# Summary of DJEanlz:
summary(DJEanlz)
```

    ##              Length Class        Mode
    ## v.norm        4     EList        list
    ## ex.norm      17     MArrayLM     list
    ## dje.out      26     data.frame   list
    ## dje.sig      26     data.frame   list
    ## logFC.plot    3     recordedplot list
    ## volcano.plot  9     gg           list
    ## model.fit     4     data.frame   list
    ## group.par    14     -none-       list

``` r
# Load splicing factor expression as trait data:
SF <- system.file("extdata", "SF.expr.rds", package = "DJExpress")
SF.exp <- readRDS(SF)

# Change format of colnames for visualization:
colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))] <- paste0("TCGA_",
seq(1,length(colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))]), 1))

# Set Group1 (normal tissue) to exclude from the analysis
Group1 <- colnames(DJEanlz$v.norm$E)[grep("SRR", colnames(DJEanlz$v.norm$E))]

# Run JCNAprepare (it takes some minutes):
Jprep <- JCNAprepare(analyze.out=DJEanlz, Group1 = Group1,
traitData = SF.exp, abline.threshold=600, input.type = "DJEanalyze.out")
```

    ## Allowing parallel execution with up to 2 working processes.
    ##  Flagging junctions and samples with too many missing values...
    ##   ..step 1

    ## pickSoftThreshold: will use block size 390.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 390 of 114639
    ##    ..working on genes 391 through 780 of 114639
    ##    ...     
    ##    ..working on genes 114271 through 114639 of 114639
    ##    Power SFT.R.sq   slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1 0.000349  0.0636          0.883 25400.00  24400.00  41600
    ## 2      2 0.359000 -1.2500          0.900  8670.00   7730.00  21200
    ## 3      3 0.687000 -1.7300          0.937  3680.00   2970.00  12700
    ## 4      4 0.790000 -1.9200          0.948  1790.00   1290.00   8270
    ## 5      5 0.819000 -2.0100          0.949   964.00    612.00   5750
    ## 6      6 0.836000 -2.0300          0.957   558.00    313.00   4170
    ## 7      7 0.846000 -2.0300          0.963   343.00    168.00   3120
    ## 8      8 0.843000 -2.0400          0.965   220.00     94.70   2400
    ## 9      9 0.835000 -2.0600          0.965   147.00     55.70   1880
    ## 10    10 0.839000 -2.0500          0.971   102.00     34.10   1500
    ## 11    12 0.840000 -2.0400          0.977    52.80     14.20    990
    ## 12    14 0.828000 -2.0400          0.976    29.90      6.83    680
    ## 13    16 0.820000 -2.0100          0.963    18.20      3.66    482
    ## 14    18 0.798000 -1.9200          0.906    11.80      2.13    350
    ## 15    20 0.944000 -1.6000          0.936     8.09      1.32    259

``` r
# Summary of JCNAprepare output:
summary(Jprep)
```

    ##            Length Class        Mode
    ## sample.den      3 recordedplot list
    ## datExpr    114639 data.frame   list
    ## datTraits     360 data.frame   list
    ## NetTop          3 recordedplot list
    ## sft             2 -none-       list

We had defined an abline.threshold of 60, to allow the removal of the
outlier sample observed in the sample dendrogram:

<img src="ReadFig/tutorial_sampleclust.png" width="571" />

**Jprep$NetTop** shows the 2 plots generated during the analysis of
network topology. The lowest power for which the scale-free topology fit
index curve flattens out upon reaching a high value (here around 0.90)
is 18. This is the power selected for the first round of network
construction (**JCNA1pass()** function).

``` r
# Network topology analysis plots:
Jprep$NetTop
```

<img src="ReadFig/tutorial_nettop.png" width="571" />

## 4.2 1-pass JCNA

Once data is pre-processed, **JCNA1pass()** function (which is a wrapper
of the **blockwiseModules()** function in WGCNA) will construct the
junction network and identify modules of expression. For this,
correlation matrices (e.g. using Pearson, Spearman or the default
biweight midcorrelation) are built for all pair-wise junctions. The full
network is specified by a weighted adjacency matrix calculated using the
soft threshold power determined by **JCNAprepare()**.

Users should explore the multiple parameters that can be adjusted for
**JCNA1pass()** (e.g. minimum module size, module detection sensitivity,
cut height of the hierarchical clustering dendrogram for module
definition, etc) in order to fine-tuning network construction. We have
kept default parameter values also suggested by WGCNA guidelines to work
well in a variety of settings.

``` r
# Run JCNA1pass (this will take some minutes. You can try more threads specified in nThreads, in case they are available in your system):
J1pass <- JCNA1pass(Jprep, cor.method = "bicor", nThreads = 2)
```

    ## [1] "Constructing the junction network and modules"
    ## Allowing parallel execution with up to 2 working processes.
    ##  Calculating module eigengenes block-wise from all genes
    ##    Flagging genes and samples with too many missing values...
    ##     ..step 1
    ##  ....pre-clustering genes to determine blocks..
    ##    Projective K-means:
    ##    ..k-means clustering..
    ##    ..merging smaller clusters...
    ## Block sizes:
    ## gBlocks
    ##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
    ## 4999 4999 4998 4997 4993 4991 4990 4988 4985 4969 4967 4964 4960 4942 4903 4896 
    ##   17   18   19   20   21   22   23   24 
    ## 4862 4851 4838 4710 4663 4167 4043 2964 
    ##  ..Working on block 1 .
    ##     TOM calculation: adjacency..
    ##     ..will use 2 parallel threads.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 1 into file JCNA_blockwiseTOM-block.1.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##      ..removing 6 genes from module 2 because their KME is too low.
    ##  ..Working on block 2 .
    ##     TOM calculation: adjacency..
    ##     ..will use 2 parallel threads.
    ##      Fraction of slow calculations: 0.000000
    ##     ..connectivity..
    ##     ..matrix multiplication (system BLAS)..
    ##     ..normalization..
    ##     ..done.
    ##    ..saving TOM for block 2 into file JCNA_blockwiseTOM-block.2.RData
    ##  ....clustering..
    ##  ....detecting modules..
    ##  ....calculating module eigengenes..
    ##  ....checking kME in modules..
    ##  ...      
    ##      ..removing 123 genes from module 1 because their KME is too low.
    ##      ..removing 104 genes from module 2 because their KME is too low.
    ##  ..merging modules that are too close..
    ##      mergeCloseModules: Merging modules whose distance is less than 0.25
    ##        Calculating new MEs...
    ## [1] " modules identified:32"
    ## 
    ##     0     1     2     3     4     5     6     7     8     9    10    11    12 
    ## 59055 24919  8625  2498  1564  1423  1356  1352  1262  1180  1174  1036   987 
    ##    13    14    15    16    17    18    19    20    21    22    23    24    25 
    ##   958   743   719   677   625   621   517   410   338   309   302   293   278 
    ##    26    27    28    29    30    31 
    ##   268   255   244   229   218   204 
    ## [1] "Label 0 (grey color) is reserved for genes outside of all modules"

    ## [1] "Calculating module membership values and trait correlations"

    ## [1] "Defining junction significance for each trait"
    ## [1] "done"

``` r
# Summary of JCNA1pass output:
summary(J1pass)
```

    ##                      Length Class        Mode     
    ## net                      10 -none-       list     
    ## datExpr              114639 data.frame   list     
    ## datTraits               360 data.frame   list     
    ## moduleColors         114639 -none-       character
    ## module.den               24 -none-       list     
    ## moduleTraitCor        11520 -none-       numeric  
    ## moduleTraitPvalue     11520 -none-       numeric  
    ## juncModuleMembership     32 data.frame   list     
    ## juncMMPvalue             32 data.frame   list     
    ## ModuleTrait               3 recordedplot list     
    ## Junctrait               720 data.frame   list     
    ## sig.cors.trait          360 -none-       list

``` r
# Junction module dendrogram (first 4 dendrogram blocks shown):
J1pass$module.den[[1]]
J1pass$module.den[[2]]
J1pass$module.den[[3]]
J1pass$module.den[[4]]
```

<img src="ReadFig/tutorial_moduleden1.png" width="571" /><img src="ReadFig/tutorial_moduleden2.png" width="571" /><img src="ReadFig/tutorial_moduleden3.png" width="571" /><img src="ReadFig/tutorial_moduleden4.png" width="571" />

``` r
# Assignment of junctions to respective module:
J1pass$net$colors[1:10] # only first 10 displayed
```

    ##   chr1:14830:14969:2   chr1:17369:17605:2   chr1:17743:17914:2 
    ##                 "23"                  "0"                  "0" 
    ## chr1:569158:569317:2 chr1:763156:764382:1 chr1:878439:878632:1 
    ##                 "18"                  "0"                  "4" 
    ## chr1:878758:879077:1 chr1:879189:879287:1 chr1:880181:880421:2 
    ##                  "4"                  "4"                  "0" 
    ## chr1:880181:880436:2 
    ##                  "0"

``` r
# correspondent module color:
WGCNA::labels2colors(J1pass$net$colors)[1:10]
```

    ##  [1] "grey60"      "turquoise"   "turquoise"   "greenyellow" "turquoise"  
    ##  [6] "white"       "white"       "white"       "turquoise"   "turquoise"

``` r
# Module Eigengenes per sample:
J1pass$net$MEs[,c(1:10)] # only first 10 modules displayed
```

    ##                 ME20         ME16        ME15         ME14          ME19
    ## TCGA_1   0.103025422 -0.085979768  0.07872600  0.036201625  0.1155179449
    ## TCGA_3   0.077388088  0.208277206  0.28396200  0.075029346  0.1339219076
    ## TCGA_4   0.283525364 -0.015779129  0.24892425  0.240541474  0.1659999945
    ## TCGA_5   0.004878349  0.256686160  0.29543648  0.165513127  0.1793970536
    ## TCGA_6   0.196332433 -0.202724060  0.06576807  0.219980665  0.1606143575
    ## TCGA_7  -0.309061016 -0.146727796 -0.30270714 -0.518139586 -0.5193446470
    ## TCGA_8  -0.180901898 -0.083101393  0.12028732  0.014987299 -0.0009810184
    ## TCGA_9  -0.014158443 -0.156754417 -0.30826127 -0.339741467 -0.4613021783
    ## TCGA_10  0.471977840  0.249560020 -0.16770460  0.064724664 -0.1646598532
    ## TCGA_11 -0.086727123 -0.110125043 -0.23698113 -0.384640759 -0.3233076214
    ## TCGA_12  0.130058683  0.040692105 -0.13100943  0.123838866  0.1587360867
    ## TCGA_13 -0.177380898 -0.375057512  0.06719161 -0.383660548 -0.3372767844
    ## TCGA_14  0.074352071 -0.070708327  0.25641499  0.053958173  0.1744532425
    ## TCGA_15  0.166798830 -0.117336995 -0.09967809  0.155346101  0.0646724857
    ## TCGA_16  0.124458900 -0.003326440  0.25037229  0.145776218 -0.0056271056
    ## TCGA_17 -0.177184502  0.372111905 -0.13085418  0.013332421  0.0622467349
    ## TCGA_18 -0.010139358  0.018419138  0.14658309  0.033578241 -0.0502907902
    ## TCGA_19 -0.105347359  0.100940621 -0.18893045 -0.089993877  0.0279343516
    ## TCGA_20  0.063090995  0.370859516 -0.10886590  0.138195004  0.1251715822
    ## TCGA_21  0.247349442 -0.177080300 -0.03528816  0.110462422 -0.0121397858
    ## TCGA_22 -0.094210008 -0.098454899  0.26669427  0.140130354  0.1560703358
    ## TCGA_23 -0.134667800 -0.072609008  0.17590567 -0.003969142  0.0559605664
    ## TCGA_24 -0.237916877  0.005787073  0.02423824 -0.010766956  0.1576507066
    ## TCGA_25 -0.145077975  0.272782176 -0.14672240 -0.069686683  0.0925133214
    ## TCGA_26 -0.329115346  0.158809785 -0.14792298 -0.126770944 -0.0029033820
    ## TCGA_27 -0.136027994 -0.177009294 -0.17842783 -0.057076877 -0.0169949330
    ## TCGA_28  0.004267095 -0.248977965 -0.20414304  0.057011195  0.0806757604
    ## TCGA_29 -0.064642694 -0.065824790  0.08029586  0.046692531  0.0394803395
    ## TCGA_30  0.255055780  0.152651429  0.02669649  0.149147113 -0.0561886726
    ##                ME21         ME28         ME24           ME8          ME9
    ## TCGA_1   0.02011512  0.070213520  0.070963544  0.1005742246  0.175879038
    ## TCGA_3   0.14533509  0.240513389  0.097109630  0.0413671717  0.138486295
    ## TCGA_4   0.11171969  0.055470603 -0.045384832  0.0008061462  0.149309270
    ## TCGA_5   0.11419077 -0.046389587 -0.071892523 -0.0333749086  0.124228831
    ## TCGA_6   0.09573965  0.087567816  0.176708955  0.0502289070  0.203514489
    ## TCGA_7  -0.56731235 -0.514299457 -0.509578918 -0.5065169274 -0.514144174
    ## TCGA_8   0.04682959  0.176652437  0.122147771  0.1492983690  0.132499407
    ## TCGA_9  -0.40115667 -0.472828249 -0.475456370 -0.4795429081 -0.502983071
    ## TCGA_10  0.03066148 -0.080665162  0.187103585 -0.1929783169 -0.075595198
    ## TCGA_11 -0.26099680 -0.221260234 -0.222182747 -0.2732343291 -0.262071880
    ## TCGA_12  0.24218073 -0.001809426  0.178738966  0.1210553752  0.113177819
    ## TCGA_13 -0.41560879 -0.361331695 -0.398709962 -0.3412629796 -0.303770508
    ## TCGA_14  0.16489960  0.031372683  0.207118923  0.1264930861  0.165070491
    ## TCGA_15 -0.01011077 -0.029484605  0.166182866  0.0733883529  0.064718440
    ## TCGA_16  0.01513443  0.032763378 -0.002736111  0.0047487797 -0.009480872
    ## TCGA_17 -0.02383060  0.109283526 -0.064485131  0.0399616752 -0.085063472
    ## TCGA_18  0.07927667  0.059840559  0.099922532  0.0697602976  0.011851582
    ## TCGA_19 -0.04270465 -0.004187904 -0.001851495  0.1159348878  0.014957875
    ## TCGA_20  0.11938908  0.142667343  0.009276971 -0.0185853821 -0.035615425
    ## TCGA_21  0.01726009 -0.060542289 -0.105072928 -0.0929620723 -0.152614392
    ## TCGA_22  0.18527830  0.084729323  0.198922796  0.1649739534  0.115474758
    ## TCGA_23  0.08295565  0.192813193  0.102252084  0.1844101536  0.117218602
    ## TCGA_24  0.08978370  0.021838165  0.063634696  0.1573683600  0.062455533
    ## TCGA_25  0.11098264  0.221709377  0.079938572  0.1175522919  0.030759915
    ## TCGA_26 -0.13839989  0.190518526 -0.023554770  0.0783302282  0.070103393
    ## TCGA_27 -0.05795956 -0.013003750  0.035733110  0.1301201803  0.101919311
    ## TCGA_28  0.03837291  0.034808037  0.046159493  0.1750222699  0.134001677
    ## TCGA_29  0.08933600  0.173935633  0.101934375  0.1373147219  0.124860247
    ## TCGA_30  0.11863890 -0.120895148 -0.022943083 -0.1002516083 -0.109147980

``` r
# Junction module membership values:
J1pass$juncModuleMembership[c(1:10),c(1:10)] # only first 10 rows and columns displayed
```

    ##                             MM20         MM16          MM15         MM14
    ## chr1:14830:14969:2   -0.14877223  0.052944771 -0.0198539391 -0.153526309
    ## chr1:17369:17605:2   -0.07993977 -0.248403998 -0.1247984181 -0.039483325
    ## chr1:17743:17914:2   -0.06385915 -0.185576948 -0.1324204839  0.049227883
    ## chr1:569158:569317:2 -0.11662347 -0.023285794  0.0007658819 -0.138827147
    ## chr1:763156:764382:1 -0.11415083  0.072718117  0.1349391358  0.317648078
    ## chr1:878439:878632:1 -0.08693247 -0.163699774  0.0773272097 -0.007002024
    ## chr1:878758:879077:1 -0.03800737 -0.189611988  0.1140912037 -0.042286393
    ## chr1:879189:879287:1  0.08201547 -0.177359464  0.1636950030  0.191416081
    ## chr1:880181:880421:2  0.06301739  0.096558628  0.4203080058  0.409683730
    ## chr1:880181:880436:2 -0.04915496  0.005130104  0.3382846892  0.375097592
    ##                             MM19          MM0        MM21         MM28
    ## chr1:14830:14969:2   -0.04000515  0.346179640 -0.14458737  0.041006498
    ## chr1:17369:17605:2   -0.01793684 -0.130997017 -0.30675197  0.183418347
    ## chr1:17743:17914:2   -0.04453809 -0.257697292 -0.16682654  0.225098647
    ## chr1:569158:569317:2 -0.02478559  0.423370092 -0.03723287 -0.003926982
    ## chr1:763156:764382:1  0.50333991  0.008196651  0.31873508  0.182771097
    ## chr1:878439:878632:1  0.10963058 -0.133608529  0.01334184  0.123534518
    ## chr1:878758:879077:1  0.06231297 -0.100335353  0.02817631  0.046707983
    ## chr1:879189:879287:1  0.17892070 -0.097559013  0.04887074  0.129816343
    ## chr1:880181:880421:2  0.62576224  0.368134688  0.34923326  0.531837366
    ## chr1:880181:880436:2  0.65508869  0.158788228  0.23439304  0.519129483
    ##                             MM24        MM8
    ## chr1:14830:14969:2   -0.33596049 -0.2822884
    ## chr1:17369:17605:2    0.02556507  0.1961402
    ## chr1:17743:17914:2    0.23872819  0.2891339
    ## chr1:569158:569317:2 -0.29098799 -0.2397877
    ## chr1:763156:764382:1  0.32636260  0.3269458
    ## chr1:878439:878632:1  0.38115907  0.5768248
    ## chr1:878758:879077:1  0.42411405  0.5074444
    ## chr1:879189:879287:1  0.45894553  0.5130710
    ## chr1:880181:880421:2  0.26577456  0.3033289
    ## chr1:880181:880436:2  0.24762612  0.4620172

``` r
# Junction-to-trait significance (GS, correlation coefficient) and associated P-value (p.GS):
J1pass$Junctrait[c(1:10),c(1:10)]
```

    ##                          GS.ACIN1   p.GS.ACIN1    GS.AGGF1 p.GS.AGGF1
    ## chr1:14830:14969:2    0.003528657 0.9855059914 -0.02646876 0.89159018
    ## chr1:17369:17605:2    0.161718535 0.4019685800 -0.15925833 0.40926172
    ## chr1:17743:17914:2   -0.014423959 0.9408016343  0.04374008 0.82175054
    ## chr1:569158:569317:2 -0.024231807 0.9007053778 -0.06793467 0.72622462
    ## chr1:763156:764382:1  0.201223382 0.2952333639 -0.13403070 0.48820542
    ## chr1:878439:878632:1  0.224990026 0.2406223560 -0.31449073 0.09660253
    ## chr1:878758:879077:1  0.129715207 0.5024414455 -0.33822422 0.07271972
    ## chr1:879189:879287:1  0.226587708 0.2372101839 -0.18279777 0.34255685
    ## chr1:880181:880421:2  0.543446928 0.0023134152  0.25834020 0.17602036
    ## chr1:880181:880436:2  0.609123130 0.0004533334  0.23732484 0.21511629
    ##                           GS.AQR    p.GS.AQR   GS.ARGLU1 p.GS.ARGLU1
    ## chr1:14830:14969:2   -0.33093484 0.079509429  0.19997543  0.29830121
    ## chr1:17369:17605:2   -0.10506178 0.587547789  0.38183554  0.04095541
    ## chr1:17743:17914:2    0.06747463 0.728008859  0.22980615  0.23043486
    ## chr1:569158:569317:2 -0.35154750 0.061475357  0.15507589  0.42182918
    ## chr1:763156:764382:1  0.10487526 0.588214866 -0.05413338  0.78032269
    ## chr1:878439:878632:1 -0.02684557 0.890056038 -0.29338520  0.12243145
    ## chr1:878758:879077:1 -0.03644917 0.851101775  0.01034415  0.95752797
    ## chr1:879189:879287:1  0.04844585 0.802928058 -0.17567310  0.36201492
    ## chr1:880181:880421:2  0.35379534 0.059720121 -0.03337284  0.86354720
    ## chr1:880181:880436:2  0.47477966 0.009255848 -0.22068007  0.24998911
    ##                           GS.BAG2  p.GS.BAG2
    ## chr1:14830:14969:2   -0.003491630 0.98565806
    ## chr1:17369:17605:2   -0.385933921 0.03866051
    ## chr1:17743:17914:2   -0.341306410 0.06998722
    ## chr1:569158:569317:2  0.006052153 0.97514320
    ## chr1:763156:764382:1  0.169537762 0.37928292
    ## chr1:878439:878632:1 -0.172823404 0.36997683
    ## chr1:878758:879077:1 -0.313905330 0.09725837
    ## chr1:879189:879287:1 -0.282704082 0.13729338
    ## chr1:880181:880421:2 -0.170888366 0.37544122
    ## chr1:880181:880436:2 -0.231039250 0.22787373

``` r
# Module-trait significant associations found:
J1pass$sig.cors.trait[1:10] # only first 10 displayed
```

    ## $ACIN1
    ## ME19 ME28  ME8  ME9 ME10 ME25 
    ##    5    8   10   11   25   31 
    ## 
    ## $AGGF1
    ## ME28  ME7  ME2 
    ##    8   22   23 
    ## 
    ## $AQR
    ## ME19 ME21 ME28 ME17 ME23  ME2 
    ##    5    7    8   20   21   23 
    ## 
    ## $ARGLU1
    ## ME16  ME3 
    ##    2   15 
    ## 
    ## $BAG2
    ## ME16 ME12  ME3  ME5  ME1 
    ##    2   13   15   26   29 
    ## 
    ## $BCAS1
    ## ME20 ME28 ME10  ME5 
    ##    1    8   25   26 
    ## 
    ## $BCAS2
    ## ME16  ME9 ME23  ME2  ME4  ME1 
    ##    2   11   21   23   24   29 
    ## 
    ## $BUB3
    ## ME8 ME9 ME7 ME2 ME4 ME5 
    ##  10  11  22  23  24  26 
    ## 
    ## $BUD13
    ## ME14 ME19 ME21 ME18 ME22 ME17 ME23  ME2 ME26  ME1 
    ##    4    5    7   17   18   20   21   23   28   29 
    ## 
    ## $BUD31
    ## ME20 ME14  ME0  ME8 ME12  ME7  ME2 ME10  ME5  ME1  ME6 
    ##    1    4    6   10   13   22   23   25   26   29   30

``` r
# Module-trait correlation coefficients and P-values:
J1pass$moduleTraitCor[,c(1:10)]
```

    ##             ACIN1        AGGF1          AQR       ARGLU1         BAG2
    ## ME20 -0.094187007 -0.079875341  0.103057540 -0.084712531  0.184627441
    ## ME16  0.126069973  0.125829436  0.338210955 -0.376360135  0.474339873
    ## ME15  0.243069048  0.226191344  0.225742639  0.149470677 -0.341869191
    ## ME14  0.223856908  0.117045011  0.345740870 -0.060061133  0.104749586
    ## ME19  0.684845147  0.031569628  0.506309064 -0.151149183  0.027515334
    ## ME0   0.247170248 -0.180359608 -0.140926106  0.069882844  0.142406195
    ## ME21  0.174891789  0.031564261  0.383094901 -0.124989088  0.160144243
    ## ME28  0.471644299  0.368091825  0.398830803 -0.077019991 -0.287583202
    ## ME24  0.103256854 -0.004519532  0.250607016 -0.194724351 -0.302337422
    ## ME8   0.449308783  0.026993582  0.305783000 -0.129593259 -0.359518739
    ## ME9   0.573457777 -0.213310151  0.098449971  0.061249011 -0.293339957
    ## ME29  0.317518097 -0.082450248  0.133917960  0.077739147 -0.287553535
    ## ME12  0.175492079  0.033631153  0.124934572 -0.039409627 -0.417383505
    ## ME13 -0.243129196  0.010254087 -0.031415695  0.187689029 -0.187542424
    ## ME3  -0.144160127  0.011575137 -0.120241864  0.394503834 -0.558334627
    ## ME30  0.048156448  0.024933861 -0.248329470  0.090333536 -0.288592452
    ## ME18 -0.054355839 -0.041905047 -0.351361550  0.194504785  0.012240492
    ## ME22  0.069960063 -0.012259707 -0.254955535  0.096372081  0.008950175
    ## ME31  0.093180213  0.047151134 -0.119823182 -0.208025937 -0.192216657
    ## ME17 -0.078960714 -0.346169303 -0.377136762 -0.061396350  0.168675149
    ## ME23  0.183746287 -0.242457340 -0.384934483 -0.001347678 -0.035949986
    ## ME7   0.204561585 -0.419069811 -0.171207053 -0.111378829  0.106850776
    ## ME2  -0.151710835 -0.480170237 -0.547714102  0.073629160 -0.050249492
    ## ME4   0.079198346 -0.323731363 -0.233958620 -0.074651633 -0.070472876
    ## ME10  0.421242304  0.084588457  0.273460762 -0.156744454  0.011293134
    ## ME5   0.352272741  0.146359616  0.003781174  0.114389694 -0.512887900
    ## ME11 -0.150156070 -0.094566937 -0.173637718  0.226521429 -0.254513104
    ## ME26 -0.153276274  0.017114139 -0.345141121  0.201171836 -0.023817382
    ## ME1  -0.001826104 -0.194361100 -0.284137744  0.210736823 -0.421367891
    ## ME6  -0.201732338  0.157294691  0.092513744 -0.033744113 -0.124718915
    ## ME25 -0.618527210  0.031132130 -0.262089771  0.162372896 -0.021662140
    ## ME27  0.110741164 -0.021406347  0.143499899 -0.033098912  0.300609740
    ##             BCAS1        BCAS2         BUB3        BUD13        BUD31
    ## ME20 -0.419622847  0.196918572  0.311627003  0.366074398  0.633608568
    ## ME16  0.294374914  0.407285851  0.277238396  0.325643598 -0.004556284
    ## ME15 -0.056297280 -0.018853025 -0.183092638  0.137846779  0.215915883
    ## ME14 -0.196672710  0.249566053  0.132162469  0.589420512  0.393537815
    ## ME19 -0.001844834 -0.050150402 -0.015933364  0.374126680  0.047405024
    ## ME0  -0.296625838 -0.065840262  0.104166821  0.008623718  0.510934892
    ## ME21 -0.014169618  0.194779564  0.173661546  0.438139067  0.201513264
    ## ME28  0.630038368 -0.021985600 -0.191515260  0.153129652 -0.281317205
    ## ME24  0.086872982  0.004656436 -0.330895444  0.321509843  0.091009461
    ## ME8   0.310105871 -0.350021573 -0.448492175  0.005127489 -0.526701842
    ## ME9  -0.068441911 -0.533008217 -0.442312760 -0.005521309 -0.037926863
    ## ME29  0.076625305 -0.321609412 -0.149093392 -0.114575132 -0.301348449
    ## ME12 -0.214239209  0.052423974 -0.018823537  0.352189726  0.418334202
    ## ME13 -0.202445807  0.061287979  0.220185646  0.081040463  0.077675299
    ## ME3  -0.124530472 -0.134770123 -0.357502805  0.085074630 -0.070560625
    ## ME30 -0.024384725 -0.091762975 -0.233417798 -0.048112615  0.034745785
    ## ME18 -0.166103309 -0.050645621  0.268136808 -0.401105917  0.187740583
    ## ME22 -0.276463273 -0.177649631  0.084863326 -0.415328354  0.161064141
    ## ME31 -0.286970779 -0.224568605 -0.194717581 -0.365205753 -0.077371549
    ## ME17 -0.053164012 -0.175700192 -0.161224907 -0.396723765 -0.328096962
    ## ME23 -0.068907490 -0.413144300 -0.246844917 -0.573898924 -0.223313467
    ## ME7   0.055707445 -0.253342992 -0.455579518 -0.242942736 -0.416963911
    ## ME2  -0.202510101 -0.477939050 -0.490436353 -0.563212301 -0.372929614
    ## ME4  -0.219591529 -0.428881826 -0.621465729 -0.228480460 -0.334942099
    ## ME10  0.767819952  0.107936312 -0.053112533  0.022359808 -0.492583261
    ## ME5   0.560563951 -0.215846983 -0.385250506 -0.359050018 -0.423819031
    ## ME11 -0.081650063 -0.039801443 -0.063740583 -0.257595906  0.022354225
    ## ME26  0.113043957 -0.185370942  0.007698684 -0.600202777 -0.261603896
    ## ME1   0.129299444 -0.426522523 -0.312920131 -0.525158845 -0.393912865
    ## ME6   0.295845740  0.029692149  0.076870658 -0.274300673 -0.523043190
    ## ME25  0.067959415  0.181555281  0.282259617 -0.212227145 -0.081945309
    ## ME27  0.095065855  0.299743346  0.317378696  0.148468043 -0.057875557

``` r
J1pass$moduleTraitPvalue[,c(1:10)]
```

    ##             ACIN1       AGGF1         AQR     ARGLU1        BAG2        BCAS1
    ## ME20 6.269782e-01 0.680427163 0.594733109 0.66217577 0.337663672 2.344606e-02
    ## ME16 5.146256e-01 0.515434632 0.072731651 0.04418979 0.009330264 1.211169e-01
    ## ME15 2.038900e-01 0.238053674 0.239010959 0.43900152 0.069496989 7.717669e-01
    ## ME14 2.430620e-01 0.545401325 0.066196095 0.75694777 0.588664537 3.065168e-01
    ## ME19 4.165264e-05 0.870857244 0.005071119 0.43381989 0.887330124 9.924220e-01
    ## ME0  1.961252e-01 0.349143409 0.465892128 0.71868430 0.461173335 1.181661e-01
    ## ME21 3.641877e-01 0.870879018 0.040238937 0.51826604 0.406626979 9.418437e-01
    ## ME28 9.797457e-03 0.049455751 0.032103125 0.69128507 0.130349755 2.495617e-04
    ## ME24 5.940169e-01 0.981436574 0.189777683 0.31142897 0.110917570 6.540841e-01
    ## ME8  1.447904e-02 0.889453522 0.106707996 0.50284672 0.055427688 1.015969e-01
    ## ME9  1.145748e-03 0.266555756 0.611393153 0.75228801 0.122491788 7.242589e-01
    ## ME29 9.326398e-02 0.670688902 0.488574688 0.68854462 0.130391183 6.927907e-01
    ## ME12 3.625177e-01 0.862500909 0.518449988 0.83915801 0.024274870 2.644290e-01
    ## ME13 2.037747e-01 0.957897453 0.871481763 0.32957122 0.329956007 2.922477e-01
    ## ME3  4.556143e-01 0.952479024 0.534402087 0.03419454 0.001646008 5.198145e-01
    ## ME30 8.040826e-01 0.897843239 0.193968018 0.64120287 0.128946135 9.000819e-01
    ## ME18 7.794419e-01 0.829117368 0.061622345 0.31198560 0.949750769 3.891541e-01
    ## ME22 7.183860e-01 0.949671989 0.181952658 0.61896958 0.963247389 1.465620e-01
    ## ME31 6.306823e-01 0.808096494 0.535836588 0.27886248 0.317823012 1.312069e-01
    ## ME17 6.838985e-01 0.065838443 0.043719043 0.75171062 0.381748454 7.841637e-01
    ## ME23 3.400149e-01 0.205066015 0.039210414 0.99446415 0.853118941 7.224562e-01
    ## ME7  2.871254e-01 0.023648563 0.374538055 0.56515293 0.581166358 7.740965e-01
    ## ME2  4.320935e-01 0.008383528 0.002101724 0.70425729 0.795741613 2.920912e-01
    ## ME4  6.829960e-01 0.086686376 0.221886700 0.70033689 0.716405667 2.523923e-01
    ## ME10 2.286123e-02 0.662641611 0.151178489 0.41679009 0.953635524 1.164466e-06
    ## ME5  6.090471e-02 0.448693575 0.984468892 0.55461710 0.004440737 1.562082e-03
    ## ME11 4.368816e-01 0.625582693 0.367691278 0.23735109 0.182738316 6.737097e-01
    ## ME26 4.273016e-01 0.929785490 0.066699295 0.29535968 0.902395479 5.593150e-01
    ## ME1  9.924990e-01 0.312350189 0.135225783 0.27250417 0.022816377 5.038238e-01
    ## ME6  2.939879e-01 0.415135735 0.633139198 0.86204344 0.519177954 1.191827e-01
    ## ME25 3.484392e-04 0.872632410 0.169608780 0.40004122 0.911191531 7.261287e-01
    ## ME27 5.673959e-01 0.912236185 0.457702661 0.86465698 0.113074294 6.237520e-01
    ##            BCAS2         BUB3        BUD13        BUD31
    ## ME20 0.305900350 0.0998428295 0.0508130483 0.0002244111
    ## ME16 0.028314805 0.1453868677 0.0847348517 0.9812856454
    ## ME15 0.922671465 0.3417653977 0.4757900829 0.2606189580
    ## ME14 0.191685021 0.4943430557 0.0007666807 0.0346761332
    ## ME19 0.796136011 0.9346192565 0.0455661021 0.8070823501
    ## ME0  0.734358885 0.5907515946 0.9645870908 0.0046204734
    ## ME21 0.311289102 0.3676245238 0.0174414920 0.2945236123
    ## ME28 0.909870734 0.3196258237 0.4277492097 0.1393152805
    ## ME24 0.980874356 0.0795473977 0.0889962956 0.6386985766
    ## ME8  0.062689657 0.0146804006 0.9789399791 0.0033316789
    ## ME9  0.002910364 0.0162808211 0.9773228473 0.8451359171
    ## ME29 0.088891774 0.4401708163 0.5539711784 0.1121483577
    ## ME12 0.787099313 0.9227920571 0.0609698160 0.0239201286
    ## ME13 0.752135284 0.2510787900 0.6760143446 0.6887877681
    ## ME3  0.485786960 0.0569110216 0.6608169422 0.7160669897
    ## ME30 0.635911399 0.2229877131 0.8042575500 0.8579887585
    ## ME18 0.794165438 0.1596195953 0.0310456526 0.3294359743
    ## ME22 0.356552571 0.6616097607 0.0250564260 0.4039012826
    ## ME31 0.241527791 0.3114461289 0.0514064441 0.6899449300
    ## ME17 0.361939719 0.4034259856 0.0331082514 0.0822803726
    ## ME23 0.025909389 0.1967335992 0.0011334066 0.2442378786
    ## ME7  0.184827593 0.0130079868 0.2041324910 0.0244327978
    ## ME2  0.008735785 0.0069130149 0.0014672273 0.0463176574
    ## ME4  0.020261167 0.0003203793 0.2332097527 0.0757192085
    ## ME10 0.577309172 0.7843677856 0.9083430096 0.0066348573
    ## ME5  0.260774818 0.0390358606 0.0557698348 0.0219553550
    ## ME11 0.837579809 0.7425419531 0.1773130609 0.9083657979
    ## ME26 0.335687457 0.9683837143 0.0005775188 0.1704301537
    ## ME1  0.021036779 0.0983696846 0.0034423895 0.0344885153
    ## ME6  0.878479426 0.6918546064 0.1498766885 0.0035993082
    ## ME25 0.345903947 0.1379390158 0.2690488387 0.6725945238
    ## ME27 0.114167502 0.0934157679 0.4421126249 0.7655429590

In order to explore module-trait significant associations found by
**JCNA1pass**, users can use the **JCNAModTrait()** function to plot
Module membership vs. Junction significance.

Here we explore the association between module No. 4 (turquoise) and
DDX1 splicing factor expression:

``` r
jMT <- JCNAModTrait(J1pass, trait ="DDX1",module = "4", cor.method = "bicor")

jMT$MMplot
```

<img src="
    
    /tutorial_MMvsJCplot.png" width="855" />

***jMT$MMplot*** shows an interactive scatterplot of Junction
Significance (JS) for DDX1 expression vs. Module Membership (MM) in the
brown module. There is a significant correlation between JS and MM,
indicating that junctions in this co-expression module are significantly
associated to the expression of this particular splicing factor. It also
illustrates how junctions highly significantly associated with a trait
are often also the most important (central) elements of modules
associated with the trait.

## 4.3 2-pass JCNA

Users have the option to continue into a second round of network
construction, that is specifically conceived to identify and remove
junction from the network, whose association to external trait data is
not splicing-specific but rather representing the association between
the trait and the total expression of the respective gene. This is
particularly relevant, since a considerable number of co-expressing
junctions are expected to cluster into single modules as a result of
intrinsic associations at the gene expression level without any specific
association to splicing.

For 2-pass JCNA, gene expression-based networks including correlations
with a user-selected sample trait are calculated using
**JCNAgenePrepare()** function. The required input data are:

1.  **JCNA1pass** output object
2.  Path to gtf file (used to asign junctions to genes)
3.  Normalized gene expression data

The absolute value of junction significance, which represents the
correlation coefficient between a given junction and the selected trait
is plotted as a function of the corresponding gene significance.
Junctions outside of the distribution by ≥2 standard deviations (showing
no correlation between junction and gene significance for trait) are
kept for network re-construction using **JCNA2pass()** function:

``` r
# Indicate path to gtf file:
gtf0 <- "~/Downloads/gencode.v32lift37.annotation.gtf" # gtf downloaded from GENCODE website

# load gene expression data:
gexp <- system.file("extdata", "genExpr.rds", package = "DJExpress")
filgenExpr <- readRDS(gexp)

# Change format of rownames for match:
rownames(J1pass$datTraits) <- paste0("TCGA_",
seq(1,nrow(J1pass$datTraits), 1))


# Run JCNAgenePrepare:
JgPrep <- JCNAgenePrepare(pass1.out = J1pass, genExpr = filgenExpr,
                          gtf=gtf0, networkType = "unsigned", cor.method = "bicor")
```

    ## [1] "Calculating co-expression networks from gene expression"
    ## Allowing parallel execution with up to 2 working processes.
    ##  Flagging junctions and samples with too many missing values...
    ##   ..step 1
    ##   ..Excluding 7 genes from the calculation due to too many missing samples or zero variance.
    ##   ..step 2
    ## Removing genes: CT45A3, CT45A4, CT45A5, CT45A6, GAGE13, GAGE2E, MAGEB1

    ## pickSoftThreshold: will use block size 2504.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 2504 of 17863
    ##    ..working on genes 2505 through 5008 of 17863
    ##    ..working on genes 5009 through 7512 of 17863
    ##    ..working on genes 7513 through 10016 of 17863
    ##    ..working on genes 10017 through 12520 of 17863
    ##    ..working on genes 12521 through 15024 of 17863
    ##    ..working on genes 15025 through 17528 of 17863
    ##    ..working on genes 17529 through 17863 of 17863
    ##    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
    ## 1      1 0.000835 -0.134          0.969 3710.000  3.68e+03 5630.0
    ## 2      2 0.272000 -1.480          0.944 1180.000  1.14e+03 2600.0
    ## 3      3 0.597000 -1.960          0.941  471.000  4.28e+02 1470.0
    ## 4      4 0.760000 -2.180          0.938  216.000  1.83e+02  939.0
    ## 5      5 0.828000 -2.270          0.928  111.000  8.61e+01  649.0
    ## 6      6 0.856000 -2.250          0.918   61.600  4.33e+01  473.0
    ## 7      7 0.885000 -2.170          0.924   36.700  2.30e+01  358.0
    ## 8      8 0.901000 -2.060          0.928   23.100  1.28e+01  279.0
    ## 9      9 0.919000 -1.950          0.937   15.200  7.45e+00  222.0
    ## 10    10 0.926000 -1.860          0.942   10.500  4.46e+00  179.0
    ## 11    12 0.955000 -1.680          0.963    5.430  1.73e+00  122.0
    ## 12    14 0.966000 -1.560          0.974    3.120  7.40e-01   87.7
    ## 13    16 0.947000 -1.530          0.964    1.930  3.39e-01   66.4
    ## 14    18 0.941000 -1.490          0.966    1.260  1.64e-01   51.9
    ## 15    20 0.912000 -1.480          0.951    0.866  8.21e-02   41.3
    ## 16    22 0.892000 -1.470          0.949    0.615  4.31e-02   33.4
    ## 17    24 0.891000 -1.450          0.956    0.450  2.32e-02   27.4
    ## 18    26 0.890000 -1.440          0.961    0.337  1.30e-02   22.8
    ## 19    28 0.895000 -1.420          0.968    0.258  7.36e-03   19.1
    ## 20    30 0.897000 -1.410          0.970    0.201  4.23e-03   16.2

    ## [1] "Constructing the gene network and modules"
    ## Allowing parallel execution with up to 2 working processes.

    ## [1] "Defining gene significance for each trait"
    ## [1] "Defining gene-Trait vs junction-Trait correlation"

``` r
# summary of JCNAgenePrepare output:
summary(JgPrep)
```

    ##                  Length Class        Mode
    ## Gene.sampletree    3    recordedplot list
    ## Gene.hc            3    recordedplot list
    ## Gene.NetTop        3    recordedplot list
    ## Gene.net          10    -none-       list
    ## Gene.ModuleTrait   3    recordedplot list
    ## Genetrait        720    data.frame   list
    ## GeneToJunct       11    data.frame   list
    ## JunctGeneTrait   360    -none-       list

``` r
J2pass <- JCNA2pass(pass1.out = J1pass, GenePrepare.out = JgPrep,
                     workDir= getwd(),
                     trait="DDX1", cor.method = "bicor")
```

    ## [1] "Calculating co-expression networks from gene expression"
    ## Allowing parallel execution with up to 2 working processes.
    ##  Flagging junctions and samples with too many missing values...
    ##   ..step 1
    ##  Flagging junctions and samples with too many missing values...
    ##   ..step 1

    ## [1] "Choosing the soft-thresholding power"
    ## pickSoftThreshold: will use block size 5257.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5257 of 5257
    ##    Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
    ## 1      1   0.0859  0.862          0.902 1230.000  1.23e+03 1810.00
    ## 2      2   0.1970 -0.751          0.956  433.000  4.18e+02  890.00
    ## 3      3   0.5500 -1.250          0.988  188.000  1.71e+02  517.00
    ## 4      4   0.7120 -1.530          0.995   93.600  7.95e+01  330.00
    ## 5      5   0.7980 -1.660          0.999   51.200  4.06e+01  225.00
    ## 6      6   0.8420 -1.770          0.998   30.300  2.25e+01  160.00
    ## 7      7   0.8710 -1.860          0.998   19.000  1.32e+01  118.00
    ## 8      8   0.8940 -1.910          0.998   12.600  8.20e+00   90.20
    ## 9      9   0.9050 -1.940          0.995    8.740  5.37e+00   70.40
    ## 10    10   0.9100 -1.960          0.990    6.300  3.73e+00   56.30
    ## 11    12   0.9140 -1.990          0.980    3.610  1.96e+00   38.00
    ## 12    14   0.8980 -2.000          0.973    2.300  1.09e+00   27.00
    ## 13    16   0.8980 -1.920          0.974    1.580  6.46e-01   20.00
    ## 14    18   0.9080 -1.820          0.987    1.160  4.03e-01   15.40
    ## 15    20   0.9280 -1.700          0.982    0.886  2.62e-01   12.10
    ## 16    22   0.9430 -1.610          0.980    0.701  1.73e-01   10.10
    ## 17    24   0.9330 -1.640          0.979    0.570  1.16e-01    9.31
    ## 18    26   0.9400 -1.630          0.989    0.473  8.03e-02    8.59
    ## 19    28   0.9530 -1.600          0.991    0.399  5.61e-02    7.94
    ## 20    30   0.9480 -1.600          0.984    0.341  3.97e-02    7.36

    ## [1] "Constructing gene effect-independent junction network and modules"
    ## Allowing parallel execution with up to 2 working processes.
    ## [1] " modules identified:31"
    ## 
    ##         black          blue         brown          cyan     darkgreen 
    ##            47           143           139            27            23 
    ##      darkgrey    darkorange       darkred darkturquoise         green 
    ##            23            21            23            23            50 
    ##   greenyellow          grey        grey60     lightcyan    lightgreen 
    ##            33          3989            27            27            25 
    ##   lightyellow       magenta  midnightblue        orange          pink 
    ##            25            41            27            23            43 
    ##        purple           red     royalblue   saddlebrown        salmon 
    ##            36            48            24            20            32 
    ##       skyblue     steelblue           tan     turquoise         white 
    ##            20            20            32           175            20 
    ##        yellow 
    ##            51 
    ## [1] "Label 0 is reserved for genes outside of all modules"

    ## [1] "Calculating module membership values and trait correlations"

    ## Registered S3 method overwritten by 'quantmod':
    ##   method            from
    ##   as.zoo.data.frame zoo

    ## [1] "Calculating adjacency and saving Cytoscape tables"
    ## Allowing parallel execution with up to 2 working processes.
    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.
    ## [1] "done"

``` r
# summary of JCNA2pass output:
summary(J2pass)
```

    ##                  Length   Class        Mode     
    ## Junct.sampletree        3 recordedplot list     
    ## net                    10 -none-       list     
    ## MEs                    31 data.frame   list     
    ## module.den              2 -none-       list     
    ## JunctMM                31 data.frame   list     
    ## MMPvalue               31 data.frame   list     
    ## JunctTS                 5 data.frame   list     
    ## JSPvalue                5 data.frame   list     
    ## MMvsJunctSig           31 -none-       list     
    ## modTOM           15912121 -none-       numeric  
    ## modGenes.2           3989 -none-       character
    ## modGenes             3989 -none-       character
    ## modNames               31 -none-       character
    ## Cytoscape.input        31 -none-       list     
    ## VisANT.input           31 -none-       list

``` r
# Re-constructed network dendrogram:
J2pass$module.den
```

<img src="ReadFig/2passdendro.png" width="571" />

``` r
# Examples of junction modules significantly associated to DDX1 splicing factor expression:
J2pass$MMvsJunctSig$purple
```

<img src="ReadFig/purplemod.png" width="855" />

``` r
# Junctions in purple module (first 10 shown):
J2pass$net$colors[which(J2pass$net$colors=="purple")[1:10]]
```

    ##   chr2:74708214:74708349:2   chr2:99193607:99203938:1 
    ##                   "purple"                   "purple" 
    ## chr2:122120877:122122669:2 chr3:126741185:126741594:1 
    ##                   "purple"                   "purple" 
    ##   chr6:31927897:31927996:1   chr6:31930869:31931189:1 
    ##                   "purple"                   "purple" 
    ##   chr6:31933791:31934485:1 chr7:100279859:100279944:2 
    ##                   "purple"                   "purple" 
    ## chr7:100280106:100280212:2 chr7:100280854:100280926:2 
    ##                   "purple"                   "purple"

**JCNA2pass** output includes tables in the format required to be loaded
in network visualization tools such as Cytoscape and VisANT which are
saved as text files in the working directoy:

``` r
# Network edges data of purple module formatted for Cytoscape:
J2pass$Cytoscape.input$purple$edgeData[c(1:10),]
```

    ##    fromNode   toNode     weight  direction fromAltName toAltName
    ## 1   CCDC142   PLXNA1 0.02750298 undirected     CCDC142    PLXNA1
    ## 2   CCDC142   SKIV2L 0.03256128 undirected     CCDC142    SKIV2L
    ## 3   CCDC142 SKIV2L.1 0.03324313 undirected     CCDC142    SKIV2L
    ## 4   CCDC142 SKIV2L.2 0.02877763 undirected     CCDC142    SKIV2L
    ## 5   CCDC142   GIGYF1 0.03382747 undirected     CCDC142    GIGYF1
    ## 6   CCDC142 GIGYF1.1 0.02298797 undirected     CCDC142    GIGYF1
    ## 7   CCDC142 GIGYF1.2 0.02046708 undirected     CCDC142    GIGYF1
    ## 8   CCDC142 GIGYF1.3 0.03512255 undirected     CCDC142    GIGYF1
    ## 9   CCDC142 GIGYF1.4 0.03946027 undirected     CCDC142    GIGYF1
    ## 10  CCDC142 GIGYF1.6 0.02237523 undirected     CCDC142    GIGYF1

``` r
# Network nodes data of purple module formatted for Cytoscape:
J2pass$Cytoscape.input$purple$nodeData[c(1:10),]
```

    ##    nodeName altName nodeAttr[nodesPresent, ]
    ## 1   CCDC142 CCDC142                   purple
    ## 2    INPP4A  INPP4A                   purple
    ## 3    CLASP1  CLASP1                   purple
    ## 4    PLXNA1  PLXNA1                   purple
    ## 5    SKIV2L  SKIV2L                   purple
    ## 6  SKIV2L.1  SKIV2L                   purple
    ## 7  SKIV2L.2  SKIV2L                   purple
    ## 8    GIGYF1  GIGYF1                   purple
    ## 9  GIGYF1.1  GIGYF1                   purple
    ## 10 GIGYF1.2  GIGYF1                   purple

``` r
# Network visualization using Cytoscape:
```

<img src="ReadFig/cytoscapeex.png" width="1897" />

***Note***: On Cytoscape 3.8.0, you can import table network using
**File** -&gt; **Import** -&gt; **Network from File** and selecting the
module edge file. On the import dialogue box, select the following
settings:

**fromNode** column as the Source Node  
**toNode** column as the Target Node  
**weight** column should be left as an Edge Attribute  
**direction** column as interaction type  
**fromAltName** as Source Node Attribute  
**toAltName** as Target Node Attribute.

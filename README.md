# LINEAGE

## 1 Introduction

Lineage analysis based on scRNA-seq could provide key insights into the fate of individual cells and inform therapeutic intervention for diseases. However, such analysis is limited by several technical challenges. On top of considerable computational resources, these analyses also require specific types of matching data. Therefore, a computational algorithm for label-free lineage analysis based on endogenous markers is desired.

We chose mitochondrial RNA variant as endogenous makers for lineage analysis in consideration of relatively small mitochondrial genome, mitochondrial genome polymorphism and heritable mitochondrial variation. In order to select the clonal feature without pre-existing knowledge, we defined “cross entropy” as the cell distribution similarity among subspaces. These subspaces are obtained by hierarchical clustering according to the variant frequency dynamic patterns. Then we designed a lineage analysis algorithm called Label-free IdeNtification of Endogenous informAtive sinGle cell mitochondrial RNA mutation for lineage analysis and clonal Evolution (LINEAGE) by integrating the low cross-entropy subspaces identification with a consensus clustering method.

We have provided an example data and its validated cell lineage information in “TF1_clones.rda” in the *data* directory, which might be helpful for users to go through this manual.

## 2 Install

LINEAGE can be successfully installed on MacOS*, Windows, CentOS, but not on Ubuntu Linux*.

*\*: checked by R-hub builder (https://builder.r-hub.io/).*

There are 2 options to install LINEAGE. 

**• Option 1: using the devtools package to install directly from GitHub**

```{r option1, eval=FALSE}
# install.packages("devtools")
devtools::install_github("songjiajia2018/LINEAGE")
```

**• Option 2: Download the source code from github, and install as follow**

The source code package is provided in the *source_pkg* directory.
```{r option2, eval=FALSE}
install.packages("the_directory_that_contain_the_package/LINEAGE_0.99.1.tar.gz",repos=NULL,type="source")
```

## 3 Example data

An example data is available in the source codes.

**TF1_clones.rda** is a list containing the input mitochondrial genotype matrix (TF1_clones\$data) and the experimental validated cell lineage information (TF1_clones\$rlabel).

## 4 Analysis

We have provided an example data **“TF1_clones.rda”** containing a mitochondrial genotype matrix of TF1_clones and its validated cell lineage information in the *data* dictionary of LINEAGE package, which can be used for test. It should be noted that **the result has a certain randomness** because of the randomness from clustering and dimension reduction processes.

#### 4.0 preprocessing: generate mitochondrial genotype matrixes
4.0.1 After alignment, new bam files consisting of MtDNA records, which were extracted from the alignment result with SAMtools, were obtained. 
```{r mtbam, eval=FALSE}
STAR --runThreadN {RUNNING_THREADS_NUMBER} --genomeDir {GENOME_PATH} --outFileNamePrefix {PREFIX} --sjdbGTFfile {GTF} --outSAMunmapped Within --readFilesIn {FASTQ1} {FASTQ2}
samtools view -bS {STAR_OUT_SAM} > {OUT_BAM}
samtools sort {INPUT_BAM} -o {OUT_SORTED_BAM}
samtools view -h {INPUT_SORTED_BAM} {REGION: MT, chrM, et al.} > {FINAL_SAM}
samtools view -bS {FINAL_SAM} >{FINAL_BAM}
```
4.0.2 The total number of reads aligned to per allele on each site of mitochondrial genome were counted using *https://github.com/songjiajia2018/ppl*. The final output is a rds file containing mitochondrial varients frequency.
```{r mtmatrix, eval=FALSE}
python ppl/ppl2_run.py -p -m -r --input {FILELIST} --input-filelist
```
(i) Input is a file list: The input should be a csv table with the sorted bam files and their output prefixes. Here is an example:
```{r input, eval=FALSE}
SRR3562459_2.bam,SRR3562459
SRR3562814_2.bam,SRR3562814
SRR3563095_2.bam,SRR3563095
SRR3563458_2.bam,SRR3563458
```
The output contains five mutation files for each bam file:
```{r output, eval=FALSE}
SRR3562459.A.txt SRR3562459.coverage.txt SRR3562459.C.txt SRR3562459.G.txt SRR3562459.T.txt
SRR3562814.A.txt SRR3562814.coverage.txt SRR3562814.C.txt SRR3562814.G.txt SRR3562814.T.txt
SRR3563095.A.txt SRR3563095.coverage.txt SRR3563095.C.txt SRR3563095.G.txt SRR3563095.T.txt
SRR3563458.A.txt SRR3563458.coverage.txt SRR3563458.C.txt SRR3563458.G.txt SRR3563458.T.txt
```
(ii) Parameters:
```{r options, eval=FALSE}
-p: calling variations from bam files and generating five txt files with 'A', 'T', 'C', 'G', 'coverage' suffixes, respectively
-m: merging all the 'A', 'T', 'C', 'G', 'coverage' txt files in a directory to five files with ".gz" format for the following rds generation
-r: generating rds file for downstream analysis from a directory containing merged 'A', 'T', 'C', 'G', 'coverage' txt files
--name: The program will create a directory named by the argument. Sample processing and downstream analysis will be performed under the directory. (default: PROJECT_MITO)
--qalign: specify minimum alignment quality required to be considered (default: 30)
--maxBP: specify maximum length of mtDNA genome (default: 16569, for mt.fa)
--reference: specify the mtDNA reference (default:./ppl/mito_reference/mt.fa)
```
4.0.3 Extract mitochondrial genotype matrix (mtMatrix.txt), where a column represented a single cell and a row represented variants frequency of a specific mitochondrial genotype, from a rds file using *https://github.com/songjiajia2018/ppl/mtMatrix.R*
```{r mtmatrix, eval=FALSE}
Rscript mtMatrix.R {RDS_FILE}
```

#### 4.1 mode1: run parallel iterative optimization
```{r mode1, eval=FALSE}
data("TF1_clones")
data=TF1_clones$data
result=lineage(data = data, repeats = 30, thread = 20)
```

#### 4.2 mode2: run non-parallel iterative optimization
```{r mode2}
data("TF1_clones")
data=TF1_clones$data
result=lineage(data = data, repeats = 30, thread = NULL)
```

#### 4.3 trace the lineage tree
```{r lineage_tree}
data("TF1_clones")
data=TF1_clones$data
result=lineage(data = data, repeats = 30, thread = 20)

# The inferred clone labels are embedded in the returned result
# as result$label. We also provide lineage_tree function to trace # the lineage tree of the returned result.
hc=lineage_tree(result)
str(hc, max.level = 1)
```

#### 4.4 visualization of the result
```{r plots}
data("TF1_clones")
data=TF1_clones$data
label=TF1_clones$rlabel    #reference clone labels
result=lineage(data = data, repeats = 30, thread = 20)

plots0=traceplot(result, label)  #plots with reference clone labels
plots=traceplot(result, result$label)  #plots with inferred clone labels
```
```{r plots_output, eval=FALSE}
# 2d visualization cluster results with clonal labels
# colored in reference clone labels 
print(plots0$d2d)
# or inferred labels
print(plots$d2d)
# Heatmap results with markers information across cells and color bar is
# colored in reference clone labels
print(plots$heatmap)
```

#### 4.5 results
```{r result, eval=FALSE}
# The recommended result is embedded in the returned result as result$best.
best=list(result=result$best, plots=plots)
```

## 5 Suggestions

*lineage*: Considering the randomness from clustering and dimension reduction processes, an iteration process with at least **30 repeats(default)** is recommended to guarantee a more stable and reliable cell clustering/clone tracing result.

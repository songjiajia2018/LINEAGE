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

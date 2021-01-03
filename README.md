## Independent Component Analysis based gene co-expression Network inference (ICAnet) to decipher functional modules for better single cell clustering and batch integration
![ICAnet](https://github.com/WWXkenmo/ICAnet/blob/master/figure-m1.png)

## ICAnet
We introduced independent component analysis (ICA) into single cell clustering to decompose the gene expression matrix into a number of independent components. Each independent component was characterized by a co-expression pattern and was associated with certain meaningful biological pathway. Such concept enables ICAnet to identify shared gene co-expression module across different batches of datasets. Based on the idea that different batches of scRNA-seq datasets derived from the same cell type don’t exhibit exactly the same gene expression patterns but the key co-expression module usually tends to keep similar, ICAnet pairs the same sub-population of cells among different batches, regardless of their library type, sequencing platform or other influences. These features of ICAnet make it performs better in cell clustering and integrative analysis with different batches of scRNA-seq datasets.

## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Integrating Cell Line Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/ICAnet_tutorial2.html)

* [Integrating Multiple Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/Pancreas_Tutorial2.html)

* [Using ICAnetTF to identify TF-regulons in Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/MouseBrain_TF2.html)

## Installation

1. Install [R](https://www.r-project.org/)  (>= 3.5)
2. Install [Rstudio](https://www.rstudio.com/products/rstudio/download/) (recommended)
3. Install all the packages required (see Dependency session)
4. Download ICAnet.tar.gz
5. Use the following R commands.
```
install.packages("~/ICAnet.tar.gz",repos=NULL, type="source",INSTALL_opts=c("--no-multiarch"))
```
## Dependency
Packages from Bioconductor: AUCell, RcisTarget, MineICA, STRINGdb, isva, clusterProfiler(required by RSCORE), genesorteR(required by RSCORE)
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## Required
BiocManager::install(c("AUCell", "STRINGdb", "MineICA","RcisTarget","isva","clusterProfiler"))
BiocManager::install("mahmoudibrahim/genesorteR") 
```
Packages from CRAN: Seurat, cluster, coop, fastICA, ica, igraph, isva, pheatmap, rARPACK, networkD3, doMC,propr(required by RSCORE), network(required by RSCORE), intergraph(required by RSCORE)
```
install.packages(c("Seurat", "cluster", "coop", "fastICA", "ica", "igraph", "isva", "pheatmap", "rARPACK", "networkD3"))
install.packages("doMC", repos="http://R-Forge.R-project.org")
install.packages(c("propr", "network","intergraph"))
```
Packages from Github: RSCORE
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("wycwycpku/RSCORE")
```
## Network Dependency
ICAnet required PPI network or cisTarget feather file as input.
### 1.PPI network (Required by ICAnet)
The PPI network could be downloaded through getPPI_String/getPPI_Biogrid in RSCORE
```
library(RSCORE)
PPI.network_STRING <- getPPI_String(data, species=9606) #(PPI network of STRING Database)
PPI.network_BioGrid <- getPPI_Biogrid(data, species=9606) #(PPI network of BioGRID Database)
```
In which the 9606 is the NCBI taxon-Id for Homo sapiens, for the taxon-id of other species, see https://string-db.org/cgi/input.pl?input_page_active_form=organisms
Meanwhile, the data could be Seurat object or matrix with gene symbol as row names.
### 2.TF-motif & motif ranking (Required by ICAnetTF)
ICAnetTF required TF-motif binding information and motif gene annotation matrix (.feather), both provided by RcisTarget, the user could download feather file from  https://resources.aertslab.org/cistarget/ with following commond
```
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc8nr/gene_based/mm9-500bp-upstream-10species.mc8nr.feather
```
## Tutorial Data Instruction
### 1.Integrating Cell Line Single-Cell RNA-seq Dataset
The normalized single cell cell line dataset could be downloaded from 
https://github.com/WWXkenmo/MouseGerm/blob/master/cell_line_exp.RDS, also, the annotation of this dataset could be downloaded from 
https://github.com/WWXkenmo/MouseGerm/blob/master/ID_for_cell_line.txt
### 2.Integrating Multiple Single-Cell RNA-seq Dataset
Three pancreas dataset could be downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/
Including baron-human.rds, muraro.rds and segerstolpe.rds
### 3.Using ICAnetTF to identify TF-regulons in Single-Cell RNA-seq Dataset
The mouse brain gene expression dataset could be downloaded from GEO:GSE60361, and the motif annotation dataset (feather file) could be downloaded from https://resources.aertslab.org/cistarget/, meanwhile, the cell annotation dataset could be downloaded from https://github.com/WWXkenmo/MouseGerm/blob/master/mouse_brain_cell_annotation.csv

## Independent Component Analysis based gene co-expression Network inference (ICAnet) to decipher functional modules for better single cell clustering and batch integration
![ICAnet](https://github.com/WWXkenmo/ICAnet/blob/master/figure-m1.png)

# ICAnet
We introduced independent component analysis (ICA) into single cell clustering to decompose the gene expression matrix into a number of independent components. Each independent component was characterized by a co-expression pattern and was associated with certain meaningful biological pathway. Such concept enables ICAnet to identify shared gene co-expression module across different batches of datasets. Based on the idea that different batches of scRNA-seq datasets derived from the same cell type don’t exhibit exactly the same gene expression patterns but the key co-expression module usually tends to keep similar, ICAnet pairs the same sub-population of cells among different batches, regardless of their library type, sequencing platform or other influences. These features of ICAnet make it performs better in cell clustering and integrative analysis with different batches of scRNA-seq datasets.

# Use
ICAnet containing three main functions
> ica <- ICAcomputing(obj, seurat.obj, batch, center=TRUE,scale=FALSE,RMT=TRUE, nc.vec=NULL, ICA.type="FastICA", nbIt=10, alg.type="deflation",
                          fun="logcosh",maxit =500, tol=10^-6, funClus="hclust",row.norm=FALSE,boostrap=FALSE,method="average",two.stage=TRUE)
* **obj**
a Seurat or matrix data type object, containing scRNA-seq dataset
* **seurat.obj**
binary variable, TRUE: the obj is Seurat object
* **batch**
if object is not Seurat object, please give a vector which annotate batch information of each cell in obj. We also annotate that for Seurat.obj, the user is required to add a variable named "batch" in their metadata, which containing batch information.
* **center**
whether to perform centerelization on gene expression space.
* **scale**
whether to perform scaling( zero mean and unit variance) on gene expression space
* **RMT**
whether to use RMT to estimate the number of independent component
* **nc.vec**
if user doesn't want to use RMT to estimate the number of component, please input the expect number of independent component for each batch
* **ICA.type**
ICAnet use two different implement to perform ICA, one is JADE, the other is fastICA. The fastICA implement is based on R package "MineICA"
* **alg.type,fun,funClus,method**
The parameters required in MineICA
* **two.stage**
whether to perform ICA+RMT to estimate the number of independent component, we note that this stategy will return the estimate which more close to the number of intrinsic cell types.

> basal <- CrossBatchGrouping(ica.pooling,k.max=(ncol(ica.pooling)-1),plot=TRUE,cor="pearson",W.top=2.5,filtering=TRUE,threshold=30)
* **ica.pooling**
The independent components of all batches
* **k.max**
The maximum number of clusters
* **plot**
whether to plot the cross batch expression programs grouping matrix
* **cor**
the correlation meric used when measure the connection among the different expression programs
* **W.top**
The threshold used to define activated genes
* **filtering**
whether to filter the independent components
* **threshold**
The minimum number of activated genes required for each independent component

> obj  <- RunICAnet(obj,ica.filter, W.top=2.5,PPI.net=NULL, PPI="String",species=9606,score=600,max.step=10,small.size=3,
                      aucMaxRank=3000,cores=6)
* **obj**
The single cell gene expression dataset, required Seurat object
* **ica.filter**
The basal programs
* **W.top**
The threshold used to define activated genes
* **PPI.net**
PPI network used for detect module
* **max.step**
The maximum random walk steps, if this parameters getting larger, the size of module getting bigger.
* **small.size**
The minimum genes required for each module
* **aucMaxRank**
The maximum number of genes used for calculating auc value.

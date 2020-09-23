## Independent Component Analysis based gene co-expression Network inference (ICAnet) to decipher functional modules for better single cell clustering and batch integration
![ICAnet](https://github.com/WWXkenmo/ICAnet/blob/master/figure-m1.png)
* [ICAnet tutorial](http://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/ICA_tutoral.nb.html)
## ICAnet
We introduced independent component analysis (ICA) into single cell clustering to decompose the gene expression matrix into a number of independent components. Each independent component was characterized by a co-expression pattern and was associated with certain meaningful biological pathway. Such concept enables ICAnet to identify shared gene co-expression module across different batches of datasets. Based on the idea that different batches of scRNA-seq datasets derived from the same cell type don’t exhibit exactly the same gene expression patterns but the key co-expression module usually tends to keep similar, ICAnet pairs the same sub-population of cells among different batches, regardless of their library type, sequencing platform or other influences. These features of ICAnet make it performs better in cell clustering and integrative analysis with different batches of scRNA-seq datasets.

## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Integrating Cell Line Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/ICAnet_tutorial.html)

* [Integrating Multiple Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/Pancreas_Tutorial.html)

* [Using ICAnetTF to identify TF-regulons in Single-Cell RNA-seq Dataset](https://htmlpreview.github.io/?https://github.com/WWXkenmo/ICAnet/blob/master/vignettes/MouseBrain_TF.html)

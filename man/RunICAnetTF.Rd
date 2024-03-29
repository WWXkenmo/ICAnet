% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICAnetTF.R
\name{RunICAnetTF}
\alias{RunICAnetTF}
\title{Using independent components to perform TF-target enrichment to detect TF-regulon for single cell clustering}
\usage{
RunICAnetTF(
  obj,
  ica.filter,
  W.top.TFs = 2.5,
  W.top.genes = 2.5,
  Motif_Net = NULL,
  TF_motif_annot = NULL,
  small.size = 3,
  aucMaxRank = 3000,
  cores = 6,
  verbose = FALSE,
  nMC = 100,
  cutoff = 0.01
)
}
\arguments{
\item{obj}{a Seurat object}

\item{ica.filter}{the filtered/unfiltered ica-components set}

\item{W.top.TFs}{the threshold to determine the activated TFs, the TFs which has absolute attributes value large than threshold*standard derivation from mean are the activated TFs (default: 2.5)}

\item{W.top.genes}{the threshold to determine the activated genes, the genes which has absolute attributes value large than threshold*standard derivation from mean are the activated genes (default: 2.5)}

\item{Motif_Net}{the matrix which indicating the boolean network of motif-target}

\item{TF_motif_annot}{Annotation to human/mouse/fly transcripton factors for the motifs in each motif collection (e.g. mc8nr or mc9nr)}

\item{small.size}{integer number to determine the minimum size of module. The module which has the number of gene member less than this value will be filtered}

\item{aucMaxRank}{Integer number. The number of highly-expressed genes to include when computing AUCell}

\item{cores}{the number of cores used for computation}

\item{verbose}{a boolean variable, whether show the running process (default: FALSE)}

\item{nMC}{the number of permutations, which is used for calculate the pvalue of each module (default: 100)}

\item{cutoff}{the significant level to filter the TF-regulons (default: 0.01)}
}
\value{
a Seurat object which contain the "IcaNet_TF" assay
}
\description{
ICAnet used independent components to perform TF-target enrichment test, and select significant TF-regulons for the following analysis
}

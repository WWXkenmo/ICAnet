###Run Cell Type Specific Modules
#' @title Module markers of identity classes
#' @description Finds markers (differentially activated modules) for identity classes
#'
#' @param obj a Seurat object
#' @param identity the identity id
#' @param assay which assay shall be tested (Default: "IcaNet")
#' @param MQC whether to control the statistical significance level of module (Default: TRUE)
#' @param pvalue integer value to indicate the threshold to filter the module (Default: 0.05)
#' @param rank.order rank the results based on which feature (Default: "myAUC")
#'
#' @return a matrix object which containing the significant module information
#' @export
#' @importFrom Seurat FindMarkers
#'
FindMarkerModule <- function(obj, identity,assay="IcaNet",MQC=TRUE,pvalue=0.05, rank.order="myAUC"){
diffmodule <- FindMarkers(obj, ident.1 =identity, test.use="roc",only.pos=TRUE,logfc.threshold = 0.01,assay=assay)
diffmodule$Module_significance <- obj@misc$modN.score[rownames(diffmodule)]
diffmodule <- diffmodule[order(diffmodule$myAUC,decreasing=TRUE),]
if(MQC) diffmodule <- diffmodule[diffmodule$Module_significance <= pvalue,]
diffmodule
}

###plot molecular network
#' @title Plot interactive molecular network
#'
#' @param gene_sets the module needed to plot the interaction network
#' @param network the molecular network
#' @param ica.score the ica-components
#' @param nodeSize the character value indicating node size is showed according to which information ("ica.score" or "degree")
#' @param fontSize the integer value indicating the size of font (Default: 15)
#' @param charge the integer value indicating the degree of expanding of network (Default: -100)
#' @param width width of the graph
#' @param height height of the graph
#' @param size_scale the size of the node (Default: 1)
#'
#' @return a interactive plot of molecular network
#' @export
#' @importFrom networkD3 forceNetwork
#' @importFrom networkD3 JS
#' @importFrom Matrix summary
#' @examples
plot_module <- function(gene_sets, network,ica.score, nodeSize="ica.score",fontSize = 15, charge=-100,width=500, height=500,size_scale=1){

subnet <- network[gene_sets,gene_sets]
subnet <- as.matrix(summary(subnet))
if(nodeSize=="ica.score"){importance <- ica.score[gene_sets]}else{if(nodeSize=="degree"){importance <- rowSums(network[gene_sets,gene_sets])}else{
cat("Error! Node size is not well-defined")
}}
NodeInfor <- data.frame(name = c(gene_sets),
                        group = c(rep("module",length(gene_sets))),
                        size = importance*size_scale)
NodeInfor$group <- as.character(NodeInfor$group)
NodeInfor$group <- as.factor(as.character(NodeInfor$group))
colnames(subnet) <- c("Source","Target","Value")
subnet[,1] <- subnet[,1]-1
subnet[,2] <- subnet[,2]-1
subnet <- as.data.frame(subnet)


forceNetwork(Links = subnet,Nodes = NodeInfor,Source="Source",Target="Target",
                       Value="Value",NodeID="name",Group = "group",Nodesize="size",
                       opacityNoHover = 2,opacity=2,fontFamily="Arial",fontSize = fontSize,charge = charge,
                       radiusCalculation = JS("d.nodesize"),
                       colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),width = width,
                       height = height)
}


#' @title Using SVD to decompose the module activity matrix
#'
#' @param obj a Seurat object
#' @param nu the number of singular vector
#' @param power integer value indicating the smoothness of signal on singular vector
#'
#' @return a Seurat object which containing reduction object "Module_SVD"
#' @export
#' @importFrom  coop tpcor
#' @importFrom rARPACK svds
#' @importFrom Seurat CreateDimReducObject
#' @importFrom Seurat GetAssayData
#' @examples
RunModuleSVD <- function(obj, nu=30,power=0.5){
moduleExp_scale <- coop::tpcor(t(as.matrix(GetAssayData(obj))))
svd.m <- rARPACK::svds(t(moduleExp_scale),k=nu,nu=nu)
svd.reduction <- (svd.m$u) %*% diag((svd.m$d)^power)
rownames(svd.reduction) <- colnames(moduleExp_scale)
obj[['Module_SVD']] <- CreateDimReducObject(embeddings = svd.reduction,key="svd_")
obj
}

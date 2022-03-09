#' @title Identify shared basal components
#' @description Identify shared basal components through using whiten matrices derived from latent semantic indexing
#'
#' @param data_list a list Seurat objects
#' @param var.names the features genes, e.g highly variable features learned from *FindVariableFeatures*
#' @param reduct.dim the reducted dimensions, default: 30
#' @param method truncted SVD solver. "irlba" or "rARPACK"
#' @param fastICA boolean variable,if TRUE, using FastICA, otherwise use JADE to perform ICA.
#' @param corr the correlation method which determine the correlation between different pair of source components (default: spearman)
#' @param W.top the threshold to determine the activated genes, the genes which has absolute attributes value large than threshold*standard derivation from mean are the activated genes (default: 2.5)
#'
#'
#' @return a list object which containing merged seurat object and basal components
#' @examples
#' \dontrun{
#' list <- SharedBasalComponents(PBMC_list,features,30)
#' objs <- list$mergedObject
#' basalcomp <- list$basalcomp
#' }
#' @importFrom cluster pam
#' @importFrom pheatmap pheatmap
#'
#' @export
#'
SharedBasalComponents <- function(data_list,var.names,reduct.dim,center=TRUE,method = "irlba",fastICA=FALSE,corr="spearman",W.top=2.5){
  ##scaling
  writeLines("Determine the reference atlas...")
  data_list <- ReferenceIntegrate(data_list,reduct.dim,var.names)
  
  for(i in 1:length(data_list)){
  if(dim(data_list[[i]]@assays$RNA@scale.data)[1]==0){
    data_list[[i]] <- ScaleData(data_list[[i]],scale.max = Inf,verbose=FALSE,features = var.names)
  }
  }
  
  writeLines("Perform scaling...")
  scaleDat <- list()
  if(!center) dat <- apply(data_list[[1]]@assays$RNA@scale.data[var.names,],2,scale)
  if(center) dat <- apply(data_list[[1]]@assays$RNA@scale.data[var.names,],2,scale,center=TRUE,scale=FALSE)
  if(method == "rARPACK") ref.svd <- rARPACK::svds(dat,k=reduct.dim)
  if(method == "irlba") ref.svd <- irlba::irlba(dat, nv = reduct.dim,center = TRUE)
  
  ##define encoder function
  writeLines("Define encoder...")
  encode.m <- diag(1/ref.svd$d) %*% t(ref.svd$u)
  decode.m <- ref.svd$u %*% diag(ref.svd$d)
  ##perform scaling to all matrix, calculate cor after projection and concatenate together
  writeLines("Identify shared co-variation...")
  cor <- NULL
  sampleName <- NULL
  for(i in 1:length(data_list)){
    if(!center) scale.data <- apply(data_list[[i]]@assays$RNA@scale.data[var.names,],2,scale)
    if(center){scale.data <- apply(data_list[[i]]@assays$RNA@scale.data[var.names,],2,scale,center=TRUE,scale=FALSE)}
    proj <- encode.m %*% scale.data
    scaleDat[[i]] <- scale.data
    cor.m <- ref.svd$v%*%proj
    sampleName <- c(sampleName,colnames(data_list[[i]]))
    cor <- cbind(cor, cor.m)
  }
  
  ##perform SVD and concatenate vector
  writeLines("Calculating loading matrix...")
  if(method == "rARPACK") loading <- rARPACK::svds(cor,k=reduct.dim)$v
  if(method == "irlba") loading <- irlba::irlba(cor, nv =reduct.dim)$v
  rownames(loading) <- sampleName
  ###
  writeLines("Project gene expression matrix on shared low dimensional space...")
  projData <- list()
  whitenData <- list()
  for(i in 1:length(data_list)){
    sampleName <- colnames(data_list[[i]])
    loading_m <- loading[sampleName,]
	basis <- qr.Q(qr = qr(x = loading_m))
    proj_m <- as.matrix(scaleDat[[i]]) %*% basis
    rownames(proj_m) <- var.names
    whitenData[[i]] <- proj_m
    
    ##correct data according to the shared latent space
    correctData <- data_list[[i]]@assays$RNA@scale.data[var.names,]
    rownames(correctData) <- var.names
    colnames(correctData) <- colnames(data_list[[i]])
    data_list[[i]][['CorrectExp']] <- CreateAssayObject(data = correctData)
    DefaultAssay(data_list[[i]]) <- "CorrectExp"
  }
  rm(projData,cor,scaleDat,cor.m,correctData)
  
  writeLines("Calculating basal components...")
  ica_pooling <- ICAcomputingCCA(whitenData,fastICA)
  ica.filter <- CrossBatchGrouping(ica_pooling,cor = corr,Unique.Preservation = Unique.Preservation,W.top=W.top)
  basalcomp <- BasalComponents(ica_pooling,ica.filter$cluster)
  
  res <- list()
  obj <- MergeSeurat(data_list)
  rownames(loading) <- colnames(obj)
  colnames(loading) <- paste0("lsi_",1:ncol(loading),sep="")
  obj[["LSI"]] <- CreateDimReducObject(embedding = loading,key="lsi_")
  
  res$mergedObject <- obj
  res$basalcomp <- basalcomp
  res
}

MergeSeurat <- function(data_list){
  obs_m <- NULL
  GEP <- NULL
  for(i in 1:length(data_list)){
    obs_m <- rbind(obs_m,data_list[[i]]@meta.data)
    data_m <- GetAssayData(data_list[[i]])
    colnames(data_m) <- paste0("Batch-",i,"-",colnames(data_m),sep="")
    GEP <- cbind(GEP,data_m)
  }
  obj <- CreateSeuratObject(GEP)
  rownames(obs_m) <- colnames(obj)
  obj@meta.data = obs_m
  obj
}

ICAcomputingCCA <- function(whitenData,fastICA = FALSE){
  ### computing ICA
  ICA <- NULL
  for(i in 1:length(whitenData)){
    if(!fastICA){
	require("ica")
	writeLines(paste("Using JADE implementation for batch-",i,"...",sep=""))
    ica.m <- icajade(whitenData[[i]],nc=ncol(whitenData[[i]]),center = FALSE)$S
    colnames(ica.m) <- paste0("Batch-",i,"_IC",1:ncol(ica.m),sep="")
    ICA <- cbind(ICA,ica.m)
  }else{
	require("fastICA")
    writeLines(paste("Using FastICA implementation for batch-",i,"...",sep=""))
    ica.m <- fastICA(whitenData[[i]],ncol(whitenData[[i]]), alg.typ = "parallel", fun = "logcosh")$S
	colnames(ica.m) <- paste0("Batch-",i,"_IC",1:ncol(ica.m),sep="")
	ICA <- cbind(ICA,ica.m)
  }
  }
  ICA
}

Min_Max <- function(vec){
  ql = min(vec)
  qu = max(vec)
  vec_t = (vec - ql)/(qu - ql)
  vec_t
}

BasalComponents <- function(ica_pooling,cluster){
  comp <- apply(ica_pooling,2,Min_Max)
  comp_sub <- NULL
  basal <- NULL
  for(i in names(table(cluster))){
     comp_sub <- as.matrix(comp[,cluster == i])
	 basal_comp <- rep(1,nrow(comp_sub))
	 for(j in 1:ncol(comp_sub)){
	    basal_comp <- basal_comp * comp_sub[,j]
     }
	 basal_comp <- basal_comp - mean(basal_comp)
	 basal <- cbind(basal, basal_comp)
   }
  basal
}

ReferenceIntegrate <- function(dataList,reduct.dim,var.names){
   ### using the overall variance to determine the reference
   var <- NULL
   for(i in 1:length(dataList)){
     svd <- rARPACK::svds(GetAssayData(dataList[[i]])[var.names,],k=reduct.dim)
	 var <- c(var,sum(svd$d^2))
   }
   order <- order(var,decreasing=TRUE)
   dataList <- dataList[order]
   
   writeLines(paste("Identify batch-",order[1]," as the reference",sep=""))
   dataList
}
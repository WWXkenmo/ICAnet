#' Title Cross Batch ICA Source Component Grouping
#' Group the expression programs (source components) generated from different batch of scRNA-seq dataset through PAM clustering.
#'
#' @param ica.pooling A matrix object that including all ICA source components generate from *ICAcomputing*,
#' @param k.max the maximum clusters number of ica-components (default: the number of ica-components - 1)
#' @param plot whether plot cross batch ica source component correlation map
#' @param cor the correlation method which determine the correlation between different pair of source components
#' @param W.top the threshold to determine the activated genes, the genes which has absolute attributes value large than threshold*standard derivation from mean are the activated genes (default: 2.5)
#' @param filtering whether filter the components which their 'activated gene' number lower than a specific values
#' @param threshold the threshold number to determine which component need to be filtered (default: 30)
#' @param Unique.Preservation whether to preserve the cluster with only one component included (default: TRUE)
#'
#' @return a list object which containing the filtered ica-components matrix.
#' @importFrom cluster pam
#' @importFrom pheatmap pheatmap
#'
#' @export
#'
CrossBatchGrouping <- function(ica.pooling,k.max=(ncol(ica.pooling)-1),plot=TRUE,cor="pearson",W.top=2.5,filtering=TRUE,threshold=30,Unique.Preservation=TRUE){

    raw.ica.pooling <- ica.pooling

    for(i in 1:ncol(ica.pooling)){
        vec<- ica.pooling[,i]
        vec[abs(vec)<W.top*var(vec)] <- 0
        ica.pooling[,i] <- vec
    }
    filter <- NULL; for(i in 1:nrow(ica.pooling)) filter <- c(filter,sum(ica.pooling[i,]!=0))
    ica.pooling <- ica.pooling[filter>0,]

    ModuleSize <- NULL;for(i in 1:ncol(ica.pooling)) ModuleSize <- c(ModuleSize,sum(ica.pooling[,i]!=0))
    if(filtering){ica.pooling <- ica.pooling[,ModuleSize>threshold]}

    if(plot){pheatmap(abs(cor(ica.pooling,method=cor)),main="Cross Batches Module Grouping")}

    distP.o <- as.dist( 1-abs(cor(ica.pooling,method=cor)));
    asw.v <- vector();
    for(k in 2:k.max){
        pam.o <- pam(distP.o,k,stand=FALSE);
        asw.v[k-1] <- pam.o$silinfo$avg.width
    }
    k.opt <- which.max(asw.v)+1;
    pam.o <- pam(distP.o,k=k.opt,stand=FALSE);

    ica.filter <- raw.ica.pooling[,pam.o$medoids]
    cat(paste("Identify",k.opt," patterns",sep=""))
    cluster <- pam.o$clustering

	if(!Unique.Preservation) ica.filter <- ica.filter[,as.numeric(names(table(cluster))[table(cluster)>1])]
    # return the result
    ica <- list();
    ica$ica.filter <- ica.filter
    ica$dist <- abs(cor(ica.pooling,method=cor))
    ica$cluster <- cluster
    ica
}

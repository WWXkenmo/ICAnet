#' Title Estimate the number of components
#' an Random Matrix theory based component number estimation
#'
#' @param data.m a matrix object which need to be estimated its number of components
#' @param default boolean variable to determine whether restrict the maximum number of decomposed components
#' @param svd.num an integer value to indicate maximum number of decomposed components
#'
#' @return a integer value which indicating the estimated number of ica components
#' @export
#' @importFrom coop tpcor
#' @importFrom rARPACK svds
#'
EstNumModule <- function(data.m, default=TRUE, svd.num=NULL){
    require(coop)
    require(rARPACK)

    data.m <- data.m[rowSums(data.m)>0,]

    M <- apply(data.m,2,function(X){ (X - mean(X))/sqrt(var(X))});
    M <- t(M)

    sigma2 <- var(as.vector(M));
    Q <- nrow(data.m)/ncol(data.m);
    ns <- ncol(data.m);
    lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
    lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
    C <- coop::pcor(t(M))

	if(default){
    if(ncol(data.m)>50){
        eigen.o <- svds(C,k=50)$d
    }else{
        eigen.o <- svd(C)$d
    }
    }else{
	    eigen.o <- svds(C,k=svd.num)$d
    }

    intdim <- length(which(eigen.o > lambdaMAX))
	if(default){
	if(intdim==50){
	cat('warning! the number of significant component is too large')
	}
	}else{
	if(intdim==svd.num){
	cat('warning! the number of significant component is too large')
	}
	}

    rm(M,data.m)

    #return
    intdim
}

#' Title Module Significance Test
#'
#' @param Data a Z-score transformed matrix
#' @param statl.v the ica-components
#' @param geneSets the module sets which need to be tested its significance
#' @param nMC the number of permutations
#' @param verbose a boolean variable, whether show the running process (default: FALSE)
#'
#' @return a integer value which indicating the significance of modules
#' @export
#'
#' @examples
ModuleSignificanceTest <- function(Data,statl.v,geneSets,nMC,verbose=FALSE){
    modN.v <- vector();
    for(v in 1:length(geneSets)){
        corM <- as.matrix(Data)[geneSets[[v]],]
        corM <- svd(corM)$d
        modN.v[v] <- max(corM)/min(corM)
    }

    if(verbose) print("Starting Monte Carlo Runs");
    ntop = length(geneSets);
    modNmc.m <- matrix(nrow=ntop,ncol=nMC);
    for(m in 1:ntop){
        for(run in 1:nMC){
            permN.idx <- sample(1:length(statl.v),length(geneSets[[m]]),replace=FALSE)
            corM <- svd(as.matrix(Data)[names(statl.v)[permN.idx],])$d
            Data_sampled <- max(corM)/min(corM)
            modNmc.m[m,run] <- Data_sampled;
        }
        if(verbose) print(paste("Done for seed/module ",m,sep=""));
    }

    modNpv.v <- rep(1,ntop);
    for(v in 1:ntop){
        modNpv.v[v] <- length(which(modNmc.m[v,] > modN.v[v]))/nMC;
    }

    modNpv.v
}

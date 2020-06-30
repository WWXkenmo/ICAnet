CrossBatchGrouping <- function(ica.pooling,k.max=(ncol(ica.pooling)-1),plot=TRUE,cor="pearson",W.top=2.5,filtering=TRUE,threshold=30){
    # ica.pooling: Independent component generate from previous steps
    # k.max: the maximum number of clusters
    # plot: whether to plot the component association matrix
    # cor: correlation methods
    # W.top: the threshold to define activated genes, so that to compute component correlation
    # filter: whether to filter the component with small number of activated genes
    # threshold: the minimum number of activated genes.
    
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
    
    if(plot){pheatmap(abs(cor(ica.pooling,method=cor)),main="Cross Batches Component Grouping")}
    
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
    # return the result
    ica <- list();
    ica$ica.filter <- ica.filter
    ica$dist <- abs(cor(ica.pooling,method=cor))
    ica$cluster <- cluster
    
    # Return basal programs
    ica
}
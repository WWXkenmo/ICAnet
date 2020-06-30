ICAcomputing <- function(obj, seurat.obj, batch, center=TRUE,scale=FALSE,RMT=TRUE, nc.vec=NULL, ICA.type="JADE", nbIt=10, alg.type="deflation",
         fun="logcosh",maxit =500, tol=10^-6, funClus="hclust",row.norm=FALSE,boostrap=FALSE,method="average",two.stage=TRUE){
    
    require("MineICA")
    require("ica")
    require("Seurat")
    require("rARPACK")
    
    #obj: dataset means gene-cell expression matrix from different batches.
    #seurat.obj: a seurat object
    #batch: a batch indicate vector
    #center: if your gene expression matrix need to be centered
    #scale:   ..... scaled
    #RMT: random matrix theory estimates
    #ICA.type: ICA type can be choosed as FastICA or JADE algorithm
    #The dataset is required to large than 30 samples
    
    ica.pooling <- NULL
    ind <- 1
    if(seurat.obj==TRUE){
        RMT.num <- NULL
        for(i in names(table(obj$batch))){
            print(paste("batch ",ind," Indepdent Component Analysis",sep=""))
            
            data <- as.matrix(GetAssayData(obj))[,obj$batch==i]
            if(ncol(data)>30){	
                print("emmm...scaling...")
                M <- t(apply(data,1,scale,scale=FALSE))
                print("Done Scaling")
                
                if(RMT){
                    print('Using RMT to estimate number of module')
                    num <- EstNumModule(data)
                    if(num == 0) num <- num+1;
                    print(paste('RMT estimate',num,'expression programm',sep=" "))
                    if(two.stage){
                        print("Running 2th-step of RMT...(Recommend for integrating different platform)")
                        res <- fastICA(M, n.comp = num);
                        mixture <- as.matrix(t(res$A))
                        num <- EstDimRMT(mixture)  
                        num <- num$dim
                        print(paste('Two step RMT estimate',num,'expression programm',sep=" "))
                    }
                }else{
                    num <- nc.vec[ind]
                }
                
                
                if(center){M <- t(apply(data,1,scale,scale=FALSE))}
                if(scale){M <- t(apply(data,1,scale,scale=TRUE))}
                
                if(ICA.type=="FastICA"){res <- clusterFastICARuns(X=M, nbComp=num, alg.type=alg.type, nbIt=nbIt,  funClus=funClus, method=method)}
                if(ICA.type=="JADE"){
                    if(num==1){res <- icajade(M, nc=num+1);
                    sc_data <- as.matrix(res$S[,1])
                    }else{
                        res <- icajade(M, nc=num);
                        sc_data <- as.matrix(res$S)
                    }
                }
                
                # name for result ica.pooling matrix
                id <- NULL;
                for(j in 1:num){ id <- c(id, paste(i,"-",j,sep=""))}
                colnames(sc_data) <- id
                ica.pooling <- cbind(ica.pooling, sc_data)
                
                paste("Done ",ind,"th batch",sep="")
                ind <- ind+1;
                RMT.num <- c(RMT.num,ncol(sc_data))
            }
        }
    }
    
    # return
    res <- list()
    res$ica.pooling <- ica.pooling
    res$RMT.num <- RMT.num
    res
}

ICAcomputing <- function(obj, seurat.obj, batch, center=TRUE,scale=FALSE,RMT=TRUE, nc.vector=NULL, ICA.type="FastICA", nbIt=10, alg.type="deflation",
                   fun="logcosh",maxit =500, tol=10^-6, funClus="hclust",row.norm=FALSE,boostrap=FALSE,method="average"){
	
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
	#nc.vector: number of components need to be estimate of each dataset
	#ICA.type: FastICA or JADE
	
    ica.pooling <- NULL
	ind <- 1
	if(seurat.obj==TRUE){
	for(i in names(table(obj$batch))){
	print(paste("batch ",ind," Indepdent Component Analysis",sep=""))
	
	data <- as.matrix(GetAssayData(obj))[,obj$batch==i] 
	print("emmm...scaling...")
	M <- t(apply(data,1,scale,scale=FALSE))
	print("Done Scaling")
	
	if(RMT){
	num <- EstNumModule(data)
	
	print(paste('RMT estimate',num,'expression programm',sep=" "))
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
	}
	}
	
	# return
	res <- list()
	res$ica.pooling <- ica.pooling
	res
}

crossBatchMapping <- function(ica.pooling,k.max=(ncol(ica.pooling)-1),plot=TRUE,cor="pearson",W.top=2.5,filtering=TRUE,threshold=30){
    require(pheatmap)
	
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
	
	if(plot){pheatmap(cor(ica.pooling,method=cor),main="Cross Batches Module Mapping")}

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
	
	# return the result
	ica.filter
}



RunICAnet <- function(obj,ica.filter, W.top=2.5,PPI.net=NULL, PPI="String",species=9606,score=600,max.step=10,small.size=3,
                        aucMaxRank=3000,cores=6){
    
	PPI.network <- getPPI_String(obj, species = species, score_threshold = score)
    netGenes <- intersect(rownames(PPI.network),rownames(ica.filter))
    hs_network_matrix <- PPI.network[netGenes,netGenes]
	geneSets <- list()
    for(c in 1:ncol(ica.filter)){
       ica.score <- ica.filter[,c]
       # iter select significant gene signitures
  
       ica.score[abs(ica.score)<W.top*var(ica.score)] <- 0
       ica.score <- ica.score[netGenes]
  
    if(sum(ica.score!=0)>0){
       #filter
       ica.score <- ica.score
       ica.score <- abs(ica.score)
    
       m <- rep(as.matrix(ica.score),length(ica.score))
       m <- matrix(m,nrow=length(ica.score),ncol=length(ica.score))
       cor.m <- m+t(m)
       cor.m <- cor.m/(2*max(ica.score))
    
       ##trimming the correlation network
       cor.m.filter <- as.matrix(ica.score) %*% t(as.matrix(ica.score))
       cor.m.filter[cor.m.filter!=0] <- 1
       cor.m <- cor.m.filter * cor.m
       rm(cor.m.filter)
    
       Data_network_final <- cor.m * as.matrix(hs_network_matrix)
       rm(cor.m,m)
       filtered_row <- rowSums(as.matrix(Data_network_final)) > 0
       Data_network_final <- Data_network_final[filtered_row,filtered_row]
    
       network_trim <- igraph::graph_from_adjacency_matrix(Data_network_final,weighted = TRUE,mode="undirected")
       ## network decomposition
       gene_sets_all <- list()
       run_label = c()
       for (rw_step in 1:max.step){
          set.seed(74156)
          network_cluster <- igraph::walktrap.community(network_trim,steps=rw_step)
          gene_sets_all <- c(gene_sets_all,communities(network_cluster))
          run_label <- c(run_label,rep(rw_step,length(communities(network_cluster))))
       }
    
    ######
    temp_len <- c()
    k = 0;
    for (i in 1:length(gene_sets_all)){
      temp <- as.character(unlist(gene_sets_all[i]))
      if(length(temp)>=small.size){
        k <- k+1;
        temp_len[k] <- length(temp)
        geneSets[[paste('ICA.net',c,run_label[i],k,temp_len[k],sep = ".")]] <- temp
      }
    }
    ##############
  }
}
    
	cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(obj)), nCores = cores)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=aucMaxRank, nCores =cores)
    cells_AUC_matrix <- getAUC(cells_AUC)
    
    # loading subnetwork activity matrix
	obj[["IcaNet"]] <- CreateAssayObject(counts = as.matrix(cells_AUC_matrix))
    DefaultAssay(obj) <- "IcaNet"
    obj<- ScaleData(obj, assay = "IcaNet")
	
	# return seurat obj
	obj
}





EstNumModule <- function(data.m){
  require(coop)
  require(rARPACK)
  
  data.m <- data.m[rowSums(data.m)>0,]
  data.m <- t(data.m)

  M <- apply(data.m,2,function(X){ (X - mean(X))/sqrt(var(X))});
 
  sigma2 <- var(as.vector(M));
  Q <- nrow(data.m)/ncol(data.m);
  ns <- ncol(data.m);
  lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
  lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
  C <- coop::tpcor(M)
  eigen.o <- svds(C,k=30)$d
 
  intdim <- length(which(eigen.o > lambdaMAX))
  rm(M,data.m)

  #return 
  intdim
}

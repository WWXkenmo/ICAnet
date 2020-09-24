#' Title Using independent components to construct weighted molecular network to detect function module for single cell clustering
#' ICAnet used independent components to construct weighted PPI, and running walk-trap algorthm to detect module on it. The resulted modules are used for the following analysis
#'
#' @param obj a Seurat object
#' @param ica.filter the filtered/unfiltered ica-components set
#' @param W.top the threshold to determine the activated genes, the genes which has absolute attributes value large than threshold*standard derivation from mean are the activated genes (default: 2.5)
#' @param PPI.net a matrix object which indicating the PPI network, a boolean network is required
#' @param max.step Integer number. The maximum step run in the walk-trapped based community detect.
#' @param small.size integer number to determine the minimum size of module. The module which has the number of gene member less than this value will be filtered
#' @param nMC the number of permutations, which is used for calculate the pvalue of each module (default: 100)
#' @param aucMaxRank Integer number. The number of highly-expressed genes to include when computing AUCell
#' @param cores the number of cores used for computation
#' @param ModuleSignificance the boolean variable to indicate whether perform module significant test (default: FALSE)
#' @param scale the boolean variable to indicate whether perform scaling over each batch of scRNA-seq gene expression data
#'
#' @return a Seurat object which contain the "IcaNet" assay
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph walktrap.community
#' @importFrom igraph communities
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat Misc
#' @importFrom AUCell AUCell_buildRankings
#' @importFrom AUCell AUCell_calcAUC
#' @importFrom AUCell getAUC
#'
RunICAnet <- function(obj,ica.filter, W.top=2.5,PPI.net=NULL, species=9606,score=600,max.step=10,small.size=3,
                            nMC=100,aucMaxRank=3000,cores=6,ModuleSignificance=TRUE,scale=TRUE){

    PPI.network <- PPI.net

    netGenes <- intersect(rownames(PPI.network),rownames(ica.filter))
    hs_network_matrix <- PPI.network[netGenes,netGenes]
    geneSets <- list()
    if(ModuleSignificance) modN.score <- NULL
    for(c in 1:ncol(ica.filter)){
        ica.score <- ica.filter[,c]
        # iter select significant gene signitures

        ica.score[abs(ica.score)<W.top*sqrt(var(ica.score))] <- 0
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

			###rescaling
			if(ModuleSignificance){
			Scale.data <- NULL
			for(i in names(table(obj$batch))){
			Data <- as.matrix(GetAssayData(obj))[,obj$batch %in% i]
			statl.v <- abs(ica.score[ica.score!=0])
			Data <- Data[names(statl.v),]
	        Data <- t(apply(Data,1, function(x) {(x - mean(x)) / sqrt(var(x))}))
			Data[is.nan(Data)] <- 0
			Scale.data <- cbind(Scale.data,Data)
			}
			}
            temp_len <- c()
            k = 0;
            for (i in 1:length(gene_sets_all)){
                temp <- as.character(unlist(gene_sets_all[i]))
                if(length(temp)>=small.size){
                    k <- k+1;
                    temp_len[k] <- length(temp)
					if(ModuleSignificance){
					score <- ModuleSignificanceTest(Scale.data,statl.v,list(temp),nMC,verbose=FALSE)
					modN.score <- c(modN.score, score)
					}
                    geneSets[[paste('ICAnet',c,k,temp_len[k],sep = "-")]] <- temp
                }
            }
            ##############
        }
    }

	if(scale){
	cat("\nRunning Scaling on each batch(Improve Batch Corrections)")
	scale.data <- NULL
	for(i in names(table(obj$batch))){
	   data <- obj[,obj$batch %in% i]
	   data <- ScaleData(data);
	   eval(parse(text = paste("data <- data@assays$",DefaultAssay(data),"@scale.data",sep="")))
	   scale.data <- cbind(scale.data, data)
	   cat(paste("\nDone Batch ",i,"\n",sep=""))
	}
	obj[['scale.data']] <- CreateAssayObject(scale.data)
	DefaultAssay(obj) <- "scale.data"
	}

    final_geneSets <- geneSets[!duplicated(geneSets)];print(length(geneSets));if(ModuleSignificance) print(length(modN.score))
    if(ModuleSignificance) names(modN.score) <- names(geneSets)
    if(ModuleSignificance) modN.score <- modN.score[names(final_geneSets)]
    print(paste('num:',length(final_geneSets),sep=""))
    cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(obj)), nCores = cores)
    cells_AUC <- AUCell_calcAUC(final_geneSets, cells_rankings, aucMaxRank=aucMaxRank, nCores =cores)
    cells_AUC_matrix <- getAUC(cells_AUC)

    # loading subnetwork activity matrix
    obj[["IcaNet"]] <- CreateAssayObject(counts = as.matrix(cells_AUC_matrix))
    DefaultAssay(obj) <- "IcaNet"
    obj<- ScaleData(obj, assay = "IcaNet")
    Misc(obj, slot = 'IcaNet_geneSets') <- final_geneSets
    if(ModuleSignificance) Misc(obj, slot = 'modN.score') <- modN.score
    # return seurat obj
    obj
}

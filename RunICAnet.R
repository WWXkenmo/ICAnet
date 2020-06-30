RunICAnet <- function(obj,ica.filter, W.top=2.5,PPI.network=NULL, max.step=10,small.size=3,aucMaxRank=3000,cores=6){
    #obj: Seurat object
    #ica.filter: the basal programs or independent components.
    #W.top: the threshold to define activated genes
    #PPI.network: the protein-protein interaction network
    #max.step: the maximum steps to perform random walk with trapping clustering
    #small.size: the minimum gene number for a module
    
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
    final_geneSets <- NULL
    geneSets <- unique(geneSets)
    for(i in 1:length(geneSets)){
        final_geneSets[[paste('Ica.net',i,sep="")]] <- geneSets[[i]]
    }
    print(paste('num:',length(final_geneSets),sep=""))
    cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(obj)), nCores = cores)
    cells_AUC <- AUCell_calcAUC(final_geneSets, cells_rankings, aucMaxRank=aucMaxRank, nCores =cores)
    cells_AUC_matrix <- getAUC(cells_AUC)
    
    # loading subnetwork activity matrix
    obj[["IcaNet"]] <- CreateAssayObject(counts = as.matrix(cells_AUC_matrix))
    DefaultAssay(obj) <- "IcaNet"
    obj<- ScaleData(obj, assay = "IcaNet")
    Misc(obj, slot = 'IcaNet_geneSets') <- final_geneSets
    # return seurat obj
    obj
}
#' Title Convert feather file into motif-gene boolean network
#' This function is provided for convert feather file into motif-gene boolean network, the feather file containing the motif binding informtion for human/mouse/fly species, the files could be downloaded from cisTarget databases https://resources.aertslab.org/cistarget/
#'
#' @param feather.readPath a character value to indicate the path which saved the feather file
#' @param cutoff a cutoff to determine how to filter the unwanted Motif-gene links (default: 1)
#'
#' @return a matrix object which containing the boolean network of TF binding motif-gene
#' @export
#' @importFrom  RcisTarget importRanking
#' @examples
TF_Net_Generate <- function(feather.readPath,cutoff=1){
cat('Loading TF,motif annotation dataset...')
motif_ranking <- importRankings(feather.readPath)
motif_ranking <- motif_ranking@rankings
motif_ranking <- as.data.frame(motif_ranking)
rownames(motif_ranking) <- motif_ranking$features
motif_ranking <- motif_ranking[,-1]
motif_ranking <- as.matrix(motif_ranking)

###Apply Zscore transform to the matrix
cat('\nTrimming the interaction...')
motif_ranking_Zscore <- t(apply(motif_ranking,1,scale))
colnames(motif_ranking_Zscore) <- colnames(motif_ranking)
## Using Zscore to identify the significant binding targets
motif_ranking_Zscore[motif_ranking_Zscore>cutoff] <- 1
motif_ranking_Zscore[motif_ranking_Zscore!=1] <- 0


###Apply Zscore transform to the matrix
motif_ranking_Zscore_col <- (apply(motif_ranking,2,scale))
rownames(motif_ranking_Zscore_col) <- rownames(motif_ranking)
## Using Zscore to identify the significant binding targets
motif_ranking_Zscore_col[motif_ranking_Zscore_col>cutoff] <- 1
motif_ranking_Zscore_col[motif_ranking_Zscore_col!=1] <- 0

###Summarize
cat('\nGenerating TF-Gene-Net...')
motif_gene_net <- motif_ranking_Zscore * motif_ranking_Zscore_col
###Convert to TF-Gene Net

cat('\nDone')

motif_gene_net
}

#' Title Using independent components to perform TF-target enrichment to detect TF-regulon for single cell clustering
#' ICAnet used independent components to perform TF-target enrichment test, and select significant TF-regulons for the following analysis
#'
#' @param obj a Seurat object
#' @param ica.filter the filtered/unfiltered ica-components set
#' @param W.top.TFs the threshold to determine the activated TFs, the TFs which has absolute attributes value large than threshold*standard derivation from mean are the activated TFs (default: 2.5)
#' @param W.top.genes the threshold to determine the activated genes, the genes which has absolute attributes value large than threshold*standard derivation from mean are the activated genes (default: 2.5)
#' @param Motif_Net the matrix which indicating the boolean network of motif-target
#' @param TF_motif_annot Annotation to human/mouse/fly transcripton factors for the motifs in each motif collection (e.g. mc8nr or mc9nr)
#' @param small.size integer number to determine the minimum size of module. The module which has the number of gene member less than this value will be filtered
#' @param aucMaxRank Integer number. The number of highly-expressed genes to include when computing AUCell
#' @param cores the number of cores used for computation
#' @param verbose a boolean variable, whether show the running process (default: FALSE)
#' @param nMC the number of permutations, which is used for calculate the pvalue of each module (default: 100)
#' @param cutoff the significant level to filter the TF-regulons (default: 0.01)
#'
#' @return a Seurat object which contain the "IcaNet_TF" assay
#' @export
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat Misc
#' @importFrom AUCell AUCell_buildRankings
#' @importFrom AUCell AUCell_calcAUC
#' @importFrom AUCell getAUC
#'
#' @examples
RunICAnetTF <- function(obj,ica.filter, W.top.TFs=2.5, W.top.genes=2.5, Motif_Net=NULL, TF_motif_annot=NULL,small.size=3,
                          aucMaxRank=3000,cores=6,verbose=FALSE,nMC=100,cutoff=0.01){

        cat('Running on TF-Gene Network')
        netGenes <- intersect(colnames(Motif_Net),rownames(ica.filter))
        hs_network_matrix <- Motif_Net[,netGenes]
        geneSets <- list()
        TF_activity_all <- moduleNames <- significance <- modularity <- NULL
        motif_list <- intersect(names(table(TF_motif_annot$motif)),rownames(Motif_Net))
        TF_motif_annot <- subset(TF_motif_annot,TF_motif_annot$motif %in% motif_list == TRUE)
        Motif_Net <- Motif_Net[motif_list,]
        TF_list <- names(table(TF_motif_annot$TF))

        for(c in 1:ncol(ica.filter)){
            ica.score.genes <- ica.score.TFs <- ica.filter[,c]
            ica.score.genes[abs(ica.score.genes)<W.top.genes*sqrt(var(ica.score.genes))] <- 10^-3
            ica.score.TFs[abs(ica.score.TFs)<W.top.TFs*sqrt(var(ica.score.TFs))] <- 10^-3
            ica.score.genes <- ica.score.genes[netGenes]
            ica.score.TFs <- ica.score.TFs[netGenes]

            seeds <- intersect(names(ica.score.TFs)[ica.score.TFs!=10^-3], TF_list)
            activator_genes <- names(ica.score.genes)[ica.score.genes!=10^-3&ica.score.genes>0]
            repressor_genes <- names(ica.score.genes)[ica.score.genes!=10^-3&ica.score.genes<0]

            #if(sum(ica.score!=0)>0){
            cat(paste('\nRunning on the ',c,'st component',sep=""))
            ica.score <- ica.filter[,c]
            ica.score <- (ica.score[netGenes])


            if(length(seeds)!=0){
                cat(paste('\nIdentify ',length(seeds),' activated TF',sep=""))
                gene_sets_all <- list();
                ind <- 1
                id <- NULL
                TF_activity <- NULL
                for(v in 1:length(seeds)){
                    seeds.TF <- subset(TF_motif_annot,TF_motif_annot$TF %in% seeds[v] == TRUE)$motif
                    if(length(seeds.TF)>=1){
                        for(k in 1:length(seeds.TF)){
                            gene_sets_all[[ind]] <- colnames(hs_network_matrix)[hs_network_matrix[seeds.TF[k],]!=0]
                            ind <- ind + 1
                            id <- c(id, paste(seeds[v],"+",seeds.TF[k],sep=""))
                            TF_activity <- c(TF_activity, ica.score[seeds[v]])
                        }
                    }
                }
                names(gene_sets_all)<- id;

                cat('\nperform modularity test')
                modN.v <- vector();
                for(v in 1:length(gene_sets_all)){
                    modN.v[v] <- mean(ica.score[gene_sets_all[[v]]])
                }
                names(modN.v) <- id;

                print("Starting Monte Carlo Runs");
                ntop <- length(gene_sets_all)
                modNmc.m <- matrix(0,nrow=ntop,ncol=nMC);

                for(m in 1:ntop){
                    permN.idx <- array(0,dim=c(length(gene_sets_all[[m]]),nMC))
                    for(run in 1:nMC){
                        permN.idx[,run] <- sample(1:length(ica.score),length(gene_sets_all[[m]]),replace=FALSE);
                    }
                    for(run in 1:nMC){
                        tmpEM.v <- mean(ica.score[permN.idx[,run]])
                        modNmc.m[m,run] <- tmpEM.v;
                    }
                    if(verbose) print(paste("Done for",id[m],sep=""));
                }

                modNpv.v <- rep(1,ntop);
                for(v in 1:ntop){
                    modNpv.v[v] <- length(which(abs(modNmc.m[v,]) > abs(modN.v[v])))/nMC;
                }

                temp_len <- c()
                k <- 0
                for (i in which(modNpv.v<=cutoff)){
                    temp <- as.character(unlist(gene_sets_all[i]))
                    if(modN.v[i]>0) temp <- intersect(temp,activator_genes)
                    if(modN.v[i]<0) temp <- intersect(temp,repressor_genes)
                    if(length(temp)>=small.size){
                        k <- k+1;
                        temp_len[k] <- length(temp)
                        geneSets[[paste('ICAnet',c,temp_len[k],names(gene_sets_all)[i],sep = "-")]] <- temp
                        moduleNames <- c(moduleNames, names(gene_sets_all)[i])
                        modularity <- c(modularity, modN.v[i])
                        significance <- c(significance, modNpv.v[i])
                        TF_activity_all <- c(TF_activity_all,TF_activity[i])
                    }
                }
            }
        }

    metadata <- cbind(modularity,significance,TF_activity_all)
    rownames(metadata) <- names(geneSets)
    final_geneSets <- geneSets[!duplicated(geneSets)]
    metadata <- metadata[names(final_geneSets),]
    print(paste('num:',length(final_geneSets),sep=""))
    cells_rankings <- AUCell_buildRankings(as.matrix(GetAssayData(obj)), nCores = cores)
    cells_AUC <- AUCell_calcAUC(final_geneSets, cells_rankings, aucMaxRank=aucMaxRank, nCores =cores)
    cells_AUC_matrix <- getAUC(cells_AUC)
    rownames(cells_AUC_matrix) <- names(final_geneSets)
	######prcessing module informations
	module_infor <- strsplit(names(final_geneSets),split="[+]")
	module_gene <- motif_list <- NULL; for(i in 1:length(module_infor)){module_gene <- c(module_gene, module_infor[[i]][1]); motif_list <- c(motif_list, module_infor[[i]][2])}
	module_gene <- strsplit(module_gene,split="[-]")
	gene <- moduleSize <- NULL; for(i in 1:length(module_gene)){gene <- c(gene, module_gene[[i]][4]);moduleSize <- c(moduleSize, module_gene[[i]][3])}
	metadata <- cbind(gene, moduleSize, motif_list, metadata)
	significance <- as.numeric(as.matrix(metadata[,5]))

    # loading subnetwork activity matrix
    obj[["IcaNet_TF"]] <- CreateAssayObject(counts = as.matrix(cells_AUC_matrix))
    DefaultAssay(obj) <- "IcaNet_TF"
	names(significance) <- rownames(obj)

    obj<- ScaleData(obj, assay = "IcaNet_TF")
    Misc(obj, slot = 'IcaNet_geneSets_TF') <- final_geneSets
	Misc(obj, slot = 'modN.score') <- significance
    Misc(obj, slot = 'IcaNet_geneSets_TF_moduleInfor') <- metadata
    # return seurat obj
    obj
}

#' @title Using Independent Component Analysis to decompose scRNA-seq gene expression matrix
#'
#' @description We introduced independent component analysis (ICA) into single cell clustering to decompose the gene expression matrix into a number of independent components. Each independent component was characterized by a co-expression pattern and was associated with certain meaningful biological pathway.
#'
#' @param obj a seurat object which contain the gene expression matrix of all batches, meanwhile it also contain the *batch* ID which document the batch label information
#' @param center a boolean variable to determine whether perform centralization to each batch
#' @param scale a boolean variable to determine whether perform scaling to each batch
#' @param RMT a boolean variable to determine whether perform random matrix theory based number of components estimation
#' @param nc.vec a numrical vector to determine the number of independent components for each batch
#' @param nc.global.mode a numeric vector to determine the number of components of all batches when running global mode.
#' @param global.mode a boolean variable to determine whether running global mode based ICA
#' @param ICA.type an optional character string given a method for computing independent components. One of "JADE" (default), "FastICA", or "MineICA": can be abbreviated
#' @param nbIt The number of iterations of FastICA in MineICA
#' @param funClus The clustering function to be used to cluster the estimates in MineICA. One of "hclust"(default),"kmeans", "pam", "agnes" can be abbreviated
#' @param two.stage Whether running ICA+RMT based components number estimation
#' @param two.stage.ICA Which type of ICA methods shall be applied in ICA+RMT estimate. One can select "FastICA" or "JADE"
#'
#' @return a list object which containing the ica-components
#' @export
#' @importFrom ica icajade
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat ScaleData
#' @importFrom fastICA fastICA
#' @importFrom rARPACK svds
#' @importFrom MineICA clusterFastICARuns
#' @importFrom isva EstDimRMT
ICAcomputing <- function(obj, center=TRUE,scale=FALSE,RMT=TRUE, nc.vec=NULL, nc.global.mode=NULL,global.mode=FALSE,ICA.type="JADE", nbIt=10, funClus="hclust",two.stage=TRUE,two.stage.ICA="JADE"){

    #obj: the single cell gene expression datasets
    #seurat.obj: if the "obj" is a Seurat object
    #batch: if obj is a matrix object, then the a batch indicate vector need to be defined
    #center: if your gene expression matrix need to be centered
    #scale:   ..... scaled
    #RMT: random matrix theory estimates
    #ICA.type: ICA type can be choosed as FastICA, MineICA or JADE algorithm
    ica.pooling <- NULL
    ind <- 1
	ica.m.pooling <- NULL
	RMT.num <- NULL
	if(global.mode){
		    print("Running scaling on each batch")
			M <- NULL
			for(i in names(table(obj$batch))){
                data <- as.matrix(GetAssayData(obj))[,obj$batch==i]
				print(paste('Centering on ',i,sep=""))
			    M <- cbind(M,t(apply(data,1,scale,scale=FALSE)))
		    }

			print("Running ICA...")
			if(ICA.type=="FastICA"){
			X.ICA <-fastICA(M,n.comp=nc.global.mode)
			sc_data <- as.matrix(X.ICA$S)
			sc_data_M <- as.matrix(X.ICA$A)
			}
			if(ICA.type=="JADE"){
			res <- icajade(M, nc=nc.global.mode)
			sc_data <- as.matrix(res$S)
			sc_data_M <- as.matrix(res$M)
			}
			ica.pooling <- sc_data
            ica.m.pooling <- sc_data_M

			}else{
            for(i in names(table(obj$batch))){
            print(paste("batch ",ind," Indepdent Component Analysis",sep=""))

            data <- obj[,obj$batch==i]
            if(ncol(data)>30){
                if(center){
				print("emmm...centering...")
                M <- ScaleData(data,do.scale=FALSE,features=rownames(data))
				eval(parse(text = paste("M <- M@assays$",DefaultAssay(M),"@scale.data",sep="")))
                print("Done Centering")
				}
                if(scale){
				print("emmm...scaling...")
                M <- ScaleData(data,do.scale=TRUE,features=rownames(data))
				eval(parse(text = paste("M <- M@assays$",DefaultAssay(M),"@scale.data",sep="")))
                print("Done scaling")
				}

                if(RMT){
				    print('Using RMT to estimate number of module')
					if(center) num <- EstNumModule(as.matrix(GetAssayData(data)))
					if(scale) num <- EstNumModule(M)
                    if(num == 0) num <- num+1;
                    print(paste('RMT estimate',num,'expression programm',sep=" "))
					if(two.stage){
					  print("Running 2th-step of RMT...")
					  if(two.stage.ICA=="FastICA"){
					  res <- fastICA(M, num);
					  mixture <- as.matrix(res$A)
					  num <- EstDimRMT(mixture)
					  num <- num$dim}
					  if(two.stage.ICA=="JADE"){
					  res <- icajade(M, num);
					  mixture <- as.matrix(res$M)
					  num <- EstDimRMT(mixture)
					  num <- num$dim}
					  print(paste('Two step RMT estimate',num,'expression programm',sep=" "))
					  }
                }else{
                    num <- nc.vec[ind]
                }


                if(ICA.type=="MineICA"){
				    res <- clusterFastICARuns(X=M, nbComp=num, nbIt=nbIt,  funClus=funClus, method="average")
				}
				if(ICA.type=="FastICA"){
				 X.ICA <-fastICA(M,n.comp=num)
			     sc_data <- as.matrix(X.ICA$S)
				}
                if(ICA.type=="JADE"){
                    if(num==1){res <- icajade(M, nc=num+1);
                    sc_data <- as.matrix(res$S[,1])
                    }else{
                        res <- icajade(M, nc=num);
                        sc_data <- as.matrix(res$S)
                    }
                }
			}else{
			print('Warning! the number of cells is too low')
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


    # return
    res <- list()
    res$ica.pooling <- ica.pooling
	res$ica.pooling.m <- ica.m.pooling
    res$RMT.num <- RMT.num
    res
}

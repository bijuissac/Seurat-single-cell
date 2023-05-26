#source required libraries
print("Loading Required Libraries and Functions")
print(paste0("Snapshot Folder: ",r_library))
.libPaths(r_library)
r_packages <- c("devtools","Seurat","SeuratDisk","tidyverse","Matrix","DT","webshot",
		"htmlwidgets","pbapply","pheatmap","CellChat","SoupX","harmony",
		"SingleR","celldex","enrichR","data.table","clusterProfiler",
		"enrichR","enrichplot","scds","velocyto.R","velociraptor","SingleCellExperiment")

installed_packages <- library()$results[,1]
missing_packages <- r_packages[c(grep(FALSE,r_packages %in% installed_packages))]

print(paste0("Following R packages are not in the R library :",missing_packages))
 
available_packages <- r_packages[! r_packages %in% missing_packages]

for(pkg in available_packages){
	library(pkg, character.only = TRUE)
}

#loading OrgDb
#print(OrgDb)
#library(OrgDb,character.only=TRUE)

#function to load input data
readInData <- function(runSoupX,data_dir,sample_name,inputFormat){
	if(runSoupX){
		raw_obj <- load10X(dataDir = paste0(data_dir,"/",sample_name,"/outs/"))
		raw_obj <- autoEstCont(raw_obj)
		sample.obj <- adjustCounts(raw_obj)
	}else
	{
		if(inputFormat == "MEX"){
			sample.obj <- Read10X(data.dir = paste0(data_dir,"/",sample_name,"/outs/filtered_feature_bc_matrix/"))
		} else if (inputFormat == "H5"){
			sample.obj <- Read10X_H5(filename = paste0(data_dir,"/",sample_name,"/outs/filtered_feature_bc_matrix.h5"))
		}
	}
	print(paste0(sample_name," data loaded successfully"))
	print("Creating Seurat Object if atleast 3 cells and 200 features")
	mySO <- CreateSeuratObject(counts = sample.obj, project = sample_name, min.cells = 3, min.features = 200)
	return(mySO)
}

#function to QC sample
qcSample <- function(qc_filtering_method, outputFolder, mySO, sample_name, mt.threshold, mtPattern, riboPattern){
	
	mt.threshold <- as.numeric(mt.threshold)

	#Raw cell numbers
	print("Raw Cell Numbers")
	print(table(mySO$orig.ident))

	#calculate percentage mitochondria
	mySO[["percent_mt"]] <- PercentageFeatureSet(object=mySO, pattern = mtPattern)

	#calculate ribosomal protein proportions
	mySO[["percent_ribo"]] <- colSums(mySO[grepl(riboPattern, rownames(mySO), ignore.case = TRUE),])/colSums(mySO)

	#Identify doublets
	whichDoublets <- idDoublets(mySO)

	#generate before and after QC plots
	plotFile <- paste0(outputFolder,"/plots/qcPlots/",sample_name,"_raw.png")
	plotQC <- qcPlot(mySO@assays$RNA@counts, mtPattern)
	plotQC <- plotQC + labs(caption = paste0('Doublets = ', sum(whichDoublets)))
	
	png(plotFile, width = 1500, height = 1000, res = 300)
	print(plotQC)
	dev.off()

	#QC sample based on selection
	if(qc_filtering_method == "percentile"){
		#calculate percentile
		my.percentile <- quantile(mySO@meta.data$percent_mt, probs = 0.95)
		mySO <- subset(x = mySO, subset = percent_mt < my.percentile[[1]])
		print("Percentile threshold")
		print(my.percentile)
		print("Percentile Filtered Cell Numbers")
		print(table(mySO@meta.data$orig.ident))
	}
	if(qc_filtering_method == "Seurat"){
		mySO.filtered <- subset(x = mySO, subset = ((nFeature_RNA > min_features) & (nFeature_RNA < max_features) & (percent_mt <- (mt.threshold*100)) & (percent_ribo < ribo_fraction)))
	}else
	{
		mySO.filtered <- myQC(mySO, mt.threshold)
	}	

	print("Filtered Cell Numbers")
	print(table(mySO.filtered$orig.ident))

        plotFile <- paste0(outputFolder,"/plots/qcPlots/",sample_name,"_filtered.png")
        plotQC <- qcPlot(mySO.filtered@assays$RNA@counts, mtPattern)
        plotQC <- plotQC

        png(plotFile, width = 1500, height = 1000, res = 300)
        print(plotQC) 
        dev.off()

	return(mySO.filtered)

}

#Identify Doublets
idDoublets <- function(mySO){
	num_cells <- ncol(mySO)
	barCodes <- colnames(mySO)
	mySO <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mySO@assays$RNA@counts))
	mySO <- cxds(mySO)
	mySO <- mySO$csds_score
	idDoublets <- rep(FALSE, num_cells)
	names(idDoublets) <- barCodes
	idDoublets[names(boxplot.stats(mySO)$out)] <- TRUE
	return(idDoublets)
}

#Generate Plots for QC
qcPlot <- function(mySO, mtPattern){
	librarySize <- colSums(mySO)
	num_features <- colSums(mySO!=0)
	mito_proportion <- colSums(mySO[grepl(mtPattern, rownames(mySO), ignore.case = TRUE),])/colSums(mySO)
	qc.df <- data.frame(librarySize, num_features, mito_proportion)

	QC1 <- ggplot(qc.df, aes('', librarySize)) +
		geom_violin() + theme_light() +
		geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
		xlab('') + ylab('Library Size') +
		labs(subtitle = parse(text = paste0('italic(n) ==', ncol(mySO))))

        QC2 <- ggplot(qc.df, aes('', num_features)) +
                geom_violin() + theme_light() +
                geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
                xlab('') + ylab('nFeatures') 

        QC3 <- ggplot(qc.df, aes('', mito_proportion)) +
                geom_violin() + theme_light() +
                geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
                xlab('') + ylab('mtProportion')
	
	plot_obj <- QC1 + QC2 + QC3
	return(plot_obj) 
}


#Generate TSNE and UMAP plots
plot_tsne_umap <- function(outputFolder, mySO, myPrefix, reduction_method){
	reductionFile <- paste0(outputFolder,"/plots/clustering/",myPrefix,"_",reduction_method,".png")
	png(reductionFile,res = 300, width = 2000, height = 2000)
	if(reduction_method == "umap"){
		reducPlot <- UMAPPlot(X, label = TRUE, repel = TRUE) + theme_classic() + 
		theme(panel.grid = element_blank(), panel.border = element_blank()) + xlab('UMAP 1') + ylab('UMAP 2')
	}else if(reduction_method == "tsne"){
		reducPlot <- TSNEPlot(X, label = TRUE, repel = TRUE) + theme_classic() + 
		theme(panel.grid = element_blank(), panel.border = element_blank()) + xlab('TSNE 1') + ylab('TSNE 2')
	}

	print(reducPlot)
	dev.off()
}

#Integrate Data
integrate_data <- function(mySO, mySO.list, integration_method, ncps){

	seed_val <- 1
	if(integration_method == "seurat"){
		my.anchors <- FindIntegrationAnchors(object.list = mySO.list, dims = 1:20)
		mySO <- IntegrateData(anchorset = my.anchors, dims = 1:20)
		DefaultAssay(mySO) <- "integrated"
		mySO <- RunPCA(mySO, verbose = T, features = VariableFeatures(object = mySO), npcs = npcs)
		mySO <- FindNeighbors(mySO, reduction = "pca", dims = 1:20)
		mySO <- FindClusters(mySO, resolution = "0.8", random.seed = seed_val)
		mySO <- RunTSNE(mySO, reduction = "pca", dims = 1:20, seed.use = seed_val)
		mySO <- RunUMAP(mySO, reduction = "pca", dims = 1:20, seed.use = seed_val)
	}else if(integration_method == "harmony"){
		mySO <- RunHarmony(mySO, group.by.vars = 'orig.ident', max.iter.harmony = 1e3)
		mySO <- RunPCA(mySO, verbose = T, features = VariableFeatures(object = mySO), npcs = npcs)
                mySO <- FindNeighbors(mySO, reduction = "harmony", dims = 1:20)
		mySO <- FindClusters(mySO, resolution = "0.8", random.seed = seed_val)
                mySO <- RunTSNE(mySO, reduction = "harmony", dims = 1:20, seed.use = seed_val)
		mySO <- RunUMAP(mySO, reduction = "harmony", dims = 1:20, seed.use = seed_val)
	}
	return(mySO)

}

#Dynamically QC and filter sample
myQC <- function(mySO, mt.threshold, minLSize  = 1000){
	
	if(class(mySO) == 'Seurat'){
		countMatrix <- mySO@assays$RNA@counts
  	} else {
		countMatrix <- mySO
	}
	librarySize <- colSums(countMatrix)
	countMatrix <- countMatrix[,librarySize >= minLSize]
	librarySize <- colSums(countMatrix)
	mtGenes <- grep('^MT-',toupper(rownames(countMatrix)))
	nGenes <- colSums(countMatrix != 0)
  
	genesLM <- lm(nGenes~librarySize)
	genesLM <- as.data.frame(predict(genesLM, data.frame(librarySize), interval = 'prediction'))
  
	if(isTRUE(length(mtGenes) > 0)){
		mtCounts <- colSums(countMatrix[grep('^MT-',toupper(rownames(countMatrix))),])
		mtProportion <- mtCounts/librarySize
		mtLM <- lm(mtCounts~librarySize)
		mtLM <- as.data.frame(predict(mtLM, data.frame(librarySize), interval = 'prediction'))
		selectedCells <- ((mtCounts > mtLM$lwr) & (mtCounts < mtLM$upr) & (nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (mtProportion <= mt.threshold) & (librarySize < 2 * mean(librarySize)))
	} else {
		selectedCells <- ((nGenes > genesLM$lwr) & (nGenes < genesLM$upr) & (librarySize < 2 * mean(librarySize)))
	}
  	selectedCells <- colnames(countMatrix)[selectedCells]
	if(class(mySO) == 'Seurat'){
		mySO <- subset(mySO, cells = selectedCells)
	} else {
		mySO <- countMatrix[,selectedCells]
	}
	return(mySO)

}

#source required libraries
print("Loading Required Libraries and Functions")
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
	library(pkg)
}

#loading OrgDb
library(OrgDb,character.only=TRUE)

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
qcSample <- function(outputFolder, mySO, sample_name, mtPattern, riboPattern){

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
	plotFile <- paste0(outputFolder,"/plots/qcPlots/",sample_name,"_raw")
	plotQC <- qcPlot(mySO@assays$RNA@counts)
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

        plotFile <- paste0(outputFolder,"/plots/qcPlots/",sample_name,"_filtered")
        plotQC <- qcPlot(mySO@assays$RNA@counts)
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
		geom_violin() + theme_light +
		geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
		xlab('') + ylab('Library Size') +
		labs(subtitle = parse(text = paste0('italic(n) ==', ncol(mySO)))

        QC2 <- ggplot(qc.df, aes('', num_features)) +
                geom_violin() + theme_light +
                geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
                xlab('') + ylab('nFeatures') 

        QC3 <- ggplot(qc.df, aes('', mito_proportion)) +
                geom_violin() + theme_light +
                geom_jitter(width = 0.2, color = rgb(0,0,0,0.05)) +
                xlab('') + ylab('mtProportion')
	
	plot_obj <- QC1 + QC2 + QC3
	return(plot_obj) 
}









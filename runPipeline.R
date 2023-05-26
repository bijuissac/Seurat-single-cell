#This code is to run the single cell RNA seq pipeline 
print("Starting Analysis")
print("Sourcing Run Parameters")
source("./runParameters.R")
session_info()

print("Initializing Output Folder")
try(dir.create(outputFolder), silent = TRUE)

print("Moving to Output Folder")
setwd(outputFolder)

print("Creating output directory structure")
for(subFolder in outputSubFolder){
	try(dir.create(subFolder), silent = TRUE)
}

if(run_conditions=="all" || run_conditions=="merge"){

 print("Read Sample Data for QC")
 index = 1
 for(myCondition in expCondTable$expConditions){
	myReplicates <- expCondTable$expReplicates[[index]]
	index_r = 1
	additional_samples <- c()
	for(myReplicate in myReplicates){
		mySO <- readInData(runSoupX,data_dir,myReplicate,inputFormat)

		#QC sample
		mySO.filtered <- qcSample(qc_filtering_method, outputFolder, mySO, myReplicate, mt.threshold, mtPattern, riboPattern)	

		if(index_r == 1){
			sample1 <- mySO.filtered
		}else
		{
			additional_samples <- append(additional_samples, mySO.filtered)
		}
		index_r = index_r + 1
	}

	#merge replicates
	if(index_r > 1){
		mySO.merged <- merge(sample1, y = additional_samples, add.cell.ids = myReplicates, project = outputFolder)
	}else
	{
		mySO.merged <- sample1
	}

        #save filtered RDS file
        rds_file <- paste0(outputFolder,"/savedData/",myCondition,".rds")
        saveRDS(mySO.filtered, file = rds_file)
	index = index + 1
 }

}


if(run_conditions=="all" || run_conditions=="integrate"){
	
	mySO.list <- list()
	index = 1
	mySO.combined <- c()
	for(myCondition in expCondTable$expConditions){
		print(paste0("Reading Saved RDS file : ",myCondition))
		rds_file <- paste0(outputFolder,"/savedData/",myCondition,".rds")
		mySO <- readRDS(rds_file)
		mySO <- AddMetaData(mySO, metadata = rep(myCondition, ncol(mySO)), col.name = "Condition")
		mySO.list[[myCondition]] <- mySO
		ifelse(index==1,mySO.combined <- mySO, mySO.combined <- merge(mySO.combined, y = mySO, project = outputFolder))
		index = index + 1
	}
	
	#Normalization, Variable Feature Selection and Scaling
	mySO.combined <- NormalizeData(object = mySO.combined, normalization.method = "LogNormalize", scale.factor = 10000)
	mySO.combined <- FindVariableFeatures(object = mySO.combined, selection.method = 'vst', nfeatures = 2000)
	mySO.combined <- ScaleData(mySO.combined, verbose = T)

	#Spectral Clustering
	seed_val <- 1
	mySO.combined <- FindNeighbors(mySO.combined, reduction = "pca", dims = 1:20)
	mySO.combined <- FindClusters(mySO.combined, resolution = "0.8", random.seed = seed_val)

	set.seed(seed_val)
	mySO.combined <- RunPCA(mySO.combined, verbose = T, features = VariableFeatures(object = mySO.combined), npcs = npcs)
	set.seed(seed_val)
	mySO.combined <- RunTSNE(mySO.combined, reduction = "pca", dims = 1:20, seed.use = seed_val)
	set.seed(seed_val)
	mySO.combined <- RunUMAP(mySO.combined, reduction = "pca", dims = 1:20, seed.use = seed_val)

	#Changing the Ident information for Plotting
	Idents(mySO.combined) <- mySO.combined$orig.ident
	
	#Pre-Integration TSNE and UMAP plots
	plot_tsne_umap(outputFolder, mySO.combined, "beforeIntegration", "tsne")
	plot_tsne_umap(outputFolder, mySO.combined, "beforeIntegration", "umap")

	#Integrate Data
	mySO.combined <- integrate_data(mySO.combined, mySO.list, integration_method, ncps)
	
	#Post-Integration TSNE and UMAP plots
	plot_tsne_umap(outputFolder, mySO.combined, "afterIntegration", "tsne")	
	plot_tsne_umap(outputFolder, mySO.combined, "afterIntegration", "umap")

	#save data
	rds_file <- paste0(outputFolder,"/savedData/Integrated_Conditions.rds")
	saveRDS(mySO.combined, file = rds_file)

}











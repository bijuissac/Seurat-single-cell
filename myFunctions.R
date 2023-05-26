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
	return mySO
}



#This code is to run the single cell RNA seq pipeline 
print("Starting Analysis")
print("Sourcing Run Parameters")
source("./runParameters.R")

session_info()

print("Initializing Output Folder")
try(dir.create(outputFolder), silent = TRUE)

print("Creating output directory structure")
for(subFolder in outputSubFolder){
	try(dir.create(subFolder), silent = TRUE)
}

print("Read Sample Data for QC")
index = 1
for(myCondition in expCondTable$expConditions){
	myReplicates <- expCondTable$expReplicates[[index]]
	index_r = 1
	for(myReplicate in myReplicates){
		mySO <- readInData(runSoupX,data_dir,myReplicate,inputFormat)

		#QC sample
		mySO.filtered <- qcSample(outputFolder, mySO, myReplicate, mtPattern, riboPattern)	

		#save filtered RDS file
		rds_file <- paste0(outputFolder,"/rds/",myCondition,"_",myReplicate,".rds")
		saveRDS(mySO.filtered, file = rds_file)

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
        rds_file <- paste0(outputFolder,"/rds/",myCondition,".rds")
        saveRDS(mySO.filtered, file = rds_file)

	index = index + 1
}



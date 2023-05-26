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
	for(myReplicate in myReplicates){
		mySO <- readInData(runSoupX,data_dir,myReplicate,inputFormat)
		rds_file <- paste0(outputFolder,"/rds/",myCondition,"_",myReplicate,".rds")
		saveRDS(mySO,file = rds_file)
	}
	index = index + 1
}



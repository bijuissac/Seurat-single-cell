#This file defines the running parameters for single cell RNAseq pipeline
#Input Directories
base_dir <- "/home/user/scrnaseq/"
data_dir <- paste0(base_dir,"data/")
code_dir <- paste0(base_dir,"code/")

#R library locations and functions
r_library <- paste0(base_dir,"r4.2.3/library/")
source(paste0(base_dir,"code/myFunctions.R")

#Estimate and remove cell-free mRNA contamination
runSoupX = FALSE #TRUE/FALSE

#Input Data Format
inputFormat = "MEX" #MEX/H5/H5Seurat

#Project Name
projectName = "myProject"

#Output Folder
outputFolder = paste0(base_dir,projectName)

#Conditions and Replicates
expConditions <- list("OVA","SAL")
expReplicates <- list(c("EC-OVA-1","EC-OVA-2","EC-OVA-3"),c("EC-SAL-1","EC-SAL-2","EC-SAL-3"))

#species
myOrganism <- "mouse" #mouse/human

#mitochondrial and ribosomal fraction thresholds
mt.threshold <- "0.05" #mouse = 0.05/human = 0.10
ribo_fraction = 1

#filtering method
#myQC will dynamically calculate feature thresholds while Seurat is a hard threshold
qc_filtering_method = "Seurat" #Seurat/myQC/percentile
min_features = 200
max_features = 2500
if(qc_filtering_method=="percentile"){
	mt.threshold = 1
}

#Integration method
integration_method <- "harmony" #seurat/harmony

#Default dimensions and resolutions and number of principal components
dimensions <- c(seq(10,100,10))
resolutions <- c(0.8) 
npcs <- 100

#Default Pathway parameters
show_top_pathways <- 10

#++++++++++++++++++++++++++++++++++++
#Do not edit beyond this line
#++++++++++++++++++++++++++++++++++++
#Output Folder Structure
outputSubFolder <- c("savedData", "differentialExpression", "pathwayEnrichment", "plots", "logFiles", "topMarkers",
			"pathwayEnrichment/KEGG", "pathwayEnrichment/GO_BP", "pathwayEnrichment/GO_CC", "pathwayEnrichment/GO_MF",
			"plots/featurePlots", "plots/clustering", "plots/pathways", "plots/pathways/KEGG", "plots/pathways/GO_BP",
			"plots/pathways/GO_CC", "plots/pathways/GO_MF")

#Conditions and replicates table
expCondTable <- data.table(expConditions, expReplicates)

#Organism specific database
organisms <- list()
organisms[["mouse"]] <- c("org.Mm.eg.db","mmu","^mt-","^Rpl|^Rps",mt.threshold)
organisms[["human"]] <- c("org.Hs.eg.db","hsa","^MT-","^RPL|^RPL",mt.threshold)

OrgDb <- toString(organisms[[myOrganism]][1])
organism <- toString(organisms[[myOrganism]][2])
mtPattern <- toString(organisms[[myOrganism]][3])
riboPattern <- toString(organisms[[myOrganism]][4])

#Integration method parameters
integration_methods <- list()
integration_methods[["seurat"]] <- c("integrated","pca")
integration_methods[["harmony"]] <- c("harmony","harmony")
reduction_method <- toString(integration_methods[[integration_method]][2])

#SingleR annotation datasets
singleRDatasets <- list()
singleRDatasets[["human"]] <- ImmGenData() #human immune dataset
singleRDatasets[["mouse"]] <- MouseRNAseqData()
reference_dataset <- singleRDatasets[[myOrganism]]
reference_cell_types <- reference_dataset$label.main






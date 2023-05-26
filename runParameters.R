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
expConditions <- list("Control","Treatment")
expReplicates <- list(c("C1","C2"),c("T1","T2"))

#species
myOrganism <- "mouse" #mouse/human

#mitochondrial thresholds
mt.threshold <- "0.05" #mouse = 0.05/human = 0.10

#filtering method
#
qc_filtering_method = "Seurat" #Seurat/myQC



#++++++++++++++++++++++++++++++++++++
#Do not edit beyond this line
#++++++++++++++++++++++++++++++++++++
#Output Folder Structure
outputSubFolder <- c("savedData", "differentialExpression", "pathwayEnrichment", "plots", "logFiles", "topMarkers",
			"pathwayEnrichment/KEGG", "pathwayEnrichment/GO_BP", "pathwayEnrichment/GO_CC", "pathwayEnrichment/GO_MF",
			"plots/featurePlots", "plots/clustering", "plots/pathways", "plots/pathways/KEGG", "plots/pathways/GO_BP",
			"plots/pathways/GO_CC", "plots/pathways/GO_MF")

#Organism specific database
organisms <- list()
organisms[["mouse"]] <- c("org.Mm.eg.db","mmu","^mt-",mt.threshold)
organisms[["human"]] <- c("org.Hs.eg.db","hsa","^MT-",mt.threshold)

OrgDb <- toString(organisms[[myOrganism]][1])
organism <- toString(organisms[[myOrganism]][2])











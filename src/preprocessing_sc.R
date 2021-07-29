# [1]:sc_data_dir  [2]:project_name ([3]:metadata(record the source clone))

args <- commandArgs(TRUE)
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
# for convert ensembl ids to gene symbols
#library(singleCellTK)
suppressMessages(library(org.Hs.eg.db))

#path <- args[1]
#print(sprintf("output folder: %s", dirname(path)))
#setwd(dirname(path))

setwd('/input_dir')

#----------------------------------------#
# Load in the single cell dataset        #
#----------------------------------------#
# if input is a table (.csv, .txt)
sc.data <- read.table(file = args[1], sep=",")
# if input is 10X folder
#sc.data <- Read10X(data.dir = basename(path))

#----------------------------------------------------#
# Convert ensembl IDs to gene symbols (if needed)    #
#----------------------------------------------------#
#sc.data <- convertGeneIDs(inSCE = sc.data, inSymbol = "ENSEMBL", outSymbol = "SYMBOL")
print(dim(sc.data))

#----------------------------------------#
# Build Seurat Objects                   #
#----------------------------------------#
# keep the genes that express in more than ${min.cells} cells
min.cells <- as.integer(dim(sc.data)[1]*0.001)

if(length(args) > 2){
	# read in metadata (records the clone source of each cell) (only MCF-7)
	metadata <- read.table(file = args[3], sep="\t")
	print(head(metadata))
	sc <- CreateSeuratObject(counts = sc.data, project = args[2], min.cells=min.cells, meta.data=metadata)
}else{
	sc <- CreateSeuratObject(counts = sc.data, project = args[2], min.cells=min.cells)
}
sc

# visualization for QC
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
#pdf(sprintf("%s_QC.pdf", args[2]))
#VlnPlot(sc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

#----------------------------------------#
# Filtration                             #
#----------------------------------------#
# for GSE127471 (PBMCs)
#sc <- subset(sc, subset = nFeature_RNA >= 500 & nFeature_RNA <= 3500 & percent.mt < 10 & percent.mt > 1)
# for GSE114459 (MCF-7)
sc <- subset(sc, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000 & percent.mt > 1 & percent.mt <= 15)
print("after filtration: ")
print(dim(sc))

#----------------------------------------#
# Normalization                          #
#----------------------------------------#
sc <- NormalizeData(sc, verbose=FALSE, scale.factor=1e5)

#----------------------------------------#
# Assign each cell a phase               #
#----------------------------------------#
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#----------------------------------------#
# Save processed data                    #
#----------------------------------------#
setwd('/output_dir')
saveRDS(sc, file="seurat_onlyfiltration.RDS")
head(sc[[]])
q()

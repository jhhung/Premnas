args<-commandArgs(TRUE)
suppressMessages(require(ACTIONet))
suppressMessages(require(SingleCellExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))



setwd(args[1])

sc <- readRDS(args[2])
DefaultAssay(sc) <- "RNA"


# Select cells that pass the confidence threshold
#cluster_df <- read.table(args[2], header=T,col.names=c("cluster_name","cluster"))
cluster_df <- read.table(args[3], header=F, sep=" ",col.names=c("cluster_name","cluster"))
clusters <- as.vector(cluster_df %>% select(cluster))

sc@meta.data$cell.type <- clusters


sc_highconf <- subset(x = sc, subset = cell.type!=0)
#sc_highconf <- sc
Idents(object = sc_highconf) <- "cell.type"

# remove cluster containing less than 3 cells
cluster.table <- table(sc_highconf$cell.type)
print(cluster.table)
remove.list <- list()
for (it in c(1:length(names(cluster.table)))){
    if(as.numeric(cluster.table[it])<3){
        remove.list <- c(remove.list, it)
    }
}


if(length(remove.list)>0){
    print("Removing the clusters containing less than 3 cells")
    sc_highconf <- subset(x = sc_highconf, subset = cell.type != remove.list)
    print(table(sc_highconf$cell.type))
}
type <- sc_highconf$cell.type
DefaultAssay(sc_highconf) <- "RNA"

# Identify DEGs
print("Identify DEGs")
markers <- FindAllMarkers(sc_highconf, assay="RNA", return.thresh = 0.0001, only.pos = TRUE)
sig_genes = unique(markers %>% filter(p_val_adj< 1e-2) %>% filter(avg_log2FC > 0) %>% pull(gene))
Fun = function(DF) {
        DF = DF %>% filter(p_val_adj< 1e-4) %>% arrange(desc(avg_log2FC))
        genes = DF %>% top_n(20) %>% pull(gene)
#        pdf(paste0("filtered_significant_genes_",unique(DF$cluster),".pdf"))
#        VlnPlot(sc_highconf, features=genes, pt.size=0.5,ncol=2)
#        dev.off()
        return(DF)
}
markers %>% group_by(cluster) %>% do(Fun(.))


# Nomalize data
sc_highconf <- NormalizeData(sc_highconf,normalization.method = "RC", scale.factor=1e6)

# Build scRNA-seq reference profile (CPM)
sc_ref <- data.frame(as.matrix(GetAssayData(sc_highconf, assay="RNA", slot="data")))
print(head(colSums(sc_ref)))
colnames(sc_ref) <- paste0("Cluster ",type)
sc_ref <- sc_ref[ , order(names(sc_ref))]
colnames(sc_ref) <- gsub(pattern="\\.[0-9]+", replacement="", colnames(sc_ref))
write.table(sc_ref, file=sprintf("Subpopulation-characteristic.txt"), sep="\t", quote=FALSE)

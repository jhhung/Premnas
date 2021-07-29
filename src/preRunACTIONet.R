args <- commandArgs(TRUE)
suppressMessages(require(ACTIONet))
suppressMessages(require(Seurat))
suppressMessages(require(SingleCellExperiment))
setwd(args[1])

sc <- readRDS(args[2])


#sce <- import.sce.from.Seurat(sc)
sce <- as.SingleCellExperiment(sc)
sce

class(sce)

batch_attr <- interaction(as.vector(sc$source), as.vector(sc$Phase))

ace.reduced <- reduce.and.batch.correct.ace.Harmony(sce, batch_attr)
saveRDS(ace.reduced, file=sprintf("sc-after-reduce-and-batchcorrect.RDS"))







args<-commandArgs(TRUE)
suppressMessages(require(ACTIONet))
suppressMessages(require(SingleCellExperiment))
suppressMessages(require(dplyr))

setwd(args[1])

ace = readRDS(args[2])
ace = run.ACTIONet(ace,layout_compactness=0,thread_no=30,layout_epochs = 2500)

saveRDS(ace, file = "ACTIONet-model")

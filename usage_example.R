library("rscala")
library("SingleCellExperiment")

source("./src/main/R/micro_cell_clust.R")
data = read.csv("./data/toy_dataset_not_normalized.csv", header = FALSE)

mcc.rs = initMccRscala()

results = runAutomaticTuning(mcc.rs, data)

print(results$cells.idx)

print(results$genes.idx)

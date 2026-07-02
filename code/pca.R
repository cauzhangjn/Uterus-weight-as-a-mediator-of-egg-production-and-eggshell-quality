#### read args
args <- commandArgs(TRUE)
trait <- args[1]

setwd(getwd())
eigvec <- read.table(paste0(trait, ".eigenvec"), header = F, stringsAsFactors = F)
write.table(eigvec[2:ncol(eigvec)], file = paste0(trait, ".eigenvector.xls"), sep = "\t", row.names = F, col.names = T, quote = F)

eigval <- read.table(paste0(trait, ".eigenval"), header = F)
pcs <- paste0("PC", 1:nrow(eigval))
eigval[nrow(eigval), 1] <- 0
percentage <- eigval$V1 / sum(eigval$V1) * 100
eigval_df <- as.data.frame(cbind(pcs, eigval[, 1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)
write.table(eigval_df, file = paste0(trait, ".eigenvalue.xls"), sep = "\t", quote = F, row.names = F, col.names = T)

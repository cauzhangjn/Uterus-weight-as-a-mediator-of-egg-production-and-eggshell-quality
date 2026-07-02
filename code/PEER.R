# Author: Francois Aguet
library(peer, quietly=TRUE)  # https://github.com/PMBio/peer
library(argparser, quietly=TRUE)
con <- file("BJmerged_matrix.log") 
sink(con, append=TRUE)##sink(con, append=TRUE) 的作用是将R中的输出写入到一个指定的文件中
sink(con, append=TRUE, type="message")

WriteTable <- function(data, filename, index.name) {
  datafile <- file(filename, open = "wt")
  on.exit(close(datafile))
  header <- c(index.name, colnames(data))
  writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
  write.table(data, datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

##这个函数用于帮助在 R 脚本中处理命令行输入，使得用户能够在运行脚本时指定一些参数。
p <- arg_parser("Run PEER factor estimation")
##用于向命令行参数解析器对象中添加具体的命令行参数规则。
p <- add_argument(p, "expr.file", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "n", help="Number of hidden confounders to estimate")
p <- add_argument(p, "--covariates", help="Observed covariates")
p <- add_argument(p, "--alphaprior_a", help="", default=0.001)
p <- add_argument(p, "--alphaprior_b", help="", default=0.01)
p <- add_argument(p, "--epsprior_a", help="", default=0.1)
p <- add_argument(p, "--epsprior_b", help="", default=10)
p <- add_argument(p, "--max_iter", help="", default=1000)
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default=".")
##parse_args 函数，传递参数解析器对象 p 作为参数。
argv <- parse_args(p)

##提示用户或程序运行者当前脚本正在执行的操作

cat("PEER: loading expression data ... ")
##grepl用于在字符向量中搜索模式的函数，它返回一个逻辑向量
##如果是.gz格式的结尾，那就zcat解压缩，再统计行，如果不是那就直接统计。
if (grepl('.gz$', argv$expr.file)) {
  ##system用于在系统上执行shell命令
  nrows <- as.integer(system(paste0("zcat ", argv$expr.file, " | wc -l | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
} else {
  nrows <- as.integer(system(paste0("wc -l ", argv$expr.file, " | cut -d' ' -f1 "), intern=TRUE, wait=TRUE))
}
##||：逻辑 OR 运算符
if (grepl('.bed$', argv$expr.file) || grepl('.bed.gz$', argv$expr.file)) {
  df <- read.table(argv$expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="")
  row.names(df) <- df[, 4]
  df <- df[, 5:ncol(df)]
} else {
  df <- read.table(argv$expr.file, sep="\t", nrows=nrows, header=TRUE, check.names=FALSE, comment.char="", row.names=1)
}
M <- t(as.matrix(df))
cat("done.\n")

# run PEER
cat(paste0("PEER: estimating hidden confounders (", argv$n, ")\n"))
model <- PEER()
##PEER_setNk是PEER 包中的一个函数，用于设置模型中隐藏混杂因子的数量
invisible(PEER_setNk(model, argv$n))
##M 是表达数据矩阵，表示样本和基因之间的表达值
invisible(PEER_setPhenoMean(model, M))
##设置模型中先验参数 Alpha，a，b通过命令行参数解析器获取的先验参数的值。
##这个先验参数通常影响模型中的正则化，对于控制模型的复杂度和防止过拟合很重要
invisible(PEER_setPriorAlpha(model, argv$alphaprior_a, argv$alphaprior_b))
##用于设置模型中先验参数 Eps，与数据的噪声水平相关，对于调整模型对观测数据的拟合有影响。
invisible(PEER_setPriorEps(model, argv$epsprior_a, argv$epsprior_b))
##是通过命令行参数解析器获取的最大迭代次数的值。
invisible(PEER_setNmax_iterations(model, argv$max_iter))
##是否额外提供协变量文件
if (!is.null(argv$covariates) && !is.na(argv$covariates)) {
  has.cov <- TRUE
  covar.df <- read.table(argv$covariates, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
  covar.df[] <- sapply(covar.df, as.numeric)
  cat(paste0("  * including ", dim(covar.df)[2], " covariates", "\n"))
  invisible(PEER_setCovariates(model, as.matrix(covar.df[rownames(M), ])))  # samples x covariates
} else {
  has.cov <- FALSE
}
time <- system.time(PEER_update(model))
##使用 PEER_getX 函数从 PEER 模型中获取估计的 PEER 因子矩阵 X。
##X 是一个矩阵，其行代表样本（samples），列代表 PEER 因子。
##这个矩阵表示每个样本在估计的 PEER 因子上的权重。
X <- PEER_getX(model)  # samples x PEER factors

##使用 PEER_getAlpha 函数从 PEER 模型中获取估计的 PEER 因子的权重矩阵 A。
##A 是一个矩阵，其列向量表示每个 PEER 因子的权重。
##这个矩阵反映了每个 PEER 因子在解释表达数据方差中的相对重要性。
A <- PEER_getAlpha(model)  # PEER factors x 1

##使用 PEER_getResiduals 函数从 PEER 模型中获取估计的剩余矩阵 R。
##R 是一个矩阵，其行代表基因，列代表样本。
##这个矩阵包含了模型无法解释的残差，即去除了 PEER 因子的表达数据残差。
R <- t(PEER_getResiduals(model))  # genes x samples

# add relevant row/column names
if (has.cov) {
  cols <- c(colnames(covar.df), paste0("InferredCov",1:(ncol(X)-dim(covar.df)[2])))
} else {
  cols <- paste0("InferredCov",1:ncol(X))
}
rownames(X) <- rownames(M)
colnames(X) <- cols
rownames(A) <- cols
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha
rownames(R) <- colnames(M)
colnames(R) <- rownames(M)

# write results
cat("PEER: writing results ... ")
WriteTable(t(X), file.path(argv$output_dir, paste0(argv$prefix, ".PEER_covariates.txt")), "ID")  # format(X, digits=6)
WriteTable(A, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_alpha.txt")), "ID")
WriteTable(R, file.path(argv$output_dir, paste0(argv$prefix, ".PEER_residuals.txt")), "ID")
cat("done.\n")



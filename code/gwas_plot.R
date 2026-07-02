library(CMplot)

args <- commandArgs(TRUE)
trait <- args[1]

setwd(getwd())
data<-read.table("duli_snp_number.txt", header = FALSE, nrows = 1)
k<-data[1,1]
print(k)

df <- read.table(paste0(trait, '.fastGWA'), header=T)
df <- na.omit(df[, c(2,1,3,10)])
names(df) <- c("SNP","Chromosome","Position","P-value")
df <- df[df[, 4] >0 & df[, 4] <1, ]

p_value=df[,4]
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
print(paste0('*** lambda: ', lambda))

CMplot(df, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(0.05/k, 1/k),threshold.lty=c(1,2),
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
        chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),signal.cex=c(1,1),
        signal.pch=c(19,19),file="png",file.name=trait,dpi=300,file.output=TRUE,verbose=TRUE)

CMplot(df,plot.type="q",conf.int.col=NULL,box=TRUE,file="png",dpi=300,file.name=trait,
        file.output=TRUE,verbose=TRUE,signal.pch=19,signal.cex=1.5,signal.col="red")

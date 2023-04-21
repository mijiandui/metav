library(getopt)

# get options
spec <- matrix(c(
  'help', 'h', 0, 'logical', 'help',
  'input', 'i', 1, 'character', 'input dir contains predict.list and likelihood.matrix'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# check options
if (!is.null(opt$help) | is.null(opt$input)){
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}

# main
likelihood_matrix_file <- sprintf("%s/llikelihood.matrix", opt$input)
prediction_list_file <- sprintf("%s/prediction.list", opt$input)

# check existence of file
if (!file.exists(likelihood_matrix_file)){
  cat("likehood matrix file dose not exist")
  q(status = 1)
}

if (!file.exists(prediction_list_file)){
  cat("prediction list file dose not exist")
  q(status = 1)
}

# load data
library(data.table)
library(parallel)

llikelihood = fread(likelihood_matrix_file, header=TRUE)
prediction <- fread(prediction_list_file, header=TRUE)


llikelihood = read.table("llikelihood.matrix")
means = apply(X=llikelihood,FUN=mean,MARGIN=1)
sigmas = apply(X=llikelihood,FUN=sd,MARGIN=1)

nullParam = cbind(means,sigmas)

write.table(nullParam,file="nullParameters.tsv",quote=F,sep = "\t",col.names=F)
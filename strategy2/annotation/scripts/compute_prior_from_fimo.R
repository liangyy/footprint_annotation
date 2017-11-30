library(optparse)

option_list <- list(
  make_option(c("-in", "--input"), type="character", default=NULL,
              help="input file name",
              metavar="character"),
  make_option(c("-ms", "--model_str"), type="character", default=NULL,
              help="model string from where you can find model rds file by replacing wildcards with training data name and motif name",
              metavar="character"),
  make_option(c("-dn", "--data_name"), type="character", default=NULL,
              help="the training data name in model name",
              metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(stringr)
library(dplyr)

snvs <- tryCatch(
  {
    read.table(opt$input, sep = '\t', header = FALSE)
  },
  error = function(cond) {
    return(NA)
  }
)

# snvs <- read.table(opt$input, sep = '\t', header = FALSE)
if(is.na(snvs)) {
  gz <- gzfile(opt$out, 'w')
  write.table(data.frame(), gz,
            quote = F, col.names = F, row.names = F, sep = '\t')
  close(gz)
  quit()
}
snvs$id <- 1 : nrow(snvs)
snvs$prior1 <- NaN
snvs$prior2 <- NaN
model_path <- str_replace(opt$model_str, '\\{data_name\\}', opt$data_name)
for(motif in unique(snvs$V7)) {
  model <- str_replace(model_path, '\\{motif\\}', motif)
  model.param <- readRDS(model)
  snvs.subset <- snvs %>%
    filter(V7 == motif) %>%
    select(V8, V9, id)
  prior1.x <- data.frame(Ict = 1, Pwm = snvs.subset$V8)
  prior1.y <- as.matrix(prior1.x) %*% model.param$BetaLogit
  prior2.x <- data.frame(Ict = 1, Pwm = snvs.subset$V9)
  prior2.y <- as.matrix(prior2.x) %*% model.param$BetaLogit
  snvs[snvs.subset$id, 'prior1'] <- prior1.y
  snvs[snvs.subset$id, 'prior2'] <- prior2.y
}
gz <- gzfile(opt$out, 'w')
write.table(snvs[, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'prior1', 'prior2', 'V7', 'V8', 'V9')], gz,
            quote = F, col.names = F, row.names = F, sep = '\t')
close(gz)

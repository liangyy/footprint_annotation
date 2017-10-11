library(optparse)

option_list <- list(
    make_option(c("-r", "--active_region"), type="character", default=NULL,
                help="Regions of interest along with motif name and score",
                metavar="character"),
    make_option(c("-f", "--five_prime_count"), type="character", default=NULL,
                help="Five prime count from sequencing experiment",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Name of output RDS file",
                metavar="character"),
                make_option(c("-s", "--signal"), type="character", default=NULL,
                            help="Name of output signal file",
                            metavar="character"),
    make_option(c("-e", "--extend_win"), type="numeric", default=NULL,
                help="Extended size for window (half of window size setup)",
                metavar="character"),
    make_option(c('-c', '--centipede_path'), type='character', default=NULL,
                help='The path of centipede script you want to use',
                metavar="character"),
    make_option(c('-p', '--plot_prefix'), type='character', default=NULL,
                help='Prefix of the output plot',
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
library(tidyr)
convertSignalToScore <- function(start, end, score) {
  out <- c()
  dup.ind <- duplicated(start)
  start <- start[!dup.ind]
  end <- end[!dup.ind]
  score <- score[!dup.ind]
  for(i in 1 : length(start)) {
    n <- end[i] - start[i]
    out <- c(out, rep(score[i], n))
  }
  out <- paste(out, collapse = ',')
  return(out)
}

if (opt$centipede_path == 'NO') {
  library(CENTIPEDE)
} else {
    source(opt$centipede_path)
  fitCentipede <- fitCentipede3
}

readCutSite <- function(filename) {
  cutsite.readin <- read.table(filename, sep = '\t', header = F)
  cutsite <- cutsite.readin %>%
    mutate(id = paste(V1, V2, V3, V4, V5))
  cutsite$V6 <- as.character(cutsite$V6)
  cutsite.count <- strsplit(cutsite$V6, ',')
  n <- length(cutsite.count[[1]])
  l <- lapply(cutsite.count, length)
  boundary.ind <- l != n
  cutsite.complete <- cutsite[!boundary.ind, ]
  cutsite.complete <- cutsite.complete %>%
    mutate(count.str = V6) %>%
    select(id, count.str)
  return(cutsite.complete)
}

doUnlist <- function(input) {
  n <- length(input[[1]])
  input <- unlist(input)
  class(input) <- 'numeric'
  input <- t(matrix(input, nrow = n))
  return(input)
}

pwm.readin <- read.table(opt$active_region, sep = '\t', header = F)
pwm <- pwm.readin %>%
  mutate(id = paste(V1, V2, V3, V4, V7)) %>%
  group_by(id) %>%
  summarise(pwm.score = V5[1]) %>%
  ungroup()
cutsite.forward <- readCutSite(opt$five_prime_count_f)
cutsite.backward <- readCutSite(opt$five_prime_count_b)
cutsite.backward <- cutsite.backward %>%
  mutate(count.str.back = count.str) %>%
  select(id, count.str.back)
cutsite.complete <- cutsite.forward %>%
  inner_join(pwm, by = 'id') %>%
  inner_join(cutsite.backward, by = 'id')

cutsite.complete.count.f <- strsplit(cutsite.complete$count.str, ',')
cutsite.complete.count.f <- doUnlist(cutsite.complete.count.f)
cutsite.complete.count.b <- strsplit(cutsite.complete$count.str.backward, ',')
cutsite.complete.count.b <- doUnlist(cutsite.complete.count.b)
cutsite.complete.count <- cbind(cutsite.complete.count.f, cutsite.complete.count.b)

model <- fitCentipede(
  Xlist = list(Seq = cutsite.complete.count),
  Y = as.matrix(data.frame(Ict = 1,Pwm = cutsite.complete$pwm.score)),
  DampLambda = 0.1,
  DampNegBin = 0.001,
  sweeps = 200)

model$DataLogRatio <- new.fit$NegBinLogRatio + new.fit$MultiNomLogRatio
ct <- cor.test(jitter(cutsite.complete$pwm.score), jitter(model$DataLogRatio), method = 'spearman')
cat("#Spearman_p.val = ", ct$p.value, "_rho = ", ct$estimate, "\n")
out <- list(model = model, data = cutsite.complete)
saveRDS(out, file = opt$output)

signal <- out$data[out$model$PostPr > 0.9, ]
signal.info <- unlist(strsplit(signal$id, ' '))
signal.info <- t(matrix(signal.info, nrow = 5))
signal.info <- data.frame(signal.info)
signal.info$X2 <- as.numeric(as.character(signal.info$X2)) + opt$extend_win
signal.info$X3 <- as.numeric(as.character(signal.info$X3)) - opt$extend_win
signal.info$score <- signal$pwm.score
gz1 <- gzfile(opt$signal, "w")
write.table(signal.info, gz1, sep = '\t', col.names = F, row.names = F, quote = F)
close(gz1)

png(paste0(opt$plot_prefix, '_cutsite.png'))
imageCutSites(cutsite.complete.count[ order(model$PostPr),][c(1 : 100, (dim(cutsite.complete.count)[1] - 100) : (dim(cutsite.complete.count)[1])),])
dev.off()

png(paste0(opt$plot_prefix, '_footprint.png'))
plotProfile(model$LambdaParList[[1]], Mlen = dim(cutsite.complete.count)[1] / 2 - opt$extend_win)
dev.off()

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
    make_option(c('-c', '--centipede_path'), type='character', default=NULL,
                help='The path of centipede script you want to use',
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
library(tidyr)
convertSignalToScore <- function(start, end, score) {
  out <- c()
  for(i in 1 : length(start)) {
    n <- end[i] - start[i]
    out <- c(out, rep(score[i], n))
  }
  out <- paste(out, collapse = ',')
  return(out)
}

if (opt$centipede_path == '') {
  library(CENTIPEDE)
} else {
  fitCentipede <- fitCentipede3
}

cutsite.readin <- read.table(opt$five_prime_count, sep = '\t', header = F)
cutsite <- cutsite.readin %>%
  mutate(id = paste(V5, V6, V7, V8, V9)) %>%
  group_by(id) %>%
  summarise(count = convertSignalToScore(V2, V3, V4)) %>%
  ungroup()
# cutsite <- strsplit(cutsite.readin, ';')
# cutsite <- unlist(cutsite)
# n <- length(cutsite[[1]])
# class(cutsite) <- 'numeric'
# cutsite <- t(matrix(cutsite, nrow = n))
pwm.readin <- read.table(opt$active_region, sep = '\t', header = F)
pwm <- pwm.readin %>%
  mutate(id = paste(V1, V2, V3, V4, V5)) %>%
  group_by(id) %>%
  summarise(pwm.score = V6[1]) %>%
  ungroup()
cusite <- cutsite %>%
  inner_join(pwm, )
pwm <- pwm.readin$V6


model <- fitCentipede(Xlist=list(Seq=cutsite),Y=as.matrix(data.frame(Ict=1,Pwm=pwm)),
												 DampLambda = 0.1, DampNegBin = 0.001,sweeps=200);

ct <- cor.test(jitter(pwm), jitter(model$DataLogRatio), method = 'spearman')
cat("#Spearman_p.val = ", ct$p.value, "_rho = ", ct$estimate, "\n")

saveRDS(model, file = opt$output)

library(optparse)

option_list <- list(
    make_option(c("-f", "--first"), type="character", default=NULL,
                help="First bed file. It should be paper output",
                metavar="character"),
    make_option(c("-s", "--second"), type="character", default=NULL,
                help="Second bed file. It should be strategy2 output",
                metavar="character"),
    make_option(c("-o", "--out"), type="character", default=NULL,
                help="Output file name",
                metavar="character"),
    make_option(c("-n", "--name"), type="character", default=NULL,
                help="Name shown in the row",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)

first <- read.table(opt$first, sep = '\t', header = F)
second <- read.table(opt$second, sep = '\t', header = T)

first <- first %>%
  mutate(id = paste(V1, V2, V6, V4))
second <- second %>%
  mutate(id = paste(X1, X2, X4, X5))
i <- length(intersect(first$id, second$id))
u <- length(union(first$id, second$id))
f <- length(first$id)
s <- length(second$id)
d <- data.frame(i = i, u = u, f = f, s = s, n = opt$name)
write.table(x = d, file = opt$out, sep = '\t', col.names = F, row.names = F, quote = F)

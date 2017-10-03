library(optparse)

option_list <- list(
    make_option(c("-r", "--active_region"), type="character", default=NULL,
                help="Regions of interest along with motif name and score",
                metavar="character"),
    make_option(c("-p", "--mpile"), type="character", default=NULL,
                help="Pile-up file generated from sequencing experiment",
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

if (opt$centipede_path == '') {
  library(CENTIPEDE)
} else {
  fitCentipede <- fitCentipede3
}

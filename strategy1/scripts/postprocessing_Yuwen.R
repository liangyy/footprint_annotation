snv <- read.table('./snpLists/170607_for_Yanyu_mut_info.txt', sep = '\t', header = T)
enhancer <- read.table('./strategy1/data/footprint_snp.enhancer.ASD_060717.bed.gz', sep = '\t', header = F)
footprint <- readRDS('./strategy1/to_Yuwen.rds')
both_in <- intersect(enhancer$V6, footprint$SNP.ID)
pos_in_footprint <- match(both_in, footprint$SNP.ID)
footprint$enhancer <- FALSE
footprint$enhancer[pos_in_footprint] <- TRUE
footprint$footprint_snp_active.ind <- footprint$footprint_snp.ind & footprint$enhancer
footprint$binding_variant_active.ind <- footprint$binding_variant.ind & footprint$enhancer
saveRDS(footprint, file = './strategy1/to_Yuwen_with_enhancer.rds')

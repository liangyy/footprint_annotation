mydata <- c()
file_str <- c('SRX026914-SRX027076-SRX027083-SRX027085-SRX027086-SRX027089-SRX027091-SRX040380-SRX040395-SRX062364-SRX121276-SRX121277-SRX121278-SRX201815')
file_ids <- strsplit(file_str, '-')[[1]]
for(i in file_ids){
    temp <- read.table(paste0('./data/output.', i, '.ASD_060717.bed'), sep = '\t', header = T)
    temp$filename <- i
    png(filename = paste0('./plot/', i, '_hist.png'))
    hist((temp$LLR.Ref - temp$LLR.Alt))
    dev.off()
    mydata <- rbind(mydata, temp)
}
ggplot(mydata) + geom_histogram(aes(x = LLR.Ref, fill = 'Ref'), alpha = .3) + geom_histogram(aes(x = LLR.Alt, fill = 'Alt'), alpha = .3)
ggplot(mydata.ag) + geom_point(aes(x = LLR.Ref, y = LLR.Alt, color = phenotype, shape = big_diff.ind)) + geom_abline(slope = 1, intercept = 0)
ggplot(mydata.ag) + geom_point(aes(x = entropy, y = abs(change), color = phenotype, shape = big_diff.ind))

ggplot(mydata.ag) + geom_point(aes(x = LLR.Ref, y = LLR.Alt))
mydata.ag[mydata.ag$big_diff.ind & mydata.ag$phenotype == 'dnv_proband',]
mydata.ag[mydata.ag$change > 5 & mydata.ag$entropy < 0.1,]

table(mydata.ag[mydata.ag$LLR.Ref > 20 & mydata.ag$change > 2,'phenotype'])

mydata.ag[mydata.ag$entropy > 1.5 & mydata.ag$change > 15 & mydata.ag$LLR.Ref > 18,]

hist(mydata$LLR.Ref)
mydata$change <- mydata$LLR.Ref - mydata$LLR.Alt
mydata$big_diff.ind <- mydata$change > 15
reg = glm(big_diff.ind~1 + Motif.ID, family=binomial, data=mydata)
summary(reg)
library(precrec)
sscurves <- evalmod(scores = reg$fitted.values, labels = mydata$big_diff.ind)
plot(sscurves)
auc(sscurves)
big_effect_motifs <- names(reg$coefficients[(abs(reg$coefficients) > 0.4)])
small_motifs <- mydata$Motif.ID[!mydata$big_diff.ind]
big_motifs <- mydata$Motif.ID[mydata$big_diff.ind]
small_motifs <- paste0('Motif.ID', small_motifs)
big_motifs <- paste0('Motif.ID', big_motifs)
sum(big_motifs %in% big_effect_motifs)
sum(small_motifs %in% big_effect_motifs)
mydata$big_effect_motif.ind <- paste0('Motif.ID', mydata$Motif.ID) %in% big_effect_motifs
ggplot(mydata) + geom_histogram(aes(x = LLR.Ref - LLR.Alt, fill = big_effect_motif.ind))

# SNP list
snp <- read.table('../snpLists/170607_for_Yanyu_mut_info.txt', sep = '\t', header = T)
mydata$phenotype <- snp[mydata$SNP.ID, 'Prediction']

mydata.ag <- c()

mydata$change2 <- mydata$LogRatioPrior.Ref - mydata$LogRatioPrior.Alt

motif_entropy <- readRDS('../motif_visualization/data/entropy.recalibratedMotifs.rds')
for(i in unique(mydata$SNP.ID)){
    temp <- mydata[mydata$SNP.ID == i,]
    temp <- temp[order(temp$change2, decreasing = T)[1],]
    temp$entropy <- motif_entropy[[as.character(temp$Motif.ID)]][temp$Relative.Pos]
    mydata.ag <- rbind(mydata.ag, temp)
    # print(i)
}
fisher.test(table(mydata.ag[, c('big_effect_motif.ind', 'phenotype')]))
ggplot(mydata.ag) + geom_density(aes(x = change, fill = phenotype), alpha = .3)

binom.test(sum(mydata.ag$phenotype == 'dnv_proband'), nrow(mydata.ag), p = sum(snp$Type == 'SNV' & snp$Prediction == 'dnv_proband') / sum(snp$Type == 'SNV'))

binom.test(sum(mydata.ag$big_diff.ind[mydata.ag$phenotype == 'dnv_proband']), sum(mydata.ag$big_diff.ind), p = sum(snp$Type == 'SNV' & snp$Prediction == 'dnv_proband') / sum(snp$Type == 'SNV'))


ggplot(mydata.ag) + geom_point(aes(x = LogRatioPrior.Ref, y = LogRatioPrior.Alt, color = phenotype))

table(mydata.ag[mydata.ag$LogRatioPrior.Ref > 2 & mydata.ag$LogRatioPrior.Alt < -2, 'phenotype'])

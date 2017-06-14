readMotif <- function(filename){
    con <- file(filename, 'r', blocking = FALSE)
    temp <- readLines(con)
    close(con)
    p.cut <- as.numeric(strsplit(temp[2], '\t')[[1]][4:12])
    p.llr <- as.numeric(strsplit(temp[4], '\t')[[1]][4:12])
    prior <- data.frame(Binding.Prior = p.cut, Motif.LLR = p.llr)
    pwm <- read.table(file = filename, comment.char = '#', header = T)
    return(list(prior=prior, pwm=pwm))
}

entropyDiscrete <- function(x){
    return(- sum(x * log2(x + 1e-10)))
}

plot_pwm <- function(plot1, pwm){
    png(plot1, width = 600, height = 400)
    seqLogo(pwm)
    dev.off()
}

plot_prior <- function(plot2, prior){
    png(plot2, width = 400, height = 400)
    ggplot(prior) + geom_point(aes(x = Motif.LLR, y = Binding.Prior))
    dev.off()
}
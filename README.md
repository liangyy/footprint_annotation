# footprint_annotation
Annotate SNV with the predicted effect of TF binding (based on the result of [this PLOS Genetics paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005875))

# Strategies

* Use the footprinting result per cell type DNase-seq experiment and compute the difference of motif score in reference and alternative allele. Recalibrate the result in some way (distribution of change under null hypothesis?; quantile threshold?) to define footprint SNP
* For each DNase-seq data set, re-train a CENTIPEDE model using the procedure described in the paper so that we can get posterior binding probability of reference/alternative allele for each SNP. Then footprint SNP can be called using the definition that "at least one of the two alleles has posterior probability greater than 0.99

# Pipelines

## Strategy 1 -- Motif Score




## Strategy 2 -- Posterior Probability



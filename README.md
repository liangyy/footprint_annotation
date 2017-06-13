# footprint_annotation
Annotate SNV with the predicted effect of TF binding (based on the result of [this PLOS Genetics paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005875))

# Strategies

* Use the footprinting result per cell type DNase-seq experiment and compute the difference of motif score in reference and alternative allele. Recalibrate the result in some way (distribution of change under null hypothesis?; quantile threshold?) to define footprint SNP
* For each DNase-seq data set, re-train a CENTIPEDE model using the procedure described in the paper so that we can get posterior binding probability of reference/alternative allele for each SNP. Then footprint SNP can be called using the definition that "at least one of the two alleles has posterior probability greater than 0.99

# Pipelines

## Strategy 1 -- Motif Score

* **Input**: motif files, SNP list, genome file
* **Procedure**:
  1. Get the list of SNPs that fall into footprint region (BED file from [this website](http://genome.grid.wayne.edu/centisnps/))
  2. Calculate the score change for each motif-SNP pair
* **Output**: a file containing the following columns `motif    SNP_id    off_site    binding(ref/alt)    ref.LLR    alt.LLR    footprint_source    ref.EstPrior    alt.EstPrior`

## Strategy 2 -- Posterior Probability

* **Input**: motif files with prior probability threshold in motif score scale, SNP list, genome file, DNase-seq experiment files
* **Procedure**:
  1. For each motif-DNase-seq pair, extract sequence windows (window size is 100 as impelemented [here](https://github.com/piquelab/which_gen_vars/blob/master/src/runCentipedeOnAll.R#L31)) with prior probability greater than 0.1 (as described [here](https://github.com/piquelab/which_gen_vars/blob/master/Makefile#L147))
  2. Filter out regions that contain SNPs
  3. For each motif-DNase-seq pair, train a CENTIPEDE model using extracted sequences
  4. Evaluate the performance of CENTIPEDE model by comparing motif score and posterior probability (see [here](https://github.com/piquelab/which_gen_vars/blob/master/src/runCentipedeOnAll.R#L60-L79))
  5. Run CENTIPEDE model on SNP list with the same window size and keep all SNPs that have posterior probability greater than 0.99 in either reference or alternative (see detail [here](https://github.com/piquelab/which_gen_vars/blob/master/src/runCentipedeOnAll.R#L141-L146))
* **Output**: a file containing the following columns `motif    SNP_id    off_site    binding(ref/alt)    ref.Posterior    alt.Posterior    ref.Prior    alt.Prior    footprint_source`

This module performs the same analysis as `strategy2/`
But it is designed for BAM inputs instead of FASTQ inputs

## Setup computing environment

```
source /project2/xinhe/yanyul/setup.sh
snmk  # activate a conda environment
module load java
```

## Steps in brief

This module takes an DNase-seq (or ATAC-seq) data set (a BAM file) and outputs the footprint of the transcription factor (a list of TF binding motif is required). For your reference, the details descriptions for each rule (in the order of running) in `Snakefile` are as follow.

### Pre-processing BAM file

This sub-module can take either pre-sorted BAM file or sorted BAM file (in `config.yaml` use `sort_bam: True` if it is sorted).

* For unsorted BAM:
    - `rule align_bam_sort`: Sort the input BAM file
    - `rule align_bam_index`: Index the sorted BAM
    - `rule remove_duplicated_reads`: Remove the duplicated reads from BAM and obtain a cleaned BAM
* For sorted BAM:
    - `rule remove_duplicated_reads_from_input`: Remove the duplicated reads directly from input BAM and obtain a cleaned BAM (**CAUTION: it haven't been tested**)

* `rule align_cleaned_bam_sort_index`: Sort and index the cleaned BAM

### Obtaining motif score

This sub-module extracts all candidate motif binding regions genome-wide. The definition of candidate motif binding region is the region such that:

* Motif binding p-value (computed by `fimo`) is more significant than a cutoff (in `config.yaml`, set it at `active_motifs: fimo_threshold: [some cutoff]`)
* One of the motif in the provided database has decent binding prior on this region

The complication comes from the second condition. Since the motifs (at `../recalibratedMotifs`) used is from [Which Genetics Variants in DNase-Seq Footprints Are More Likely to Alter Binding?](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005875). It also has a way to obtain prior other than motif score (*i.e.* log odds ratio). So, in this module, **prior > 0.1 is used as the cutoff to call candidate region for motif binding**

(**CAUTION: there is no way to use other cutoff but it can be added in the future**)

* `rule scan_motif_score_in_genome_prepare`: Convert motif file into MEME format (**CAUTION: only works for motifs at `../recalibratedMotifs`. Other format can be added in the future**)
* `rule scan_motif_score_in_genome`: Takes a list of motifs provided in MEME format and output all candidate region with TF binding p-value less than `fimo_threshold`
* `rule scan_motif_score_in_genome_gzip`: Compress the output file of `rule scan_motif_score_in_genome`
* `rule scan_post_filtering`: Filter out candidate region with prior less than 0.1

### Preparing input file for CENTIPEDE

This sub-module summarizes motif scores and the count of read's 5' end around candidate region. It simply takes the candidate regions output from previous sub-module and extract the count of the region and its flanking region (in `config.yaml`, set up the window size in `centipede: window_size: [some_number]`, 200 is recommended) from the BAM file.

* `rule scan_motif_score_in_genome_extend`: Extends the candidate region (*i.e.* the region size is just the length of motif) to the size of flanking window (motif length plus the number set in `config.yaml`)
* `rule pile_prepare`: Extracts reads that lie in extended candidate region in BAM file
* `rule pile_do_pile`: Obtains 5' end count for forward strand and backward strand respectively
* `rule pile_postprocess_extract`: Formats the output of previous step in good shape
* `rule pile_postprocess_formatting`: Some more formatting

### Training CENTIPEDE

This sub-module runs CENTIPEDE model for each motif-BAM pair. In brief, for each pair, the CENTIPEDE model is trained and posterior of binding is saved in outputs. It also does some plotting on which the HTML report depends.

* `rule train_centipede`

### Generating report about training

This sub-module generates report for each BAM data set which contains the training results of all motifs on that data set.

* `rule extract_good_of_fit`: Extracts the p-value of the correlation between motif score and posterior for each CENTIPEDE run. The idea is that if the correlation is not significant, what is being picked is not in good shape. Also, the correlation coefficient `rho` is also extracted
* `rule summary_rmd` and `rule summary_html`: Generates the RMD and HTML for report

### Obtaining footprint

The sub-module extracts all footprint region (candidate region of motif binding) with decent posterior for each BAM data set. Importantly, only motif with significant correlation p-value (in `config.yaml`, set up the p-value cut off at `active_motifs: output_threshold: [some cut off]`, 6e-7 is recommended). Also, **only regions with posterior greater than 0.9 are reported**.

* `rule extract_significant_motifs`: Extracts significant motifs based on correlation p-value cut off
* `rule combine_prepare`: Prepares the input for next rule
* `rule combine`: Obtains the footprint regions along with the posterior
* `rule postprocess`: Formatting

## How to use

The pipeline is glued by `snakemake` where `config.yaml` is required to provide the information of required input files along with the name of output files. An example `config.yaml` is `config.test.yaml` and all output files (one training report per data set and one footprint BED file per data set) generated by test run are at `output/recalMotifs`, `plots/recalMotifs`, and `summary/recalMotifs`. The test run is

```
snakemake --configfile config.test.yaml
```

The following describes how `config.yaml` should be set up

### Experiments

This part specifies which data sets should be analyzed (called footprint from).

```
sort_bam: False
experiments:
  my_test1:
    bam: 'path_to_bam1'
  my_test2:
    bam: 'path_to_bam2'
```

`mytest1` and `mytest2` determines the part of the name in output file. For sorted BAM, set `sort_bam: True`.

### Motifs

This part specifies which motifs are considered. To work on a set of motifs, it should be placed at a folder with motif files only and the suffix should be `.pwm`.

```
motifs:  # one motif database per config flie (the ones other than the first will be ignored)
  my_motif_db: 'path_to_motif_folder/'
```

Note that only one motif set is required, others will be discarded.

### Genome assembly

This part specifies which genome assembly is used.

```
genome_assembly: # only one genome assembly as well
  my_genome:
    fasta: 'path_to_fasta'
    size: 'path_to_chrom_sizes'
```

Note that only one genome assembly is required and others will be discarded.

### Motif parameters

```
active_motifs:
  fimo_threshold: 1e-4
  output_threshold: 6e-7
```

Please refer to [Obtaining motif score](#obtaining-motif-score) and [Obtaining footprint](#obtaining-footprint) for the usage of these two parameters

### Centipede parameters

```
centipede:
  script: 'path_to_fitCentipedeV2.R'  # at https://github.com/piquelab/which_gen_vars/blob/master/src/fitCentipedeV2.R
  window_size: 200
```

If you want to use `fitCentipede` from `CENTIPEDE` R package, set `centipede: script: 'NO'`. If instead, `fitCentipede3` is to be used, download `https://github.com/piquelab/which_gen_vars/blob/master/src/fitCentipedeV2.R` and put the path there. For `window_size`, it determines the width of the flanking window for modeling footprint. 200 bp is a reasonable choice but it can be changed for your purpose.

### For debugging

```
test_motifs: 'M00001,M00008'  # only for debugging
debug: True
```

Often, the motif database contains a huge number of motifs. To run on a small subset of these motifs, fill in the names of the motifs to be tested in `test_motifs` as string separated by ',' and set `debug: True`. Whenever `debug: False`, all motifs will be run on.

### Other parameters

```
## DO NOT CHANGE IF UNNECESSARY
alignment_params:
  ncpus: 1
  clean_cmd: 'java -jar /project2/xinhe/yanyul/softwares/picard.jar \
  MarkDuplicates \
  I=\{input\} \
  O=\{output\} \
  M=\{output\}.metric.txt \
  REMOVE_DUPLICATES=true'
## END
```

Typically, do not change this part.

## Output

Using the `config.yaml`

```
sort_bam: False
experiments:
  my_test1:
    bam: 'path_to_bam1'
  my_test2:
    bam: 'path_to_bam2'
motifs:  # one motif database per config flie (the ones other than the first will be ignored)
  my_motif_db: 'path_to_motif_folder/'
genome_assembly: # only one genome assembly as well
  my_genome:
    fasta: 'path_to_fasta'
    size: 'path_to_chrom_sizes'
active_motifs:
  fimo_threshold: 1e-4
  output_threshold: 6e-7
centipede:
  script: 'path_to_fitCentipedeV2.R'  # at https://github.com/piquelab/which_gen_vars/blob/master/src/fitCentipedeV2.R
  window_size: 200
test_motifs: 'M00001,M00008'  # only for debugging
debug: True
## DO NOT CHANGE IF UNNECESSARY
alignment_params:
  ncpus: 1
  clean_cmd: 'java -jar /project2/xinhe/yanyul/softwares/picard.jar \
  MarkDuplicates \
  I=\{input\} \
  O=\{output\} \
  M=\{output\}.metric.txt \
  REMOVE_DUPLICATES=true'
## END
```

By running `snakemake --configfile config.yaml`, the output files are

```
# BED file indicating footprint location
output/my_motif_db/my_test1__genome_my_genome__window.200.final.bed.gz
output/my_motif_db/my_test2__genome_my_genome__window.200.final.bed.gz
# HTML file reporting the training
summary/my_motif_db/my_test1__genome_my_genome__window.200.report.html
summary/my_motif_db/my_test2__genome_my_genome__window.200.report.html
plots/my_motif_db/*__my_test1__genome_my_genome__window.200_cutsite.png
plots/my_motif_db/*__my_test1__genome_my_genome__window.200_footprint.png
plots/my_motif_db/*__my_test2__genome_my_genome__window.200_cutsite.png
plots/my_motif_db/*__my_test2__genome_my_genome__window.200_footprint.png
```

The `plots/` contains the figures which are necessary for HTML to work properly. So, you need to download `plots/` along with HTML file.

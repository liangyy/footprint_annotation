This module performs the same analysis as `strategy2/``
But it is designed for BAM inputs instead of FASTQ inputs

## Setup computing environment

source /project2/xinhe/yanyul/setup.sh
snmk  # activate a conda environment
module load java

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

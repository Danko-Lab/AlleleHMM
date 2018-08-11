# AlleleHMM
AlleleHMM takes allele-specific read counts file as input. By default, it runs with 9 values of tao (0.1, 0.01,..., 1e-09), and output bed files of predicted states for each value of tao and input file.


## Input files

AlleleHMM takes allele-specific read counts file in the following formats, please notes that:
+ Please use tab delimited text file
+ Must have header at the first line and only the first line.
+ SNP position (snppos) must be sorted according to genomic location. 
+ Please see full example of counts_hmm.txt in input_file_exmaples folder.

```````
chrm    snppos  mat_allele_count        pat_allele_count        total_reads_count       state
1       565006  0       17      17      P
1       565286  46      0       46      M
1       565406  37      0       37      M
1       565419  31      0       31      M
1       565591  27      0       27      M
1       566573  0       2       2       S
1       568214  0       6       6       P
1       569094  93      0       93      M
1       569933  0       2       2       S
```````


## Usage

+ For strand-specific data such as PRO-seq, please prepare two files.
  * counts_plus_hmm.txt: allele-specific read counts file generated from plus strand
  * counts_minus_hmm.txt: allele-specific read counts file generated from minus strand
```````
python AlleleHMM.py -p counts_plus_hmm.txt -m counts_minus_hmm.txt
```````
+ For non-strand-specific data such as ChIP-seq, please prepare one file counts_hmm.txt.
```````
python AlleleHMM.py -i counts_hmm.txt
```````
  

## Output files
+ counts_hmm_regions_t1e-05.bed, counts_plus_hmm_regions_t1e-05.bed, counts_minus_hmm_regions_t1e-05.bed: genomic regions with predicted allele-specificities.
+ counts_hmm_t=1e-05_parameters.txt: Optimized transition probability and emission probablity using AlleleHMM. They are used to predict the allele-specificities of genomic regions.

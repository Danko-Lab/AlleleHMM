# AlleleHMM



## Input files

AlleleHMM takes allele-specific read counts file in the following formats, please notes that:
+ Please use tab delimited text file
+ Must have header at the first line and only the first line.
+ chrm must be interger.
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
+ For strand-specific data such as PRO-seq, please prepared two files.
  * counts_plus_hmm.txt: allele-specific read counts file generated from plus strand
  * counts_minus_hmm.txt: allele-specific read counts file generated from minus strand
+ autosome_num: Human is 22, mouse is 19 

```````
python AlleleHMM.py prefix counts_plus_hmm.txt counts_minus_hmm.txt autosome_num
```````
The following is an example for human PRO-seq data
```````
python AlleleHMM.py human_GM12878 counts_plus_hmm.txt counts_minus_hmm.txt 22
```````

For non-strand-specific data such as ChIP-seq, please prepared one file.


```````
python AlleleHMM.py counts_hmm.txt counts_plus_hmm.txt counts_minus_hmm.txt 22
```````
  

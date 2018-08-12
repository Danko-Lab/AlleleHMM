# AlleleHMM
The key goal of AlleleHMM is to identify allele-specific blocks of signal in distributed functional genomic data assuming that contiguous genomic regions share correlated allele-specific events. We developed a HMM that represents allelic bias in a distributed genomic mark using three hidden states: symmetric (S) distribution of the mark from both parental alleles (which shows no allelic bias), and maternally- (M) or paternally-biased (P) regions. AlleleHMM takes read counts corresponding to each allele as input. AlleleHMM uses these allele-specific read counts to set the parameters of the HMM using Baum Welch expectation maximization. The Viterbi algorithm is then used to identify the most likely hidden states through the data, resulting in a series of candidate blocks of signal with allelic bias.

AlleleHMM.py identify candidate allele-specific blocks using 9 values of tao (1E-01, 1E-02, ...,1E-09) by default. User can assign a specific tao using -t option.

## Usage
```````
python AlleleHMM.py [options]

options:

To get help:
-h,                    Show this brief help menu.

Required options:
For non-strand-specific data such as ChIP-seq:
-i, --input_hmm=PATH   Path to the non-strnad-specific, allele-specific read counts file (counts_hmm.txt)

For strand-specific data such as PRO-seq:
-p, --input_plus_hmm=PATH    Path to the plus-strand allele-specific read counts file (counts_plus_hmm.txt)
-m, --input_minus_hmm=PATH   Path to the minus-strand allele-specific read counts file (counts_minus_hmm.txt)

Optional operations:
-o, --output_prefix=STR      prefix for the output file. default=AlleleHMM_output
-t, --tao=FLOAT   AlleleHMM identify allele-specific blocks using 9 values of t (1E-01, 1E-02, ...,1E-09) by default.
                  User can assign a specific tao for the calculation.
```````

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


## Output files
+ AlleleHMM_output_[STRAND]_regions_t[TAO].bed: candidate blocks of signal with allelic bias in bed file format.
+ AlleleHMM_output_t=[TAO]_parameters.txt: Optimized transition and emission probablities using AlleleHMM.

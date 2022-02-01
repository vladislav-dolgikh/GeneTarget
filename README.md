# GeneTarget

Requirements
------------

PYTHON 3
numpy: `pip3 install numpy`

pandas: `pip3 install pandas`

Biopython: `pip3 install bio`

Usage
--------
The command `python3 gene_target.py -h` return:
```
  usage: gene_target.py [-h] chipseq_peaks rnaseq_table motif_recognition file_output file_PWM threshold chipseq_fasta 
  
  positional arguments:
  
  `chipseq_peaks` BED file contains ChIP-seq peaks
  
  `rnaseq_table` Table contains log2FC and p-adj values from RNA-seq processing experiment
  
  `motif_recognition` A boolean function determines whether the motif will be recognized (0 - if not, 1 - if yes)
  
  `file_output` Name of output file
  
  optional arguments:
  
  `file_PWM` File contains position weight matrix
  
  `threshold` Threshold value for position weight matrix

  `chipseq_fasta` FASTA file contains sequences of ChIP-seq peaks  
 
   `-h, --help` show this help message and exit  
   
 ```
    
 Example run without motif recognition: `python3 binner.py pwm_1.txt 0.91 arabidopsis_promoters.fas arabidopsis_promoters_coordinates.bed DEGs_matrix.txt pwm_1_output.txt`

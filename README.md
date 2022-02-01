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

  usage: gene_target.py [-h] chip_seq_peaks rnaseq_table motif_recognition file_output file_PWM threshold chip_seq_fasta 
  
  positional arguments:
  
  `chip_seq_peaks` BED file contains ChIP-seq peaks
  
  `rnaseq_table` Table contains log2FC and p-adj values from RNA-seq processing experiment
  
  `motif_recognition` A boolean function determines whether the motif will be recognized (0 - if not, 1 - if yes)
  
  `file_output` Name of output file
  
  optional arguments:
  
  `file_PWM` File contains position weight matrix
  
  `threshold` Threshold value for position weight matrix

  `chip_seq_fasta` FASTA file contains sequences of ChIP-seq peaks  
 
   `-h, --help` show this help message and exit  
   
The program can work in two modes: with and without motive recognition. If the motif recognition function is active, then the user must provide the positional weight matrix, threshold value, and ChIP-seq peak sequences in FASTA format.
 
 Example run without motif recognition: `python3 gene_target.py ein3_chipseq_4h.bed ethylene_rnaseq_full.txt 0 target_genes.txt`
 Example run with motif recognition: `python3 gene_target.py ein3_chipseq_4h.bed ethylene_rnaseq_full.txt 1 target_genes.txt pwm_1.txt 0.91 ein3_chipseq_4h.fasta`
 
 File formats
--------

### 1. ChiP-seq peaks (chip_seq_peaks argument)

A simple BED3 file contains columns with chromosome number, coordinates of start and end of the peak.

|Column name|Information|
|---|---|
|Chromosome|Chromosome number|
|Start|Peak start coordinate|
|End|Peak end coordinate|

#### Coordinate example

 ```
1	2131	3631
1	8737	10237
1	13714	15214
1	21646	23146
1	27000	28500

### Table of RNA-seq output (rnaseq_table argument)

Table contains gene ID's and log2FC and p-adj values for each RNA-seq experiment. The number of experiments can be any.














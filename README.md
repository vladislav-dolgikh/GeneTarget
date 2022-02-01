# GeneTarget

Requirements
------------

PYTHON 3

numpy: `pip3 install numpy`

pandas: `pip3 install pandas`

Biopython: `pip3 install bio`

Usage
--------
The command `python3 gene_target.py -h` returns:

  `usage: gene_target.py [-h] chip_seq_peaks rnaseq_table motif_recognition file_output file_PWM threshold chip_seq_fasta`
  
####  positional arguments:
  
  `chip_seq_peaks` BED file containing ChIP-seq peaks
  
  `rnaseq_table` Table containing log2FC and p-adj values from RNA-seq processing experiment
  
  `motif_recognition` A boolean function determines whether the motif will be recognized (0 - if not, 1 - if yes)
  
  `file_output` Name of output file
  
####  optional arguments:
  
  `file_PWM` File containing position weight matrix
  
  `threshold` Threshold value for position weight matrix

  `chip_seq_fasta` FASTA file containing sequences of ChIP-seq peaks  
 
   `-h, --help` show this help message and exit  
   
The program can work in two modes: with and without motif recognition. If the motif recognition function is active, user must provide the positional weight matrix, threshold value, and ChIP-seq peak sequences in FASTA format.
 
#### Example run without motif recognition: 
 `python3 gene_target.py ein3_chipseq_4h.bed ethylene_rnaseq_full.txt 0 target_genes.txt`
 
#### Example run with motif recognition: 
 `python3 gene_target.py ein3_chipseq_4h.bed ethylene_rnaseq_full.txt 1 target_genes.txt pwm_1.txt 0.91 ein3_chipseq_4h.fasta`
 
 File formats
--------

### 1. ChIP-seq peaks (chip_seq_peaks argument)

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
 ```

### 2. Table of RNA-seq output (rnaseq_table argument)

Table contains gene ID's and log2FC and p-adj values for each RNA-seq experiment. The number of experiments can be any.

|Column name|Information|
|---|---|
|AGI ID|Name of gene in AGI format|
|log2FC|Fold change value for gene|
|p-adj|Adjusted p-value for gene|

#### Coordinate example

 ```
AT1G01200	0.798468892	0.999963365	0.262868659	0.858116364
AT1G01210	0.24716958	0.999963365	0.299284788	0.795936156
AT1G01220	-0.130354152	0.999963365	0.243543128	0.811947466
AT1G01225	-0.398548781	0.999963365	0.080885801	0.980423396
AT1G01230	0.094237385	0.999963365	-0.401992625	0.372281038
 ```

### 3. Positional weight matrix (file_pwm argument)

A simple position weight matrix, where columns correspond to letters of genetic alphabet (A, C, G, T), and rows - to match weight in positions.

#### PWM example

|A|C|G|T|
|---|---|---|---|
|-1.62|-7.96|1.53|-0.35|
|0.33|-1.68|1.27|-7.96|
|-0.21|-7.96|1.49|-1.68|
|-7.96|-7.96|1.78|-0.88|
|-1.57|-1.68|1.74|-7.96|

### 4. Sequences of ChIP-seq (chip_seq_fasta argument)

A simple file in FASTA format. The header should have the following format: `chromosome_number-start-end`

#### ChIP-seq sequence example

 ```
>1-237697-237961
AAAAAGATTAGAGAGGAAAGATGGAGAAAAAATGGAGGAAGGTA...
 ```

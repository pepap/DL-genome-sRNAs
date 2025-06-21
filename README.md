# DL-genome-sRNAs

## miRNA annotation
The PERL program miREAP was used for miRNA annotation : [miREAP](https://github.com/liqb/mireap)

All scripts used for annotation can be found here :

- [miRNAs](https://github.com/pepap/DL-genome-sRNAs/tree/main/miRNAs)

### 01. Input file preparation : merge & map
Input sequencing files are required in a specific collapsed-fasta format, with more details available in [miREAP](https://github.com/liqb/mireap). [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/download.html) software was used to format the files. The collapsed-fasta files were then mapped using [STAR](https://github.com/alexdobin/STAR/releases/tag/2.7.7a) software. All of these steps are summarized in the BASH script `01.PREP-INPUT--map.bash`.

### 02. Extract input `smrna.fa` and `map.txt`
The required input files for the miREAP PERL script are extracted in R. All the steps are summarized in the BASH script `02.01.FORMAT-INPUT--run.bash`.

### 03. miREAP run
Custom miREAP run : `03.MIREAP--run.bash`. The resulting GFF3 file was manually curated.

## piRNA annotation
The R library PiCB was used for the identification of piRNA clusters : [PiCB](https://github.com/HaaseLab/PICB)

You can find all used scripts here :

- [piRNAs](https://github.com/pepap/DL-genome-sRNAs/tree/main/piRNAs)

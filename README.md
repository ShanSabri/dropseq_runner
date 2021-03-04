# dropseq_runner
[![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/master) 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-ff69b4.svg)](https://www.gnu.org/licenses/gpl-3.0)

`dropseq_runner` is a collection of scripts used to process raw single cell sequencing data from the Drop-sea platform to a digital gene expression matrix (DGE). 

This repository was used to generate DGEs for the analysis in:

```bash
@inproceedings{CARRARO2021,
Author = {Gianni Carraro and Justin Langerman and Shan Sabri and Zareeb Lorenzana and 
Arunima Purkayastha and Bindu Konda and Cody J. Aros and Ben A. Calvert and 
Aleks Szymaniak and Emily Wilson and Michael Mulligan and Priyanka Bhatt and 
Preethi Vijayaraj and Changfu Yao and David W. Shia and Edo Israely and 
Tammy M. Rickabaugh and Martin Mense and Scott H. Randell and Eszter K. Vladar and
Amy L. Ryan and Kathrin Plath and John Mahoney and Barry R. Stripp and Brigitte N. Gomperts},
Title = {Transcriptional analysis of Cystic Fibrosis airways at single cell resolution reveals altered epithelial cell states and composition},
Year = {2021}
}
```

## Usage

Raw sequencing files in QSEQ format from Illumina needs to be merged into one file in the input format to `src/0_demux.py`. For example,


```bash
zcat /path/to/s_2_1_*_qseq.txt.gz | gzip -c > s_2_1.qseq.txt.gz
zcat /path/to/s_2_2_*_qseq.txt.gz | gzip -c > s_2_2.qseq.txt.gz
zcat /path/to/s_2_3_*_qseq.txt.gz | gzip -c > s_2_3.qseq.txt.gz
```

Latter scripts feed from the demultiplexing output and are numerically ordered by procedure. All dependencies are contained within `beds` and `DST` but paths will need to be updated within each script, typically within the header. 


## Contributing
This repo is not being maintained though pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html) - Free software!

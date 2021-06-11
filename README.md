# RNA-3Dgenome code repository

Repository of codes used for the paper "Three-dimensional organization of chromatin associated RNAs and their role in chromatin architecture in human cells" ([pre-print](https://www.biorxiv.org/content/10.1101/2021.06.10.447969v1)).

## iMARGI

### Raw data processing

iMARGI data were processed using the publicly available [iMARGI pipeline](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker) to obtain final `BEDPE` files used in the downstream analysis.

### Data analysis and visualization

- [`make_imargi_contact_matrix.r`](./iMARGI/make_imargi_contact_matrix.r) is used to create contact matrices at several resolutions from the BEDPE file.
- [`MARGI_heatmap.r`](./iMARGI/MARGI_heatmap.r) is used to plot contact heatmaps from contact matrices, as well as other iMARGI-related genomic tracks.
- [`margi_compartment.r`](./iMARGI/margi_compartment.r) is used for iMARGI analysis on Hi-C A/B compartments.
- [`margi_tads.r`](./iMARGI/margi_tads.r) is used for iMARGI analysis on Hi-C topologically associating domains (TADs).
- [`margi_loops.r`](./iMARGI/margi_loops.r) is used for iMARGI analysis on Hi-C chromatin loops.
- Irene add your scripts here (script files to be included under the folder `iMARGI`)
- 

### caRNA domains

Zhijie please add the link to your repository and additional codes here (script files to be included under the folder `iMARGI`).


## Hi-C

### Raw data processing

Hi-C data were processed using the publicly available [4DN Hi-C processing pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline) to obtain final `.hic` files used in the downstream analysis.


### Data analysis and visualization

Final `.hic` files were analyzed with the software [Juicer](https://github.com/aidenlab/juicer), which was used to extract contact data as well as call A/B compartments, TADs and loops.

- [`dump_contact_matrix.sh`](./HiC/dump_contact_matrix.sh) is used to extract contact matrices at different resolutions.
- [`HiC_heatmap.r`](./HiC/HiC_heatmap.r) is used to plot contact heatmaps from contact matrices.
- [`compartments.sh`](./HiC/compartments.sh) is used to call A/B compartments.
- [`compartment_analysis.r`](./HiC/compartment_analysis.r) is used for analysis on A/B compartments.
- [`TADs.sh`](./HiC/TADs.sh) is used to call TADs.
- [`TAD_analysis.r`](./HiC/TAD_analysis.r) is used for analysis on topologically associating domains (TADs).
- [`loops.sh`](./HiC/loops.sh) is used to call loops and perform Aggregate Peak Analysis (APA).
- [`loop_analysis.r`](./HiC/loop_analysis.r) is used for analysis on chromatin loops.

## PLAC-seq

Zhijie please fill in this part.
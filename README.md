

NanopoReaTA - Nanopore Real-Time Transcriptional Analysis Tool
==================================================


**NanopoReaTA** is an R shiny application that integrates both preprocessing and downstream analysis pipelines for RNA sequencing data from [Oxford Nanopore Technologies (ONT)](https://nanoporetech.com/) into a user-friendly interface. This repository integrates NanopoReaTA into Epi2Me and operates on a browser based Graphical user interface- NanopoReaTA focuses on the analysis of (direct) cDNA and RNA-sequencing (cDNA, DRS) reads and guides you through the different steps up to final visualizations of results from i.e. differential expression or gene body coverage. Furthermore, NanopoReaTa can be run in real-time right after starting a run via MinKNOW, the sequencing application of ONT. 


**Currently available analysis modules:**
1. [Run Overview](#run-overview) - Experiment statistics over time
2. [Gene-wise analysis](#gene-wise-analysis) - Gene-wise analysis of expression (gene counts, gene body coverage)
3. [Differential expression analysis](#differential-expression-analysis) - Differential expression and/or usage analysis of genes (DGE) and transcripts (DTE + DTU)


# Installation
## Requirements
Hardware | Min. number
 :---: |  :---:
RAM | 64GB 
Threads | 12 
**Biological input** | 
Total number of samples | 4
Number of conditions | 2
Min. number of samples per condition | 2

## User guide
To use NanopoReaTA [Epi2Me](https://labs.epi2me.io/downloads/), [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/) must be installed.
On Windows we recommend to install [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) and subsequently install Epi2Me, Nextflow and Docker within the command line.

To download NanopoReTA, open Epi2Me and navigate to Launch and press the button Import Workflow. A pop-up window will appear in which you should copy the following link: "https://github.com/stegiopast/wf-nanoporeata"
Press download and the Workflow should be integrated in Epi2Me. 


## Usage
Before running/exploiting real experiments with NanopoReaTA, we highly recommend to test the app, first (see [Testing](#testing) for more information). NanopoReaTA operates with a backend preprocessing pipeline based on [nextflow](https://www.nextflow.io/) [3] and multiple R and python based scripts for downstream analyses. All results are visualized within the [R shiny](https://shiny.rstudio.com/) [2] based frontend on the browser, which can be reached by clicking the link in the Epi2Me report session. The link in the reports sections shows up when for every barcode/sample in the metadata file at least a single dataset chunk has been processed.

#### Example metadata file
For usage you need to create a metadata file with necessary information. If you have barcoded sample your samplnames must be "barcode01, barcode02 etc."

 Samples | Conditions | Replicates
 :---: | :---: | :---: 
 Sample1 | Cond1 | R1 
 Sample2 | Cond2 | R2 
 Sample3 | Cond1 | R1 
 Sample4 | Cond2 | R2 

### Inputs via Epi2ME

 Parameter | Datatype
 :---: | :---: 
 Path to sequencing directory | path 
 Path to a metadata tsv | path 
 Path to Reference genome fasta (Gencode) | path
 Path to Reference transcriptome fasta (Gencode) | path 
 GTF annotation (Gencode) | path 
 BED annotation (RseQC) | path 
 Number of threads | integer
 Run preprocessing | bool 
 Barcoded | bool  

After all configurations are set, the configurations will be saved as config.yaml in the defined output folder.
 
#### Settings overview
The input configurations can be finally checked by the user. If the parameters are correct, the user can start the preprocessing by clicking the **Start** button. Otherwise the user can rearrange the settings by going back to the configuration tab.


### Run Overview
The Run Overview tab shows the number of mapped reads and gene counts and visualizes the sample- and group-wise read length distribution and gene expression variability per preprocessing iteration. Additionally, the time each tool needs in each iteration is shown. All information is constantly updating when preprocessing is running.

#### Number of observations
The table in this tab shows the number of mapped genes (*minimap2* [2]), gene counts (*featureCounts* [3]) and transcriptome counts (*salmon* [6]). The counts are provided for each sample, respectively.

#### Read length distribution 
One can see the read length distributions for respective samples and conditions. The read length information is extracted directly from the fastq files (MinKNOW-defined passed reads only).

#### Gene expression variability

On the left side the number of genes detected is plotted per iteration for samples and selected conditions, respectively. The information is extracted from the output count table of *featureCounts*. 

On the right side the deviation of relative gene abundancy compared to the last iteration is plotted. This is a measure for the change of gene abundancy variability within a single sample. The latter allows an assumption whether relative abundancies have stabilized throughout the ongoing sequencing.  


### Gene-wise analysis

In the Gene-wise anaylsis tab, one is able to explore the expression levels and the gene body coverage of particular genes of interest. Be aware that at least two samples per condition have to be considered in order to use this functionality.  

<p align="center"><img src="Gifs/nanoporeata_supp2_fig2.png"  width="100%"></p>


#### Gene counts

The table on the lefthand side lists all the genes annotated in the loaded GTF file. One can search and select several genes of interest via click on the table entry. Once a gene is selected it will occur on the table at the right hand side. By clicking the submit genes button, the analysis will start. A median of ratio normalization via DESeq2 [4] will be performed and the user can plot the raw and normalized counts per condition as Dot-, Violin- or Boxplot.

<p align="center"><img src="Gifs/Selected_Genes.gif"  width="80%"></p>

#### Gene Body coverage

Here, one gene can be selected for gene body coverage analysis each time. The gene selection functions similar as in [Gene counts](#gene-counts). After the gene selection is submitted, the percentage of coverage for a gene divided into 100 percentiles is shown sample- and group-wise (=mean). The calculation is based on the RSeQC script for gene body coverage analysis (https://rseqc.sourceforge.net/) [7].


### Differential Expression Analysis
In the Differential Expression Analysis tab, the user can run three different analyses: Differential Gene Expression (DGE), Differential Transcript Expression (DTE) and Differential Transcript Usage (DTU) by clicking the respective button. Note that these analyses do not update automatically when processing will be started again and new data is generated. That means that after stopping the preprocessing pipeline again, the analyses buttons need to be pressed to analyse latest input files (like counts files). Once the analysis is completed, the user will be linked to the respective analysis output tab (may take a few minutes). 


#### Gene-level analysis (DGE with DESeq2[4])
Differential gene expression analysis will be performed and the following visualizations are shown:
- A table of all differentially expressed genes
- PCA analysis
- Volcano plot (DGE)
- Sample-2-Sample plot
- Heatmap of the top 20 differentially expressed genes based on p-adjusted

<p align="center"><img src="Gifs/PCA.gif"  width="80%"></p>
<p align="center"><img src="Gifs/Volcano Plots_DGE.gif"  width="85%"></p>

#### Transcript-level analysis
##### Differential Transcript Expression (DTE with DESeq2 [4])
Differential transcript expression analysis will be performed and the following visualizations are shown:
- A table of all differentially expressed transcripts
- PCA analysis
- Volcano plot (DTE)
- Sample-2-Sample plot
- Heatmap of the top 20 differentially expressed transcripts based on p-adjusted

<p align="center"><img src="Gifs/Volcano Plots_DTE.gif"  width="85%"></p>

##### Differential Transcript Usage (DTU with DRIMSeq [5] and DEXSeq [1])
Differential transcript usage analysis will be performed with DEXSeq and DRIMSeq.
The following visualizations are shown: 
- Tables of all differentially used transcripts' analysis results (Results of DEXSeq and DRIMSeq)
- General: To allow a general DTU overview, the log2FoldChange values from DEXSeq are plotted against the adjusted p-value (Volcano plot)
- Gene specific: One can select a gene to show the transcripts abundances within a gene of interest as boxplots per condition based on DRIMSeq's output.

<p align="center"><img src="Gifs/DTU_Boxplot.gif"  width="85%"></p>


## Test data
We provide a dataset of cDNA extracted from 2 samples of HEK293 and 2 samples of HeLa cells.

Test data is available in a seafile repository https://seafile.rlp.net/d/7a99b8b210e44eb9b70a/. The output folder structure of MinKnow is kept intact with this plattform, which enables users to directly test the functionality of NanopoReaTA.

Test data is additionally available on the ENA (European Nucleotide Archive) with the project number PRJEB61670: https://www.ebi.ac.uk/ena/browser/view/PRJEB61670.
Please note that one will have to reconstruct the folder structure of the MinKnow output using barcoded samples when dowloading files from ENA to use NanopoReaTA properly.



Examples: ***Experiment_folder**/Sample_folder/Identifier/fastq_pass/barcodeXY/run_xyz_999.fastq* with (barcode01-barcode04). Barcode01 + barcode02 are HEK293 cDNA samples and barcode03-04 are HeLa cDNA samples.  

The metadata file can be found under https://github.com/AnWiercze/NanopoReaTA/blob/master/example_conf_files/example_metadata.txt.


## Publications

This is an derivative of the software published in Oxford Bioinformatics:

Anna Wierczeiko, Stefan Pastore, Stefan Mündnich, Anne M Busch, Vincent Dietrich, Mark Helm, Tamer Butto, Susanne Gerber, NanopoReaTA: a user-friendly tool for nanopore-seq real-time transcriptional analysis, Bioinformatics, Volume 39, Issue 8, August 2023, btad492, https://doi.org/10.1093/bioinformatics/btad492

An application paper of the original NanopoReaTA has been published in eLife:

Tamer Butto, Stefan Pastore, Max Müller, Kaushik Viswanathan Iyer, Marko Jörg, Julia Brechtel, Stefan Mündnich, Anna Wierczeiko, Kristina Friedland, Mark Helm, Marie-Luise Winz, Susanne Gerber
2024, Real-time transcriptomic profiling in distinct experimental conditions, eLife13:RP98768, https://doi.org/10.7554/eLife.98768.2


## References

[1] Anders, S., Reyes, A., & Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22(10), 2008–2017. https://doi.org/10.1101/GR.133744.111

[2] Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/BIOINFORMATICS/BTY191

[3] Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923–930. https://doi.org/10.1093/BIOINFORMATICS/BTT656

[4] Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 1–21. https://doi.org/10.1186/S13059-014-0550-8/FIGURES/9

[5] Nowicka M, Robinson MD (2016). “DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics [version 2; referees: 2 approved].” F1000Research, 5(1356). doi: 10.12688/f1000research.8900.2, https://f1000research.com/articles/5-1356/v2. 

[6] Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 2017 14:4, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

[7] Wang, L., Wang, S., & Li, W. (2012). RSeQC: quality control of RNA-seq experiments. Bioinformatics, 28(16), 2184–2185. https://doi.org/10.1093/BIOINFORMATICS/BTS356



## Contact

Please open an [issue](https://github.com/stegiopast/wf-nanoporeata/issues) if you encounter any issues/troubles. 
However, please go over the previous issues (including closed issues) before opening a new issue, as your same exact question might have been already answered previously. Thank you!

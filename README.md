# Nextflow workflow: Tracking Evolutionary Dynamics in E. coli Using RNA-Seq Variants

## Introduction 
### Background
Bacteria such as Escherichia coli provide an unparalleled model system for studying evolution due to their rapid reproduction rates, high mutation rates, and the relative simplicity of their genomes. In laboratory conditions, E. coli populations adapt quickly to environmental changes, revealing the mechanisms underlying microbial evolution, adaptation, and population structure. Time-series experiments, where bacterial populations are sampled at regular intervals, offer a powerful approach to track evolutionary trajectories in real-time. These experiments allow researchers to uncover how mutations arise, persist, or disappear under selective pressures, enabling insights into processes like adaptation, clonal expansion, and lineage divergence [1].

RNA sequencing (RNA-Seq) has revolutionized our ability to investigate microbial biology, offering both quantitative insights into gene expression and qualitative insights into sequence variation. While expression-based approaches can provide valuable insights into environmental responses, they often fail to capture the stable, heritable changes that define evolutionary divergence. Expression levels are transient, influenced by environmental conditions or regulatory dynamics, and do not necessarily reflect permanent genetic changes. In contrast, mutations identified from variant calls, such as single nucleotide polymorphisms (SNPs) and insertions or deletions (indels), represent the cumulative genetic changes that drive adaptation and divergence over time [2]. These variants are particularly valuable in experimental evolution studies because they reflect heritable changes in the genome that underpin adaptive processes. By transforming RNA-Seq data into variant call files (VCFs), researchers can shift the focus from transient gene expression patterns to stable genetic mutations, enabling a more precise reconstruction of evolutionary history [3].

Phylogenetic trees, constructed from variant data, serve as a critical tool in evolutionary biology by visually representing the genetic relationships among populations or time points. When applied to time-series E. coli experiments, phylogenetic trees provide a framework for tracking how bacterial populations evolve under experimental conditions, identifying which lineages diverge more rapidly and elucidating clonal expansion patterns [4]. Moreover, this approach highlights mutations that may be under selection, offering a direct link between genotype and phenotype [5]. For example, in antibiotic resistance studies, variant-based phylogenies can reveal how specific mutations arise and spread within a population subjected to drug pressure, shedding light on the mechanisms driving resistance evolution [6].

This study presents an approach for analyzing RNA-Seq data from a time-series E. coli experiment, focusing on extracting genetic variants and constructing phylogenetic trees. By doing so, we aim to capture the evolutionary dynamics of the population, identify key mutations driving adaptation, and provide a framework for understanding how microbial populations respond to environmental challenges. This approach is particularly relevant in microbial evolution studies, where insights into lineage relationships and mutational landscapes are critical for advancing our understanding of adaptation, population structure, and evolutionary processes [7].

### Purpose
This pipeline aims to construct phylogenetic trees based on mutational changes derived from RNA seq data to conduct an evolutionary analysis in time-based experiments. 

##  Usage
### Installation
Installing this pipeline requires [conda](https://docs.anaconda.com/miniconda/), [git](https://github.com/git-guides/install-git), and [nextflow](https://www.nextflow.io/docs/latest/install.html).

Navigate to the directory of your choice or create a new directory and then clone the repository by running the following command in a terminal.

```
git clone https://github.com/NikolayAlabi/BIOF501.git
```
Navigate into the repository and create the conda environment form the environment file, activate it. One final step will require you to ensure all of the additional scripts required by the workflow are executable. 

```
cd BIOF501
conda deactivate #in case you had any activated environments at the time
conda env create --file environment.yml -n BIOF501
conda activate BIOF501
chmod +x bin/*
```
Now you are all ready to run the workflow.

```
nextflow main.nf
```
## Pipeline Overview
This workflow was built using `Nextflow`, a useful tool to create reproducible and scalable data analyses. Example data have been provided in the data folder to demonstrate how this pipeline runs. 
### Input Data
The example input data consists of a standard fasta file of an _E. coli_ genome along with 3 sets of fastq paired end reads of RNA seq data from a time experiment. For more information on these samples, you can refer to this [SRA Experiment](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA295606) It is important to note that these fastq samples are already trimmed and filtered by using the clipped and filtered options available when download SRA RNAseq data. Moreover, these samples have been downsampled in order to demonstrate the pipeline on a small and quick dataset. The fastq files must be gunzipped and have a pair of files to represent each sample. 

This pipeline can be ran with custom data, showing its generizalibililty for other bigger projects or with different species. Currently this pipeline can only be ran for haploid organisms, however in the future it will be able to others. 
```
#ensure fastq naming pattern matches your own
nextflow main.nf --reads 'path/to/reads_directory/*_{1,2}.trim.sub.fastq" --reference "path/to/reference_genome/ecoli.fasta" --ploidy 1
```
### Pipeline In-Depth
![alt text](https://github.com/NikolayAlabi/BIOF501/blob/main/dag3.png)

A figure depicting the workflow is shown above. The workflow is broken up into 4 major steps. 
1. The alignment of reads against the genome. The input files of this pipeline include the e.coli genome which gets indexed as well as the pairs of reads which get aligned. 
2. The next step is variant calling on the alignments to identify mutations. The variant call files for each sample are then merged into one file and a reference genome is added to the file to aid with tree interpretation.
3. The SNP genotypes present in the VCF file are then used to create a matrix for phylogenetic analysis.
4. A phylogenetic tree is then constructed using IQTree and visualized.

For more specific details at each step, please refer to the documentation in the `main.nf` file.

### Expected Outputs
The key files produced in the `output/` folder are:
1. `merged_altered_variants.vcf` This is a standard vcf file that can be viewed in Excel. It should contain a column per sample as well as "Reference" column.
2. `merged_altered_variants.min1.phy` This is the matrix that is used to construct the phylogenetic tree.
3. `merged_altered_variants.min1.phy.iqtree` This is the main output from the IQTree and contains details of the parameters of the construction of the tree.
4. `tree_plot.png` This is a visualization of the tree made from `merged_altered_variants.min1.phy.tree`.


## References
1. Lenski, R. E. et al. (1991). Long-term experimental evolution in Escherichia coli. American Naturalist.
2. McGettigan, P. A. (2013). Transcriptomics in the RNA-seq era. Current Opinion in Chemical Biology.
3. Rocha, E. P. C. (2008). The organization of the bacterial genome. Annual Review of Genetics.
4. Felsenstein, J. (2004). Inferring Phylogenies. Sinauer Associates.
5. Lang, G. I., & Desai, M. M. (2014). The spectrum of adaptive mutations in experimental evolution. Genetics.
6. Toprak, E. et al. (2012). Evolutionary paths to antibiotic resistance under dynamically sustained drug selection. Nature Genetics.
7. Ochman, H., & Selander, R. K. (1984). Standard reference strains of Escherichia coli from natural populations. Journal of Bacteriology.

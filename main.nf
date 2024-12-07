#!/usr/bin/env nextflow

// Define input parameters
params.reads = "data/reads_fastq_trm_downsampled/*_{1,2}.trim.sub.fastq"  // Glob pattern for paired-end reads
params.reference = "data/reference_genome/ecoli.fasta"  // Reference genome
params.ploidy = 1 //the ploidy of the organism of study

// The following process performs bwa indexing on the reference genome. This is necessary for the subsequent alignment phase
process genome_index{
    input:
        path(genome)

    output:
        tuple path("${genome.getName()}"), path("${genome.getName()}.*")
    
    //We create an empty fasta file in order for the following step to know where the index files are stored.
    script:
    """
    bwa index -p ${genome.getName()} ${genome}
    touch "${genome.getName()}"
    """
}

// Process to run BWA for a pair of reads. This is then outputted to samtools to convert to a bam file, and then gets sorted
process alignReads {
    input:
    tuple path(genome), path(index_files)
    tuple val(sample_id), path(R1), path(R2)  // Sample ID, paired reads, and reference genome

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), emit: bam // BAM file named by sample ID

    script:
    """
    # Align paired-end reads to the reference genome
    bwa mem ${genome} ${R1} ${R2} | samtools view -S -b | samtools sort -o ${sample_id}.aligned.sorted.bam
    """
}

// This process takes Generates genotype likelihoods from a BAM file, then Calls variants (SNPs and indels) from the likelihoods 
// the variant calls are then filtered, lastly, the calls are normalized and the output is compressed
process variantCalling {
    input:
    path(genome)
    tuple val(sample_id), path(bamfile)
    val ploidy
    
    output:
    tuple path("${sample_id}_final_variants.vcf.gz"), path("${sample_id}_final_variants.vcf.gz.tbi")

    script:
    """
    bcftools mpileup -O b -f ${genome} ${bamfile} | bcftools call --ploidy ${ploidy} -m -v | vcfutils.pl varFilter | bcftools norm -f ${genome} -Oz -o ${sample_id}_final_variants.vcf.gz
    bcftools index -t ${sample_id}_final_variants.vcf.gz
    """
}

//this process merges all of the variant files we produced
process mergeVCF {
    input:
    path vcfFiles  // List of VCF files to merge
    path indices
    output:
    path "merged_variants.vcf"

    script:
    """
    bcftools merge -o merged_variants.vcf -O v ${vcfFiles.join(' ')}
    """
}

process alterVCF {
    publishDir "output", mode: "copy"
    input:
    path(merged_file)
    output:
    path "merged_altered_variants.vcf"
    script:
    """
    #first we run this awk script to add a reference genome sample to the calls
    #when we build the tree, it will help put in perspective how distant other samples are
    add_reference.awk ${merged_file} > temp_alter.vcf
    #for samples that do not have a variation 
    sed -e 's/.:.:./0:255,0,10/g' temp_alter.vcf > merged_altered_variants.vcf
    """
}

//this process creates a matrix from the vcf file that can then be used to build a phylogenetic tree
process phyGeneration {
    publishDir "output", mode: "copy"
    input:
    path(merged_altered_file)
    output:
    path "merged_altered_variants.min1.phy"
    script:
    """
    vcf2phylip.py -i ${merged_altered_file} -m 1 -o Reference
    """
}

process tree {
    publishDir "output", mode: "copy"
    input:
    path(phyfile)
    output:
    path("merged_altered_variants.min1.phy.iqtree")
    path("merged_altered_variants.min1.phy.treefile")
    script:
    """
    #we let IQ tree's ModelFinder determine the best model for the data with the ASC constraints since we are using variant calls
    iqtree -s ${phyfile} -m MFP+ASC -nt AUTO
    """
}

//this lets us visualize the actual tree
process visualize{
    publishDir "output", mode: "copy"
    input:
    path(iqtree)
    path(treefile)
    output:
    path("tree_plot.png")
    script:
    """
    visualize_tree.R ${treefile} tree_plot.png
    """

}

// Workflow
workflow{
    // Reference genome as input
    reference = file(params.reference)
    bwaindex = genome_index(reference)  
    //read pairs 
    read_pairs = Channel.fromFilePairs(params.reads, flat: true)
    bam = alignReads(bwaindex, read_pairs)
    variants = variantCalling(reference, bam.bam, params.ploidy)
    //this multimap allows us to split the output of the vcf files and their indices
    variants.multiMap { vcf,tbi ->
        vcf_files: vcf
        indices: tbi
    }
    .set { splits }

    merged_vcf = mergeVCF(splits.vcf_files.collect(),splits.indices.collect())
    altered_vcf = alterVCF(merged_vcf)
    phyfile = phyGeneration(altered_vcf)
    tree_out = tree(phyfile)
    visualize(tree_out)
}

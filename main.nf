#!/usr/bin/env nextflow

// Define input parameters
params.reads = "data/reads_fastq_trm_downsampled/*_{1,2}.trim.sub.fastq"  // Glob pattern for paired-end reads
params.reference = "data/reference_genome/ecoli.fasta"  // Reference genome
params.ploidy = 1

//Input Channel: Create a channel for paired-end reads
Channel
    .fromFilePairs(params.reads, flat: true)
    .set { read_pairs }

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

// Process to run BWA for each pair of reads. This is then outputted to 
process alignReads {
    input:
    tuple path(genome), path(index_files)
    tuple val(sample_id), path(R1), path(R2)  // Sample ID, paired reads, and reference genome

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), emit: bam
    path "${sample_id}.aligned.sorted.bam"  // BAM file named by sample ID

    script:
    """
    # Align paired-end reads to the reference genome
    bwa mem ${genome} ${R1} ${R2} | samtools view -S -b | samtools sort -o ${sample_id}.aligned.sorted.bam
    """
}

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
    input:
    path(merged_file)
    output:
    path "merged_altered_variants.vcf"
    script:
    """
    add_reference.awk ${merged_file} > temp_alter.vcf
    sed -e 's/.:.:./0:255,0,10/g' temp_alter.vcf > merged_altered_variants.vcf
    """
}

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
    iqtree -s ${phyfile} -m MFP+ASC -nt AUTO
    """
}

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
    // bwaindex.view()
    
    read_pairs = Channel.fromFilePairs(params.reads, flat: true)
    // // read_pairs.view()
    bam = alignReads(bwaindex, read_pairs)
    // // bam.view()
    variants = variantCalling(reference, bam.bam, params.ploidy)
    variants.multiMap { vcf,tbi ->
        vcf_files: vcf
        indices: tbi
    }
    .set { splits }

    // splits.vcf_files.collect().view()
    // splits.indices.collect().view()
    merged_vcf = mergeVCF(splits.vcf_files.collect(),splits.indices.collect())
    altered_vcf = alterVCF(merged_vcf)
    phyfile = phyGeneration(altered_vcf)
    tree_out = tree(phyfile)
    visualize(tree_out)

    //     // Process each pair of reads with the alignReads process
    // read_pairs
    //     .map { sample, reads -> tuple(sample, reads[0], reads[1]) }
    //     .set{ mapped_reads } // Feed paired reads and BWA index to the alignment process
    // mapped_reads.view { println "Mapped Reads: $it" }
    // bwaindex.view { println "BWA Index: $it" }


    // test = alignReads(bwaindex, read_pairs)
}

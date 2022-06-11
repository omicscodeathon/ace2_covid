#!/usr/bin/env nextflow

// Script parameters

params.samples = "/home/emurungi/gitau/marion/TNBC/TNBCtrimmed/*.fastq"
params.fastq = "/home/emurungi/gitau/marion/TNBC/raw"
params.zips = "/home/emurungi/gitau/marion/ACE/results/quality-reports"
params.trimmed = "/home/emurungi/gitau/marion/TNBC/TNBCtrimmed"
params.trim = "/home/emurungi/gitau/marion/TNBC/TNBCtrimmed"
params.ref = "/home/emurungi/gitau/marion/raw"
params.sam = "/home/emurungi/gitau/marion/TNBC/sam"
params.bam = "/home/emurungi/gitau/marion/TNBC/sam/bam"

params.outDir = "/home/emurungi/gitau/marion/ACE/results"
params.rscript = "/home/emurungi/gitau/marion/import_reads.R"
params.dgescript = "/home/emurungi/gitau/marion/dge.R"
params.enrichmentdownscript = "/home/emurungi/gitau/marion/Enrichment-down.R"
params.enrichmentupscript = "/home/emurungi/gitau/marion/Enrichment-up.R"


fastq_ch = Channel.fromPath(params.samples)
trim_ch = Channel.fromPath(params.trim)
quality_ch = Channel.fromPath(params.zips)
trimmed_ch = Channel.fromPath(params.trimmed)
sam_ch = Channel.fromPath(params.sam)
bam_ch = Channel.fromPath(params.bam)
ref_ch = Channel.fromPath(params.ref)

SAMPLE="SRR10729852"


process quality_control {

	input:
	file fastq from fastq_ch

	output:
	publishDir "${params.outDir}", mode: 'copy'

	script:
	"""
	fastqc ${fastq} -o ${params.outDir}/quality-reports/
	"""
}

process multiqc {

	input:
	file zips from quality_ch

	output:
	publishDir "${params.outDir}", mode: 'copy'

	script:
	"""
	module load multiqc
	multiqc ${params.zips} --outdir ${params.outDir}/quality-reports/

	"""
}

process remove_adaptors {

	output:
	publishDir "${params.outDir}", mode: 'copy'

	script:
	"""
	module load cutadapt
	cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m15 -o ${params.trim}/${SAMPLE}_1.trimmed.fastq -p ${params.trim}/${SAMPLE}_2.trimmed.fastq ${params.fastq}/${SAMPLE}_1.fastq ${params.fastq}/${SAMPLE}_2.fastq
	"""
}


process index_ref {

	output:
	publishDir "${params.outDir}", mode: 'copy'

	script:
	"""
	module load hisat2
	hisat2-build ${params.ref}/GCF_000001405.39_GRCh38.p13_genomic.fna ${params.outDir}/GCF_000001405.39_GRCh38.p13_genomic.fna_index_hisat2

	"""
}

process align_to_ref {

	input:
	file trimmed from trimmed_ch

	output:
	publishDir "${params.sam}", mode: 'copy'

	script:
	"""
	module load hisat2
	hisat2 \
                 -x ${params.ref}/GCF_000001405.39_GRCh38.p13_genomic.fna_index_hisat2 \
                 -1 ${params.trimmed}/${SAMPLE}_1.trimmed.fastq \
                 -2 ${params.trimmed}/${SAMPLE}_2.trimmed.fastq \
                 -S ${params.sam}/${SAMPLE}.sam \
                 -p 6 \
                --summary-file ${SAMPLE}.txt \
                --new-summary

	"""
}

process sam_to_bam {

        input:
        file sam from sam_ch

        output:
        publishDir "${params.sam}", mode: 'copy'

        script:
        """
        module load samtools
        samtools view -S -b ${SAMPLE}.sam | \
        samtools sort -n -o ${params.sam}/bam/${SAMPLE}.sorted.bam

        """
}

process generate_counts {

	input:
	file bam from bam_ch
	file ref from ref_ch

	output:
	publishDir "${params.outDir}", mode: 'copy'

	script:
	"""
	module load htseq
        htseq-count \
            -f bam \
            -r pos \
            -s no \
            -t exon \
            -i gene \
            ${params.bam}/${SAMPLE}.sorted.bam \
            ${params.ref}/GCF_000001405.39_GRCh38.p13_genomic.gtf \
            > ${params.bam}/${SAMPLE}.counts.txt

	"""
}

process import_reads {

        output:
        publishDir "${params.outDir}", mode: 'copy'

        script:
        """
	Rscript --vanilla ${params.rscript}
        """

}

process dge {

        output:
        publishDir "${params.outDir}", mode: 'copy'

        script:
        """
	Rscript --vanilla ${params.dgescript}
        """

}

process enrichment_up {

        output:
        publishDir "${params.outDir}", mode: 'copy'

        script:
        """
	Rscript --vanilla ${params.enrichmentupscript}
        """

}

process enrichment_down {

        output:
        publishDir "${params.outDir}", mode: 'copy'

        script:
        """
	Rscript --vanilla ${params.enrichmentdownscript}
        """

}

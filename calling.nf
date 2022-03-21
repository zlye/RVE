#!/usr/bin/env nextflow

// ######################################
// #    GATK4 
// #    snp and indell call from 
// #    GenomicsDB files
// ######################################


params.help=""
params.exe="slurm"

if (params.help) {
        log.info " "
        log.info "This is Nextflow pipeline for genotyping SNPS in rice 3k dataset. Most parameters are stored in file call_3k.conf"
        log.info " "
        log.info "Usage: nextflow run calling.nf -c calling.config --ref --chrom --(other options)"
        log.info "Options:"
        log.info "--help\t[BOOLEAN]\tShow this help message"
        log.info "--ref\t[STRING]\tPath to the indexed referece fasta file [OBLIGATORY]"
        log.info "--chrom\t[STRING]\tPath to the file with run IDs and fastq files to be processed [OBLIGATORY]"
//	log.info "--out\t[STRING]\tPath to directory where bamfiles will be written [OBLIG    ATORY]"
        log.info "--exe\t[STRING]\tExecutor mode, -local- or -slurm- [DEFUALT: local]"
        log.info " "
        exit 1
}

// Initalize Input
CHROM = file("${params.chrom}")
REF = file("${params.ref}")
OUT = file("${params.out}")
CONDA = <your conda>

// Create Input Channel
SampleData = Channel.fromPath("${CHROM}").splitCsv(header: ['CHR'], skip: 0, by:1)

// ########### Genotype GVCF ###############
process genotype {

	label 'himem'
	conda CONDA
	errorStrategy 'terminate'
	publishDir OUT, mode: 'move'	  
        executor = 'slurm'

        input:
        val(GVCF) from SampleData

        output:
        file({ "${CALLD}" }) into result

        script:
	CHR = "${GVCF.CHR}"
	CALLD = "${GVCF.CHR}.vcf"
	DB = "../../../../${GVCF.CHR}_db"

	"""
	gatk GenotypeGVCFs -V gendb://$DB -R ${REF} -L ${CHR} -O ${CALLD}
        """
}

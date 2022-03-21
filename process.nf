#!/usr/bin/env nextflow

// ######################################
// #    Zoe Lye 2020                    #
// #    GATK4 best practices            #
// #    Pipeline that performs:         #
// #    - mapping to reference          #
// #    - sorting reads		              #
// #    - removal PCR duplicates        #
// #    - calling haplotypes            #
// #	- storing g.vcfs in GenomicsDB	  #
// #					                          #
// ######################################

params.help=""

if (params.help) {
        log.info " "
        log.info "This is Nextflow pipeline for sno processing. Most parameters are stored in file process.conf"
        log.info " "
        log.info "Usage: nextflow run process.nf -c process.conf --ref --list --chrom --(other options)"
        log.info "Options:"
        log.info "--help\t[BOOLEAN]\tShow this help message"
        log.info "--ref\t[STRING]\tPath to the indexed referece fasta file [OBLIGATORY]"
        log.info "--list\t[STRING]\tPath to the file with run IDs and fastq files to be processed [OBLIGATORY]"
        log.info "--exe\t[STRING]\tExecutor mode, -local- or -slurm- [DEFUALT: local]"
        log.info " "
        exit 1
}

// Initalize Input
DAT = file("${params.list}")
REF = file("${params.ref}")
CHR = file("${params.chrom}")
CONDA = file("$params.conda")

// Create Input Channel
SampleData = Channel.fromPath("${DAT}").splitCsv(header: ['SID','RID','P1','P2'], skip: 0, by:1, sep:",")
SampleData2 = Channel.fromPath("${CHR}").splitCsv(header: ['CHR'], skip: 0, by:1)

// ########### Run BWA #############
process bwa {
	label 'bwa'
	conda CONDA
	errorStrategy 'terminate'

        input:
	val(RAW) from SampleData

        output:
        set val(ID), val(SID), file({ "${MAP}" }) into aligned

  script:
	      SID = "${RAW.SID}"
	      ID = "${RAW.RID}"
        RG = "\'@RG\\tID:${RAW.P1}\\tSM:${RAW.SID}\\tPL:ILLUMINA\'"
        MAP = "${RAW.RID}.sam"
        P1 = file("${RAW.P1}")
        P2 = file("${RAW.P2}")

        """
        bwa mem -t 6 -M -R ${RG} ${REF} ${P1} ${P2} > ${MAP}
	"""
}


// ########### Sort SAM ###############
process sort {
	label 'short_run'
	conda CONDA

        input:
        set val(ID), val(SID), file(SAM) from aligned

        output:
        set val(ID), val(SID), file({ "${SBAM}" }) into sorted

        script:

        SBAM = "${ID}.sorted.bam"
        del_sam = SAM.getName()
         
	"""
  picard SortSam INPUT=${SAM} OUTPUT=${SBAM} SORT_ORDER=coordinate
	"""
}

// ####### MERGE FILE CHANNELS #########

sorted
	.map {id, sid, file -> tuple(sid, file) }
	.groupTuple()
	.set { merg_sid }


// ###### Remove Read Duplicates #######
process rmdup {
	label 'long_run'
	conda '/home/znl207/mysoftware/miniconda3'

        input:
         set val(SID), file(SBAM) from merg_sid

        output:
         set val(SID), file({ "${RMDUP}" }) into rmdup

        script:

	 ALL_IN = SBAM.collect { "INPUT=$it" }.join(' ')
         RMDUP = "${SID}.dedup.bam"
         RMET = "${SID}.dedup.met"

       	 """
         picard MarkDuplicates ${ALL_IN} OUTPUT=${RMDUP} METRICS_FILE=${RMET}
         """
}


// ########### Check BAM ##########
process checkbam {
	label 'short_run'
	conda CONDA

	input:
	set val(SID), file(RMDUP) from rmdup

	output:
	set val(SID), file(RMDUP) into bchecked

	script:
	"""
	picard ValidateSamFile INPUT=${RMDUP} MODE=SUMMARY
	"""
}


// ########### Index BAM ##########
process index {
	label 'short_run'
	conda CONDA
	publishDir params.bam_outdir, mode: 'copy'

        input:
        set val(SID), file(RMDUP) from bchecked

        output:
        set val(SID), file({ "${RMDUP}" }), file({ "${BAI}" }) into index

        script:

        BAI = "${SID}.dedup.bai"
        
      	"""
        picard BuildBamIndex INPUT=${RMDUP} OUTPUT=${BAI}
        """
}


// ######## HaplotypeCaller ########

process haplocall {
	label 'medmem'
	conda CONDA
	errorStrategy = 'retry'
	maxRetries = 3

 	input:
 	set val(SID), file(RMDUP), file(BAI) from index
 
 	output:
 	file({ "${CALL}" }) into (called)
	file({ "${VCI}" }) into (vcfindex)
 
 	script:
 	
 	CALL = "${SID}.g.vcf"
	VCI = "${SID}.g.vcf.idx"
 	"""
 	gatk HaplotypeCaller -ERC GVCF -R ${REF} -I ${RMDUP} -O ${CALL}	
 	"""
}

// ########## Check GVCF ###########

process checkvcf {
	label 'long_run'
	conda '/home/znl207/mysoftware/miniconda3'

	input:
	file(CALL) from called
	
	output:
	file(CALL) into vchecked

	script:

	"""
	gatk ValidateVariants --java-options "-Xmx40G" \
	--validation-type-to-exclude ALLELES -R ${REF} -V ${CALL}
	"""
}

// ######## Convert into lists ##########

vchecked
	.toList()
	.set{ listed }

vcfindex
	.toList()
	.set{ vcflist }

// ######### CombineGVCFs ##########

process combine {
 	label 'himem'
	conda CONDA

	input:           
	val(CHROM) from SampleData2
	file(CALL) from listed
	file(VCI) from vcflist

	output:
	file( "${DB}" ) 

	script:
	CHR = "${CHROM.CHR}"
	COMB_IN = CALL.collect { "--variant $it" }.join(' ')
	DB = "../../../../${CHROM.CHR}_db"

	"""
	gatk GenomicsDBImport -R ${REF} ${COMB_IN} --genomicsdb-workspace-path $DB -L ${CHR}
	""" 
}

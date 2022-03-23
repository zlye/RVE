# RVE
Nextflow and python scripts acompanying the paper:
Deleterious mutations and the role of rare genetic variation on rice gene expression


Nextflow scripts - SNP and indel calling:
calling.config
calling.nf
process.conf
process.nf
 

Processing and merging SVs:

gridss_sample_processing.py
- input a vcf file of BND variants, the output of SV callers such as GRIDSS
- parses BND into duplications, deletion, insertions and copy number variable regions (complex rearrangement regions, where breakend intervals overlap within the sample)
- this script uses functions from 'write_vcf.py', 'cnvr_functions.py'

merge_samples.py
- merges SVs VCFs of individuals to group 
- this script uses functions from 'write_merged.pys'

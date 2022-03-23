import pysam
from itertools import count
from classes.Structural_variant import StructuralVariant

def write_to_vcf(sv_list, in_vcf, outfile):
    """
    write vcf of simplified SVs
    :param sv_list: if a list of sv objects one for each sv
    :param in_vcf: the original vcf file
    :param outfile: text name of output file
    uses pysam
    """
    # make new header from original vcf file
    samp_vcf = pysam.VariantFile(in_vcf)
    new_header = pysam.VariantHeader()
    # write contigs / chromosomes
    for x in samp_vcf.header.records:
        if (x.type == "CONTIG"):
            new_header.add_record(x)
    # get sample name
    sample_name = list(samp_vcf.header.samples)[0]
    split_string = sample_name.split("_", 1)
    newname = split_string[0]
    new_header.add_sample(newname)

    # add new meta data to header: filter, genotype, info
    new_header.add_meta('FILTER', items=[('ID', 'PASS'),
                                         ('Description', 'All filters passed')])
    new_header.add_meta('FORMAT', items=[('ID', "GT"),
                                         ('Number', 1), ('Type', 'String'),
                                         ('Description', 'Genotype')])
    new_header.add_meta('VAF', items=[('ID', "VAF"),
                                      ('Number', 1), ('Type', 'Integer'),
                                      ('Description', 'variant allele fraction')])
    new_header.add_meta('INFO', items=[('ID', 'SVTYPE'),
                                       ('Number', 1), ('Type', 'String'),
                                       ('Description', "Type of structural variant")])
    new_header.add_meta('INFO', items=[('ID', 'SVLEN'),
                                       ('Number', '.'), ('Type', 'Integer'),
                                       ('Description', "difference in length between REF and ALT alleles")])
    new_header.add_meta('INFO', items=[('ID', 'END'),
                                       ('Number', 1), ('Type', 'Integer'),
                                       ('Description', 'End position of the variant described in this record')])
    new_header.add_meta('INFO', items=[('ID', 'IMPRECISE'),
                                       ('Number', '0'), ('Type', 'Flag'),
                                       ('Description', 'Imprecise structural variation')])
    new_header.add_meta('INFO', items=[('ID', 'STRANDS'),
                                       ('Number', '.'), ('Type', 'String'),
                                       ('Description', 'Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)')])
    new_header.add_meta('INFO', items=[('ID', 'CIEND'),
                                       ('Number', 2), ('Type', 'Integer'), ('Description', 'Confidence interval around END for imprecise variants')])
    new_header.add_meta('INFO', items=[('ID', 'CIPOS'),
                                       ('Number', 2), ('Type', 'Integer'), ('Description', 'Confidence interval around remote breakend POS for imprecise variants')])
    new_header.add_meta('INFO', items=[('ID', 'GRIDSSEVENT'),
                                       ('Number', "."), ('Type', 'String'), ('Description', 'GRIDSS id for event associated to breakend')])
    new_header.add_meta('INFO', items=[('ID', 'INSCHROM'),
                                       ('Number', '.'), ('Type', 'String'), ('Description', 'Breakpoint pair chromosome of interchromosome breakpoint')])
    new_header.add_meta('INFO', items=[('ID', 'INSPOS'),
                                       ('Number', '.'), ('Type', 'Integer'), ('Description', 'insertion coordinate on interchromosome breakpoint')])
    new_header.add_meta('INFO', items=[('ID', 'CNVRCTX'),
                                       ('Number', '0'), ('Type', 'Flag'),
                                       ('Description', 'interchromosome breakpoint in CNVR')])
    new_header.add_meta('INFO', items=[('ID', 'CNVRINV'),
                                       ('Number', '0'), ('Type', 'Flag'),
                                       ('Description', 'inversion breakpoint in CNVR')])
    new_header.add_meta("INFO", items=[('ID', 'CNVRTYPE'),
                                       ('Number', '.'), ('Type', 'String'),
                                       ('Description', 'structural variant types associated with CNVR')])
    new_header.add_meta("INFO", items=[('ID', 'VAF'),
                                       ('Number', "."), ('Type', 'String'),
                                       ('Description', 'variant allele fraction of gridss breakends supporting this SV')])
    new_header.add_meta("INFO", items=[('ID', 'GTBNDS'),
                                       ('Number', "."), ('Type', 'String'),
                                       ('Description', 'genotype from all gridss breakend supporting this SV ')])
    sites_identifiers = count(1)
    # open file for writing
    newfile = pysam.VariantFile(outfile, 'w', header=new_header)
    for sv in sv_list:
        # fix formatting of the VAF and GT list
        pair_gt = ",".join([str(gt) for gt in sv.pair_gt])
        pair_vaf = ",".join([str(vaf) for vaf in sv.pair_vaf])
        rec = newfile.new_record(contig=sv.chromosome,
                                 start=(sv.start - 1),
                                 stop=sv.stop,
                                 alleles=sv.alleles,
                                 qual=sv.qual,
                                 id=str(next(sites_identifiers)),
                                 filter=sv.filter,
                                 info={'SVTYPE': sv.SVTYPE,
                                       'SVLEN': sv.SVLEN,
                                       'END': sv.stop,
                                       'IMPRECISE': sv.imprecise,
                                       'STRANDS': sv.strand,
                                       'CIEND': sv.ciend,
                                       'CIPOS': sv.cipos,
                                       'GRIDSSEVENT': sv.gridss_id,
                                       'INSCHROM': sv.inschrom,
                                       'INSPOS': sv.inspos,
                                       'CNVRINV': sv.cnvr_inv,
                                       'CNVRCTX': sv.cnvr_ctx,
                                       'CNVRTYPE': sv.cnvr_types,
                                       'GTBNDS': pair_gt,
                                       'VAF': pair_vaf})
        rec.samples[newname]['GT'] = sv.genotype
        newfile.write(rec)
    newfile.close()

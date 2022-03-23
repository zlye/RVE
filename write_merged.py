
from itertools import count
import operator
import datetime


def write_merged_file(output_vcf, chrom_dict,
                      sites_interval_trees, samples_dict):
    """
    write vcf of merged data
    :param output_vcf: Name of VCF output file
    :param chrom_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param sites_interval_tress:
     - A dictionary containing chromosomes as keys
    and IntervalTree instances as values. The IntervalTree instances contain dictionaries with collapsed sites as keys and a set of non-collapsed sites as values. Sites are tuples following:
    (chrom, pos, alt, end, sv_type, inschrom, inspos).
    :param samples_dict: Dict containing sites as keys and genotype lines of samples as values.
    :return: 0 (integer)
    """
    sites_identifiers = count(1)
    sites_dict = {}
    for chrom in sites_interval_trees:
        for interval in sites_interval_trees[chrom]:
            # add the intervals to the sites dict
            sites_dict.update(interval.data)

    # sort sites if not empty
    sorted_sites = []
    if sites_dict:
        sorted_sites = sorted(sites_dict.keys())
    # should be all the sites
    header_elems = "CHROM\tPOS\tID\tSVTYPE\tEND\tINSCHROM\tINSPOS\t"
    sample_names = sorted(list(samples_dict.keys()))
    header = header_elems + ("\t".join(sample_names) + "\n")
    with open(output_vcf, "w") as output:
        output.write(header)
        for collapsed_site in sorted_sites:
            chrom = str(collapsed_site[0])
            pos = str(collapsed_site[1])
            identifier = str(next(sites_identifiers))
            alt = collapsed_site[2]
            end = str(collapsed_site[3])
            sv_type = collapsed_site[4]
            if sv_type == "CTX":
                ins_chrom = str(collapsed_site[5])
                ins_pos = str(collapsed_site[6])
            else:
                ins_chrom = "."
                ins_pos = "."
            variant_line_elems = [chrom, pos,
                                  identifier, sv_type, end,
                                  ins_chrom, ins_pos]
            # extract sample elements
            original_sites = sites_dict[collapsed_site]
            for sample_name in sample_names:
                sample_field = 0
                for original_site in original_sites:
                    if original_site in samples_dict[sample_name]:
                        sample_field_temp = samples_dict[sample_name][original_site]
                        sample_field = sum(
                            [int(i) for i in sample_field_temp.split("/")])
                variant_line_elems.append(str(sample_field))
            # create new line for variant
            variant_line = "\t".join(variant_line_elems)
            output.write(variant_line)
            output.write("\n")
    return 0


def write_to_output_vcf(output_vcf, chrom_length_dict, sites_interval_trees, samples_dict):
    """
    Produce a VCF file containing all CNV sites and all genotypes of all samples

    :param output_vcf: Name of VCF output file
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param sites_interval_tress: A dictionary containing chromosomes as keys
    and IntervalTree instances as values. The IntervalTree instances contain
    dictionaries with collapsed sites as keys and a set of non-collapsed sites
    as values. Sites are tuples
    (chrom, pos, alt, end, sv_type, inschrom, inspos).
    :param samples_dict: Dict containing sites as keys and genotype lines of samples as values.
    :return: 0 (integer)
    """
    # create identifier generator for sites
    sites_identifiers = count(1)
    sites_dict = {}
    for chrom in sites_interval_trees:
        for interval in sites_interval_trees[chrom]:
            sites_dict.update(interval.data)
    # sort sites, if not empty
    sorted_sites = []
    if sites_dict:
        sorted_sites = sorted(
            sites_dict.keys(), key=operator.itemgetter(0, 1, 3, 4, 2, 5, 6))
    # produce VCF header
    sample_names = sorted(list(samples_dict.keys()))
    sample_names_header = "\t".join(sample_names)
    vcf_header_elems = ["##fileformat=VCFv4.2",
                        "##fileDate={}".format(
                            datetime.date.today().strftime("%Y%m%d")),
                        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                        "##INFO=<ID=END,Number=.,Type=Integer,Description=\"End position of the variant described in this region\">",
                        "##INFO=<ID=INSCHROM,Number=.,Type=String,Description=\"Chromosome on which insertion site of the dispersed duplication is located\">",
                        # "##INFO=<ID=INSPOS,Number=.,Type=Integer,Description=\"Position of insertion site of the dispersed duplication\">",
                        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
                        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
                        "##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around POS for imprecise variants\">",
                        "##INFO=<ID=CIEND95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around END for imprecise variants\">",
                        "##ALT=<ID=DEL,Description=\"Deletion\">",
                        "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                        "##FORMAT=<ID=SUP,Number=1,Type=Float,Description=\"Median number of supporting reads\">",
                        "##FORMAT=<ID=RP,Number=1,Type=Float,Description=\"Median number of supporting discordantly aligned read pairs\">",
                        "##FORMAT=<ID=SR,Number=1,Type=Float,Description=\"Median number of supporting split reads\">",
                        "##FORMAT=<ID=TOOL,Number=.,Type=String,Description=\"Supporting tools\">",
                        "##FORMAT=<ID=RQ,Number=1,Type=Float,Description=\"Probabilistic score of random forest model\">",
                        "##FORMAT=<ID=DHFC,Number=1,Type=Float,Description=\"duphold depth fold-change\">",
                        "##FORMAT=<ID=DHBFC,Number=1,Type=Float,Description=\"duphold depth fold-change compared to bins with matching GC\">",
                        "##FORMAT=<ID=DHFFC,Number=1,Type=Float,Description=\"duphold depth flank fold-change compared to 1000bp left and right of event\">"]
    for chrom in chrom_length_dict:
        contig_line = "##contig=<ID={},length={}>".format(
            chrom, chrom_length_dict[chrom])
        vcf_header_elems.append(contig_line)
    vcf_header_elems.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
        sample_names_header))
    vcf_header = "\n".join(vcf_header_elems)
    with open(output_vcf, "w") as output_vcf:
        # write header
        output_vcf.write(vcf_header)
        output_vcf.write("\n")
        # loop through sites, writing information on each line for each site
        ref = "N"
        qual = "-1"
        filtered = "PASS"
        for collapsed_site in sorted_sites:
            chrom = str(collapsed_site[0])
            pos = str(collapsed_site[1])
            identifier = str(next(sites_identifiers))
            alt = collapsed_site[2]
            print("alt site:", alt)
            sv_type = collapsed_site[4]
            # if svtyper_formatting:
            #     if sv_type == "DUP:TANDEM":
            #         sv_type = "DUP"
            type_info_field = "SVTYPE={}".format(sv_type)
            info_field_elems = [type_info_field, ]
            if sv_type != "INS":
                end = "END={}".format(str(collapsed_site[3]))
                info_field_elems.append(end)
            # add insertion site for dispersed duplications
            # if sv_type == "DUP:DISPERSED":
            #     ins_chrom = "INSCHROM={}".format(str(collapsed_site[5]))
            #     info_field_elems.append(ins_chrom)
            #     ins_pos = "INSPOS={}".format(str(collapsed_site[6]))
            #     info_field_elems.append(ins_pos)
            # add dummy confidence interval values
            cipos_elem = "CIPOS=-10,10"
            info_field_elems.append(cipos_elem)
            ciend_elem = "CIEND=-10,10"
            info_field_elems.append(ciend_elem)
            cipos_95_elem = "CIPOS95=-10,10"
            info_field_elems.append(cipos_95_elem)
            ciend_95_elem = "CIEND95=-10,10"
            info_field_elems.append(ciend_95_elem)
            info_field = ";".join(info_field_elems)
            format_field = "GT:SUP:RP:SR:TOOL:RQ:DHFC:DHBFC:DHFFC"
            variant_line_elems = [chrom, pos, identifier, ref, alt, qual,
                                  filtered, info_field, format_field]
            # extract sample elements
            original_sites = sites_dict[collapsed_site]
            for sample_name in sample_names:
                sample_field = 0
                for original_site in original_sites:
                    if original_site in samples_dict[sample_name]:
                        sample_field = samples_dict[sample_name][original_site]
                if sample_field == 0:
                    sample_field_elems = ["0/0", "0", "0",
                                          "0", ".", "0", "-1", "-1", "-1"]
                    sample_field = ":".join(sample_field_elems)
                variant_line_elems.append(sample_field)
            # create new line for variant
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
    return 0

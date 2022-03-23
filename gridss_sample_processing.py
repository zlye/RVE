"""
Process vcf output from GRIDSS to resolve structural variants
(1) Genotype
(2) Cluster breakpoints to find possible CNVR
(3) Classify structural variant to: dup, del, ins, CNVR, bnd
(2) write to simplified vcf
- includes options to change filter for size
- includes options to cluster variants into CNVR - copy number variable regions (overlapping breakpoints)
"""
import argparse
import os
import sys
import pysam
import write_vcf
import cnvr_functions
from itertools import count
from collections import defaultdict
from classes.breakend import Breakend
from classes.genomic_coordinate import GenomicCoordinate
from classes.Structural_variant import StructuralVariant


def classify_bnd(breakend, cluster_VAF, cluster_GT):
    """
    classify the breakend object to :  DEL, DUP, INS, INV, CTX (interchromosomal), BND (unpaired breakend)
    : param breakend: Breakend object
    : returns Instance of Structural_variant
    """
    info_fields = list(breakend.variant_record.info)
    length = (abs(breakend.coordinate.position -
                  breakend.mate_coordinate.position))
    if breakend.variant_record.info["IMPRECISE"]:
        imprecise = True
    else:
        imprecise = False
    # UNPAIRED BREAKEND
    if breakend.strand == '.':
        SV = StructuralVariant(
            chromosome=breakend.coordinate.chromosome,
            id=breakend.identifier,
            start=breakend.coordinate.position, stop=breakend.mate_coordinate.position + 1,
            genotype=(1, 0),
            gridss_id=breakend.variant_record.info["EVENT"],
            qual=breakend.variant_record.qual,
            strand=breakend.strand,
            SVTYPE='BND', SVLEN=length,
            pair_vaf=cluster_VAF, pair_gt=cluster_GT,
            imprecise=imprecise)
        if "CIPOS" in info_fields:
            SV.cipos = breakend.variant_record.info["CIPOS"]
        if "CIEND" in info_fields:
            SV.ciend = breakend.variant_record.info["CIEND"]
    # INTERCHROMOSOMAL
    elif breakend.coordinate.chromosome != breakend.mate_coordinate.chromosome:
        SV = StructuralVariant(
            chromosome=breakend.coordinate.chromosome,
            id=breakend.identifier,
            start=breakend.coordinate.position,
            stop=breakend.coordinate.position + 1,
            genotype=breakend.genotype,
            gridss_id=breakend.variant_record.info["EVENT"],
            qual=breakend.variant_record.qual,
            strand=breakend.strand,
            SVTYPE='CTX', SVLEN=0,
            pair_vaf=cluster_VAF, pair_gt=cluster_GT,
            imprecise=imprecise,
            inschrom=breakend.mate_coordinate.chromosome,
            inspos=breakend.mate_coordinate.position)
        if "CIPOS" in info_fields:
            SV.cipos = breakend.variant_record.info["CIPOS"]
        if "CIEND" in info_fields:
            SV.ciend = breakend.variant_record.info["CIEND"]
    # INVERSIONS
    elif (breakend.strand == "++") or (breakend.strand == "--"):
        SV = StructuralVariant(
            chromosome=breakend.coordinate.chromosome,
            id=breakend.identifier,
            start=breakend.coordinate.position,
            stop=breakend.mate_coordinate.position,
            genotype=breakend.genotype,
            gridss_id=breakend.variant_record.info["EVENT"],
            qual=breakend.variant_record.qual,
            strand=breakend.strand, SVTYPE="INV",
            SVLEN=(breakend.mate_coordinate.position -
                   breakend.coordinate.position),
            imprecise=imprecise,
            pair_vaf=cluster_VAF, pair_gt=cluster_GT)
        if "CIPOS" in info_fields:
            SV.cipos = breakend.variant_record.info["CIPOS"]
        if "CIEND" in info_fields:
            SV.ciend = breakend.variant_record.info["CIEND"]
    # INSERTION
    # elif length == 1 or length == 0:
    elif len(breakend.inserted_sequence) >= length:
        SV = StructuralVariant(
            chromosome=breakend.coordinate.chromosome,
            id=breakend.identifier,
            start=breakend.coordinate.position, stop=breakend.mate_coordinate.position,
            genotype=breakend.genotype,
            gridss_id=breakend.variant_record.info["EVENT"],
            qual=breakend.variant_record.qual,
            strand=breakend.strand, SVTYPE='INS',
            SVLEN=len(breakend.inserted_sequence),
            imprecise=imprecise,
            pair_vaf=cluster_VAF, pair_gt=cluster_GT,
            alleles=('N', breakend.inserted_sequence))
    # DELETIONS and DUPLICATION
    elif breakend.coordinate.position < breakend.mate_coordinate.position:
        length = breakend.mate_coordinate.position - breakend.coordinate.position
        if breakend.strand == "+-":
            svtype = "DEL"
        elif breakend.strand == "-+":
            svtype = "DUP"
        SV = StructuralVariant(
            chromosome=breakend.coordinate.chromosome,
            id=breakend.identifier,
            start=breakend.coordinate.position, stop=breakend.mate_coordinate.position,
            genotype=breakend.genotype, gridss_id=breakend.variant_record.info["EVENT"],
            qual=breakend.variant_record.qual, strand=breakend.strand,
            SVTYPE=svtype, SVLEN=length, imprecise=imprecise,
            pair_vaf=cluster_VAF, pair_gt=cluster_GT,)
        if "CIPOS" in info_fields:
            SV.cipos = breakend.variant_record.info["CIPOS"]
        if "CIEND" in info_fields:
            SV.ciend = breakend.variant_record.info["CIEND"]
    else:
        print("WTF found a breakend I couldn't classify:{}".
              format(breakend.variant_record.info["EVENT"]))
        return None
    return SV


def structural_variant_from_cluster(breakend_cluster):
    """
    make a structural variant from each cluster
    input: breakend_cluster, a list of breakends
    return: structural variant object
    """
    structural_variants = []
    mate_ids = []
    breakend_cluster = list(breakend_cluster)
    breakend_cluster.sort(
        key=lambda x: [x.coordinate.chromosome, x.coordinate.position])
    # get the list of genotypes and gt for each bnd in the cluster
    cluster_GT = [breakend.genotype for breakend in breakend_cluster]
    cluster_VAF = [breakend.vaf for breakend in breakend_cluster]
    if len(breakend_cluster) > 2:  # cluster size > 2
        clu_dict = {}
        for breakend in breakend_cluster:
            # make dict items are list of breakends, key is the gridss
            clu_dict.setdefault(
                breakend.variant_record.info["EVENT"], []).append(breakend)
        sub_structural_variants = []
        # for each id in the dict
        for gridss_id in clu_dict.keys():
            breakend_pair = clu_dict[gridss_id]
            breakend_pair.sort(key=lambda x:
                               [x.coordinate.chromosome, x.coordinate.position])
            cluster_GT = [breakend.genotype for breakend in breakend_pair]
            cluster_VAF = [breakend.vaf for breakend in breakend_pair]
            # get structrual variant for each pair within the cluster
            structural_variant = classify_bnd(
                breakend_pair[0], cluster_VAF, cluster_GT)
            sub_structural_variants.append(structural_variant)

        # add to the whole list
        # for sv in sub_structural_variants:
        structural_variants.extend(sub_structural_variants)
    elif len(breakend_cluster) == 2:
        for breakend in breakend_cluster:
            # don't process the same breakend pair twice
            if breakend.identifier in mate_ids:
                continue
            # store mate id
            mate_ids.append(breakend.mate_id)
        structural_variant = classify_bnd(
            breakend_cluster[0], cluster_VAF, cluster_GT)
        # both end of the CTX are recorded in vcf
        if structural_variant:
            if structural_variant.SVTYPE == "CTX":
                CTX2 = classify_bnd(
                    breakend_cluster[1], cluster_VAF, cluster_GT)
                structural_variants.append(CTX2)
            structural_variants.append(structural_variant)
    elif len(breakend_cluster) == 1 and breakend_cluster[0].strand == ".":
        # is a BND
        structural_variant = classify_bnd(
            breakend_cluster[0], cluster_VAF, cluster_GT)
        structural_variants.append(structural_variant)
        # this is the issue where one pair was a pass and the other wasn't
        # else if breakend_cluster...?
    elif len(breakend_cluster) == 1 and breakend_cluster[0].strand != ".":
        structural_variants.append
    else:
        print("warning: can't classify cluster:")
        print("breakends are: {}".format(
            [breakend.variant_record.info["EVENT"] for breakend in breakend_cluster]))
        return None
    return(structural_variants)


def cluster_breakends(breakends):
    """
    input: dictionary of breakends
    cluster breakends into dictionary of breakends
    return dictionary of sets of breakends
    * mostly unncessary only using one discovery method
    """
    cluster_id_generator = count(1)
    # make empty dictionary for cluster_ids
    breakend_id_dict = defaultdict(int)
    if not breakends:
        print("error, breakends dict missing")
    # make list of breakendss
    sorted_breakends = list(breakends.values())
    sorted_breakends.sort(
        key=lambda x: [x.coordinate.chromosome, x.coordinate.position])
    # cluster_id_dict: cluster_id is key,set of breakpoint ids
    cluster_id_dict = defaultdict(set)
    cluster_id = next(cluster_id_generator)
    breakend_id = sorted_breakends[0].identifier
    breakend_mate_id = sorted_breakends[0].mate_id
    # each breakend is assigned the same cluster id
    breakend_id_dict[breakend_id] = cluster_id
    breakend_id_dict[breakend_mate_id] = cluster_id
    # breakends and the mate id assigned to cluster dict
    cluster_id_dict[cluster_id].add(breakend_id)
    cluster_id_dict[cluster_id].add(breakend_mate_id)

    # loop through the rest of the breakends
    for i in range(1, len(sorted_breakends)):
        breakend = sorted_breakends[i]
        # check if id is already clustered:
        breakend_id = breakend.identifier
        # look in breakend_id_dict for the breakend_id, if not present it will return 0
        current_cluster_id = breakend_id_dict[breakend_id]
        if current_cluster_id == 0:  # breakend hasn't assined a cluster
            # make a new cluster for itself and its mate
            new_cluster_id = next(cluster_id_generator)
            breakend_id_dict[breakend_id] = new_cluster_id
            breakend_id_dict[breakend.mate_id] = new_cluster_id
            # update cluster in dict
            cluster_id_dict[new_cluster_id].add(breakend_id)
            cluster_id_dict[new_cluster_id].add(breakend.mate_id)
            # update cluster id (current in loop)
            current_cluster_id = new_cluster_id

    # convert the dictionary of sets of breakend coordinates
    # in a dictionary of all cluster ids as kews and a set of breakends as values.
    breakend_dict = defaultdict(set)
    for breakend_id in breakends:  # for each input
        breakend = breakends[breakend_id]
        cluster_id = breakend_id_dict[breakend_id]
        breakend_dict[cluster_id].add(breakend)

    # filter for clusters without LOW QUAL / LOW QUAL List comprehension
    hi_qual_clusters = {cluster_id: v for cluster_id, v in
                        breakend_dict.items()
                        if ("PASS" in [bkend.filter for bkend in breakend_dict[cluster_id]])}
    return hi_qual_clusters


def parse_breakends(in_vcf):
    """
    input: vcf file
    output: dict of breakends
    breakends dict: key id, value:breakend class object
    parse vcf, calc variant allele fraction (VAF) and infer genotype
    """
    vcf = pysam.VariantFile(in_vcf)
    breakends = {}
    # keep track of adjencencies & BNDs to find mates
    id_generator = count(1)
    # adjacencies_dict: key(chromA, posa, chromB, posB), value id
    adjacencies_dict = defaultdict(int)
    for record in vcf.fetch():
        # initialize inserted seq
        inserted_sequence = "."
        variant_type = record.info["SVTYPE"]
        if variant_type == "BND":
            # get strasnedness
            alt_field = record.alts[0]
            if "[" in alt_field:
                alt = alt_field.split("[")
                if alt[0]:
                    strand = "+-"  # 3' to 5'
                    alt_allel = alt[0]
                elif alt[2]:
                    strand = "--"  # 5' to 5'
                    alt_allel = alt[2]
            elif "]" in alt_field:
                alt = alt_field.split("]")
                if alt[0]:
                    strand = "++"  # 3' to 3'
                    alt_allel = alt[0]
                elif alt[2]:
                    strand = "-+"  # 5' to 3'
                    alt_allel = alt[2]
            # BND type
            elif "]" or "[" not in alt_field:
                strand = "."
                alt_allel = alt_field
            inserted_sequence = alt_allel
            #print("inserted_sequence: {}".format(inserted_sequence))
            # calculate the genotypes based on VAF
            for sample in record.samples:
                try:
                    VAF = record.info['VF'] / (record.info['VF'] +
                                               record.info['REF'] + record.info['REFPAIR'])
                except ZeroDivisionError:
                    VAF = 0
                if (VAF >= .75):
                    GT = (1, 1)
                elif (VAF > .25) and (VAF < .75):
                    GT = (1, 0)
                else:
                    # not the beset solution
                    GT = (0, 0)
            # extract genomic position of breakend and mate
            breakend_chrom = record.chrom
            breakend_pos = record.pos
            # for the edge case where BND has no mate:
            if strand == ".":
                mate_chrom = breakend_chrom
                mate_pos = breakend_pos
            else:
                mate_chrom = alt[1].split(":")[0]
                mate_pos = int(alt[1].split(":")[1])
            breakend_coordinate = GenomicCoordinate(
                breakend_chrom, breakend_pos)
            mate_coordinate = GenomicCoordinate(mate_chrom, mate_pos)
            new_id = next(id_generator)
            # create breakend  object
            breakend = Breakend(new_id, breakend_coordinate,
                                mate_coordinate, strand,
                                record, GT, VAF, filter=record.filter.items()[
                                    0][0],
                                inserted_sequence=inserted_sequence)
            # add adjacency to adjacency dict
            breakend_adjacency = (breakend_chrom, breakend_pos,
                                  mate_chrom, mate_pos)
            adjacencies_dict[breakend_adjacency] = new_id
            # check if mate is already in breakend dictionary:
            # -> if the mate adjacency isn't a key in the dict
            # default dict returns zero:
            mate_adjacency = (mate_chrom, mate_pos,
                              breakend_chrom, breakend_pos)
            # this will be zero if it isn't in there
            mate_id = adjacencies_dict[mate_adjacency]
            breakends[new_id] = breakend
            if mate_id != 0:  # the pair is in the dict
                # update mate ids
                breakend.mate_id = mate_id
                breakends[mate_id].mate_id = new_id
    vcf.close()
    return breakends


def parse_cl_args(in_args):
    """
    Parse command line arguments
    :param in_args: All command line arguments
    :return)none
    """
    description = "Parse VCF files, extracting Structural variants in  VCF format"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_vcf", type=str, required=True,
                        help="Path to VCF file")
    parser.add_argument("-o", "--output_name", type=str,
                        required=True,
                        help="Name of VCF output file")
    parser.add_argument("--large_variant", dest='large_variant',
                        action='store_true',
                        help="produce additional output file of variants >100,000bp, default False")
    parser.add_argument("-maxsize", "--max_variant_size", type=float,
                        default=100000, help="maximum size of resported variants, default: 100000")
    parser.add_argument("-cnvr", "--cluster_cnvr", action="store_true",
                        help="calculate CNVR regions as overlapping dup, del, ins, inv, default: FALSE")
    args = parser.parse_args(in_args)
    return args


def main():
    # testing mode
    testing_mode = False
    if testing_mode == True:
        cluster_cnvr = True
        input_vcf = "<YOUR INPUT VCF>"
        max_variant_size = 100000
        output_name = "<YOUR OUTPUT>"
        large_variant = False

    elif testing_mode == False:
        #args = parser.parse_args()
        args = parse_cl_args(sys.argv[1:])
        cluster_cnvr = args.cluster_cnvr
        input_vcf = args.input_vcf
        max_variant_size = args.max_variant_size
        output_name = args.output_name
        large_variant = args.large_variant
    # parse VCF into breakends
    breakends = parse_breakends(input_vcf)  # dict of breakend
    breakend_clusters = cluster_breakends(breakends)
    final_sv_list = []
    for breakend_cluster in breakend_clusters.values():
        structural_variants = structural_variant_from_cluster(breakend_cluster)
        if structural_variants:
            for variant in structural_variants:
                final_sv_list.append(variant)

    # filter for size
    maxsize = max_variant_size
    sv_large = [sv for sv in final_sv_list if sv.SVLEN > maxsize]
    sv_ctx = [sv for sv in final_sv_list if sv.SVTYPE == "CTX"]

    if cluster_cnvr is True:
        sv_small = [sv for sv in final_sv_list if (
            (sv.SVLEN < maxsize) and (sv.SVTYPE != "BND") and (sv.SVTYPE != "CTX") and (sv.SVTYPE != "INS"))]
        # merge cnvrs - excluding CTX and BND
        sv_merge = cnvr_functions.merge_cnvrs(sv_small)
        # BND and CTX are treated differently
        sv_bnd = [sv for sv in final_sv_list if sv.SVTYPE == "BND"]
        sv_out = sv_merge + sv_bnd
    else:
        sv_out = [sv for sv in final_sv_list if
                  (sv.SVLEN < maxsize and sv.SVTYPE != "CTX")]
    sv_out.sort(key=lambda x: [x.chromosome, x.start])

    if testing_mode is True:
        print("raw breakends:")
        for bk in breakends.values():
            print([bk.coordinate.chromosome, bk.coordinate.position,
                   bk.mate_coordinate.position])

        print("list output")
        print([[sv.chromosome, sv.start, sv.stop, sv.id, sv.SVTYPE, sv.alleles]
               for sv in sv_out])

    # write output:
    else:
        base_output = os.path.splitext(output_name)[0]
        if large_variant is True:
            write_vcf.write_to_vcf(sv_large, input_vcf,
                                   base_output + "large_var.vcf")
        if sv_ctx:
            write_vcf.write_to_vcf(sv_ctx, input_vcf, base_output + ".ctx.vcf")
        write_vcf.write_to_vcf(sv_out, input_vcf, base_output + ".vcf")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

"""
Merge VCF files of single sample VCFs containing SVs
one merged VCF file / one bed file
"""
import argparse
import pysam
import sys
from write_merged import write_merged_file
from write_merged import write_to_output_vcf
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from statistics import mean
from math import ceil

"""
merge cnv vcfs generated from the gridss post process script
merge is reciprocal 50% by default and within 1000 bp
inputs:
1. file that is a list of vcfs
2. reference.fasta.fai - the ref fasta index 
3. output vcf
4. max fraction of reciprocal overlap
output format is similar to bed pe
excludes CTX
"""

def get_overlap(a, b):
    """
    Compute the overlap between two discrete and closed intervals (borders are inclusive) a and b
    :param a: An interval written as a list. E.g. [5,15]
    :param b: An interval written as a list. E.g. [10,20]
    :return: The size of the overlap
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def obtain_sites_and_genotypes(input_fns):
    """
    get CNV sites from a list of VCF files
    :return: set containing call sides, sites are tuples
    (chrom, pos, alt, end, sv_type, inschrom, inspos).
    dict containing sites as keys and genotype lines of samples as values.
    """
    # obtain all sites
    sites_set = set()
    # obtail all genotypes of samples
    samples_dict = {}
    for input_vcf_fn in input_fns:
        with pysam.VariantFile(input_vcf_fn) as vcf:
            # add new sample to sample dict
            sample = vcf.header.samples[0]
            samples_dict[sample] = defaultdict(int)
            for record in vcf.fetch():
                chrom = record.chrom
                pos = record.pos
                alt = record.alts[0]
                end = record.stop
                sv_type = record.info["SVTYPE"]
                if sv_type == "CTX":
                    print("error file contains CTX variants, these are being skipped")
                    continue
                else:
                    inschrom = ""
                    inspos = -1
                site_tuple = (chrom, pos, alt, end, sv_type, inschrom, inspos)
                # get gt info
                if len(list(record.samples)) > 1:
                    print("error: file with multiple samples")
                for sample in record.samples:
                    gt = record.samples[sample]['GT']
                # skip variants with 0/0 genotype probably low quality:
                if gt == (0, 0):
                    continue
                GT = "/".join(str(a) for a in gt)
                genotype_line = ":".join([GT])
                samples_dict[sample][site_tuple] = genotype_line
                sites_set.add(site_tuple)

    return sites_set, samples_dict


def get_interval_tree(sites_set, chrom_dict, reciprocal):
    """
    get interval tree where sites considered to be the same CNV are collapses
    :return: dictionary containing chromosomes as keys and IntervalTree instances containing collapsed sites as values. IntervalTree data instrances contain dicts with collapsed sites as the key
    and a set of original sites as values.
    """
    sites_interval_trees = {}
    # search for more than x percent reciprocal overlap
    for site in sites_set:
        chrom = site[0]
        pos = site[1]
        alt = site[2]
        end = site[3]
        sv_type = site[4]
        inschrom = site[5]
        inspos = site[6]
        # check chromosome present in interval tree, if not add
        if chrom not in sites_interval_trees:
            sites_interval_trees[chrom] = IntervalTree()
        # check if SVs with reciprocal overlap
        sv_interval = [pos, end]
        sv_interval_len = end - pos + 1
        interval_start = sv_interval[0] - 1
        interval_end = sv_interval[1]
        #print("this is the site interval: {}".format(sv_interval))
        if interval_start == interval_end:
            overlapping_svs = sites_interval_trees[chrom].at(interval_start)
        else:
            # get calls that overlap with the midpoint or a slightly wider range
            # aim is to reduce the search range as much as possible to speed up queries
            # get midpoint
            if sv_interval_len % 2 == 0:
                midpoint_start = mean([interval_start, interval_end])
                midpoint_end = mean([interval_start, interval_end])
            else:
                midpoint_start = mean([interval_start, interval_end]) - 0.5
                midpoint_end = mean([interval_start, interval_end]) + 0.5
            reciprocal_slop = ceil(reciprocal * sv_interval_len)
            reciprocal_interval_start = min(interval_start + reciprocal_slop,
                                            interval_end - reciprocal_slop)
            reciprocal_interval_end = max(interval_start + reciprocal_slop,
                                          interval_end - reciprocal_slop)

            if midpoint_start >= reciprocal_interval_start:
                if midpoint_start == midpoint_end:
                    overlapping_svs = sites_interval_trees[chrom].at(
                        midpoint_start)
                else:
                    overlapping_svs = sites_interval_trees[chrom].overlap(midpoint_start,
                                                                          midpoint_end)
            else:
                overlapping_svs = sites_interval_trees[chrom].overlap(
                    reciprocal_interval_start, reciprocal_interval_end)
        # check if there are any cnvs of the same type with enough reciprocal overlap, saving only the best hit
        best_hit = None
        for old_sv_interval_instance in overlapping_svs:
            overlap_found = False
            old_sv = next(iter(old_sv_interval_instance.data))
            if sv_type != old_sv[4]:
                continue  # overlap is with different type of SV
            old_sv_interval = [old_sv[1], old_sv[3]]
        # IS THIS WHAT I REALLY WANT TO CLASSIFY INS?
            if sv_type == "INS":
                if old_sv[2] == alt:
                    overlap_found = True
            else:
               # check if there is enough reciprocal overlap and breakpoints are within 1 kbp
                old_sv_interval_len = old_sv_interval[1] - \
                    old_sv_interval[0] + 1
                overlap = get_overlap(old_sv_interval, sv_interval)
                if overlap > reciprocal * sv_interval_len\
                        and overlap > reciprocal * old_sv_interval_len \
                        and abs(sv_interval[0] - old_sv_interval[0]) <= 1000 \
                        and abs(sv_interval[1] - old_sv_interval[1]) <= 1000:
                    overlap_found = True
            if overlap_found:
                #print("Overlap found")
                if best_hit is None:
                    best_hit = old_sv_interval_instance
                else:
                    # compare overlap of this hit with that of the best hit and replace best hit if necessary
                    best_hit_data = next(iter(best_hit.data))
                    best_hit_interval = [best_hit_data[1], best_hit_data[3]]
                    best_hit_overlap = get_overlap(
                        best_hit_interval, sv_interval)
                    if overlap > best_hit_overlap:
                        best_hit = old_sv_interval_instance
                    else:
                        # support is less or tied, just keep the previously found best hit
                        continue
        if best_hit:
            # update entry in svs dictionary
            #print("the best hit is: {}".format(best_hit))
            best_hit_sv = next(iter(best_hit.data))
            new_sv_original_sites = best_hit.data[best_hit_sv]
            new_sv_site = list(best_hit_sv)
            # update coordinates, taking the union of the intervals
            if sv_type == "INS":
                new_sv_site[1] = min(best_hit_sv[1], pos)
                new_sv_site[3] = max(best_hit_sv[3], end)
            else:
                new_sv_site[1] = min(best_hit_sv[1], pos)
                new_sv_site[3] = max(best_hit_sv[3], end)
                # ensure end does not exceed chrom length
                new_sv_site[3] = min(new_sv_site[3],
                                     chrom_dict[str(chrom)])

            # change new site to tuple, # add new site to original sites
            new_sv_site = tuple(new_sv_site)
            new_sv_original_sites.add(site)
            #print("new sv original sites: {} ".format(new_sv_original_sites))
            new_site_dict = {new_sv_site: set()}  # create new sv interval
            # update the dictionary
            for original_site in new_sv_original_sites:
                new_site_dict[new_sv_site].add(original_site)
            # remove old sv interval from interval tree
            sites_interval_trees[chrom].remove(best_hit)
            #print("sites interval trees remove best hit".format(sites_interval_trees))
            # add new interval
            if new_sv_site[4] == "INS":
                interval_start = max(0, pos - 11)
                interval_end = min(end + 10, chrom_dict[str(chrom)])
                sites_interval_trees[chrom].addi(
                    interval_start, interval_end, new_site_dict)
            else:
                sites_interval_trees[chrom].addi(
                    new_sv_site[1] - 1, new_sv_site[3], new_site_dict)
        else:
            # ensure that end does not exceed chromosome, add to new dict
            end = min(end, chrom_dict[str(chrom)])
            new_site = (chrom, pos, alt, end, sv_type, inschrom, inspos)
            new_site_dict = {new_site: set()}
            new_site_dict[new_site].add(new_site)

            if new_site[4] == "INS":
                # add 10 bp upstream and downstream for overlap queries
                interval_start = max(0, pos - 11)
                interval_end = min(end + 10, chrom_dict[str(chrom)])
                sites_interval_trees[chrom].addi(
                    interval_start, interval_end, new_site_dict)
            else:
                sites_interval_trees[chrom].addi(
                    new_site[1] - 1, new_site[3], new_site_dict)
    return sites_interval_trees


def merge_vcfs(chrom_dict, input_fn, output_vcf, reciprocal, vcf_format):
    """
    Merge sample VCF files of samples processed from BND 
    """
    input_fns = []
    with open(input_fn) as input_file:
        for line in input_file:
            input_fns.append(line.strip())

    # get sites and gt info
    sites_set, samples_dict = obtain_sites_and_genotypes(input_fns)
    # get interval tree in which sites considered to be the smae CNV are collapsed
    sites_interval_tree = get_interval_tree(sites_set, chrom_dict,
                                            reciprocal)

    if vcf_format == True:
        write_to_output_vcf(output_vcf, chrom_dict,
                            sites_interval_tree, samples_dict)
    if vcf_format == False:
        write_merged_file(output_vcf, chrom_dict,
                          sites_interval_tree, samples_dict)


def parse_fai(fai):
    """
    parse the fai file in to a dic
    """
    chrom_dict = {}
    with open(fai) as fai_file:
        for line in fai_file:
            line_elems = line.strip().split()
            chrom = line_elems[0]
            chrom_length = int(line_elems[1])
            chrom_dict[chrom] = chrom_length
    return(chrom_dict)


def parse_cl_args(in_args):
    """
    parse command line arguments
    """
    description = "Merge VCF files of single samples containing SVs, producing table format"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_file", type=str,
                        help="File containing path to single sample VCF file "
                             "on each line")
    parser.add_argument("-f", "--fai_fn", type=str,
                        help="Path to fasta index file of used genome")
    parser.add_argument("-o", "--output", type=str,
                        help="Name of VCF output file")
    parser.add_argument("-r", "--reciprocal", type=float, default=0.5,
                        help="Minimum fraction of reciprocal overlap needed "
                             "to collapse calls")
    parser.add_argument("-v", "--out_format", help="write output as a vcf",
                        action='store_true')
    args = parser.parse_args(in_args)
    return args


def main():
    args = parse_cl_args(sys.argv[1:])
    input_file = args.input_file
    output_file = args.output
    reciprocal = args.reciprocal
    out_format = args.out_format
    fai_fn = args.fai_fn

    chrom_dict = parse_fai(fai_fn)
    if reciprocal < 0 or reciprocal > 1:
        raise ValueError("Reciprocal overlap must be between 0 and 1")
    merge_vcfs(chrom_dict, input_file, output_file,
               reciprocal, out_format)


if __name__ == "__main__":
    main()

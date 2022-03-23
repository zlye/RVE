"""
Functions for dealing with putative CNVRs while post processing raw gridss output

"""

from classes.Structural_variant import StructuralVariant
def make_cnvr_sv(group):
    """
    input: list of structural variant object (group)
    output: structural variant object
    combine a list of sv objects into a singple sv object representing a CNVR
    """
    sv_coord = [sv.start for sv in group] + [sv.stop for sv in group]
    svlen = max(sv_coord) - min(sv_coord)
    gridss_id = ",".join([sv.gridss_id for sv in group])
    qual = min([sv.qual for sv in group])  # take min qual
    # list comprehension to un nest
    genotypes = [sv.pair_gt for sv in group]
    genotypes = [gt for sample in genotypes for gt in sample]
    # to get genotype: take most freq:
    gt_counts = {x: genotypes.count(x) for x in genotypes}
    v = list(gt_counts.values())
    k = list(gt_counts.keys())
    genotype = k[v.index(max(v))]
    strands = [sv.strand for sv in group]
    inv = True if ("--" or "++") in strands else False
    # get the sv types - alll bnd not cnvr
    svtypes = [sv.SVTYPE for sv in group]
    imprescise = True
    if all(sv == "BND" for sv in svtypes):
        svtype = "BND"
        imprescise = False
    else:
        svtype = "CNVR"
    cnvr_ctx = True if "CTX" in svtypes else False
    # get vaf_list gt_list:
    gt_all = [sv.pair_gt for sv in group]
    vaf_list = [sv.pair_vaf for sv in group]
    # make CNVR object
    cnvr = StructuralVariant(
        chromosome=group[0].chromosome,
        id=group[0].id,
        start=min(sv_coord), stop=max(sv_coord),
        genotype=genotype, gridss_id=gridss_id,
        qual=qual, strand="+-",
        SVTYPE=svtype, SVLEN=svlen, imprecise=imprescise,
        pair_vaf=vaf_list, pair_gt=gt_all,
        cnvr_inv=inv, cnvr_ctx=cnvr_ctx, cnvr_types=",".join(svtypes))
    return(cnvr)

def merge_cnvrs(structural_variants):
    """
    input list of SVs
    output modified list of SV with overlapping SVs as CNVRs
    """
    structural_variants.sort(key=lambda x: [x.chromosome, x.start])
    merge_sv = []
    # intialize group and prior range
    first_sv = structural_variants[0]
    prior_range = [first_sv.chromosome, first_sv.start, first_sv.stop]
    #prior_range = [first_sv.chromosome, first_sv.start - 20, first_sv.stop + 20]
    prior_sv = first_sv
    group = [first_sv]  # make a list
    # test overlap of ranges for consecutive structural variants
    for sv in (structural_variants[1:]):
        #test_range = [sv.chromosome, sv.start - 20, sv.stop + 20]
        test_range = [sv.chromosome, sv.start, sv.stop]
        if sv.SVTYPE == "BND":
            test_range[1] = sv.start - 25
        if sv.SVTYPE == "CTX":
            merge_sv.append(sv)
            test_range[2] = sv.stop + 25

        # IF LAST SV
        if sv == structural_variants[-1]:
            # if overlap make cnvr, append
            if test_range[0] == prior_range[0] and (prior_range[1] <= test_range[1] <= prior_range[2] or prior_range[1] <= test_range[2] <= prior_range[2]):
                group.append(sv)
                cnvr = make_cnvr_sv(group)
                merge_sv.append(cnvr)
            # No overlap, append prior and current
            else:
                if len(group) == 1:
                    merge_sv.append(prior_sv)
                    merge_sv.append(sv)
            break
        # test if overlap == True
        if test_range[0] == prior_range[0] and (prior_range[1] <= test_range[1] <= prior_range[2] or prior_range[1] <= test_range[2] <= prior_range[2]):
            if sv in group:
                # overlap but sv already in group, make CNVR SV
                cnvr = make_cnvr_sv(group)
                merge_sv.append(cnvr)
            else:
                # SV not in group, append SV to group
                group.append(sv)
                prior_range[2] = max([prior_range[2], test_range[2]])
                prior_range[1] = min([prior_range[1], test_range[1]])
        else:  # No overlap:
            if len(group) == 1:
                # no overlap sv to list
                prior_range = test_range
                merge_sv.append(prior_sv)
                prior_sv = sv
                group.clear()
                group = [sv]
            if len(group) > 1:
                # no overlap, make new cnvr SV
                cnvr = make_cnvr_sv(group)
                merge_sv.append(cnvr)
                group.clear()
                group = [sv]
                prior_range = test_range
                prior_sv = sv
    return(merge_sv)

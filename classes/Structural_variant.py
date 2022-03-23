"""
class so represent the SVs that will be written to VCF
will be written to a new record
"""
class StructuralVariant(object):
    def __init__(self,
                 chromosome: str,
                 id: str,
                 start: int,
                 stop: int,
                 genotype: tuple,
                 gridss_id: str,
                 qual: float,
                 strand: str,
                 SVTYPE: str,
                 SVLEN: int = 0,
                 cipos: tuple = (None, None),
                 ciend: str = (None, None),
                 alleles: tuple = ('N', 'N'),
                 imprecise: bool = False,
                 filter: str = "PASS",
                 pair_vaf: str = '',
                 pair_gt: str = '',
                 inschrom: str = '',
                 inspos: int = None,
                 cnvr_inv: bool = False,
                 cnvr_ctx: bool = False,
                 cnvr_types: str = None):
        self.chromosome = chromosome
        self.id = id
        self.start = start
        self.stop = stop
        self.alleles = alleles
        self.genotype = genotype
        self.gridss_id = gridss_id
        self.cipos = cipos
        self.ciend = ciend
        self.qual = qual
        self.strand = strand
        self.SVTYPE = SVTYPE
        self.SVLEN = SVLEN
        self.alleles = alleles
        self.imprecise = imprecise
        self.filter = filter
        self.pair_vaf = pair_vaf
        self.pair_gt = pair_gt
        self.inschrom = inschrom
        self.inspos = inspos
        self.cnvr_inv = cnvr_inv
        self.cnvr_ctx = cnvr_ctx
        self.cnvr_types = cnvr_types

class CTXVariant(object):
    def __init__(self,
                 chrom: str,
                 pos: int,
                 inschrom: str,
                 inspos: int,
                 genotype: tuple,
                 sample_id: str,
                 gridss_id: str):
        self.chrom = chrom
        self.pos = pos
        self.inschrom = inschrom
        self.inspos = inspos
        self.sample_id = sample_id
        self.gridss_id = gridss_id

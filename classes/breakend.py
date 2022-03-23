"""
Class to represent breakends of VCF files
"""
from classes.genomic_coordinate import GenomicCoordinate
from pysam import VariantRecord

class Breakend(object):
    def __init__(self, identifier: int,
                 coordinate: GenomicCoordinate,
                 mate_coordinate: GenomicCoordinate,
                 strand: str,
                 variant_record: VariantRecord,
                 genotype: tuple,
                 vaf: int,
                 filter: str,
                 mate_id: int = None,
                 inserted_sequence: str = "."
                 ):
        self.identifier = identifier
        self.coordinate = coordinate
        self.mate_coordinate = mate_coordinate
        self.strand = strand
        self.variant_record = variant_record
        self.genotype = genotype
        self.vaf = vaf
        self.filter = filter
        self.inserted_sequence = inserted_sequence
        self.mate_id = mate_id

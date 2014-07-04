"""Defines the Transcript object, which represents a row in a UCSC table.

"""
import itertools

from probe_generator import probe
from probe_generator.sequence import SequenceRange

_REQUIRED_FIELDS = (
    # Fields which are assumed to exist by Transcript methods.
    'name',
    'exonStarts',
    'exonEnds',
    'cdsStart',
    'cdsEnd',
    'chrom',
    )

_GENE_NAME_FIELDS = (
    # The names of the fields which might contain the name of a gene in any of
    # the supported UCSC file formats.
    'name2',
    'proteinID',
    )

class Transcript(object):
    """Represents a UCSC annotation of a transcript.

    Public attributes:

        name:        unique transcript identifier
        gene_id:     non-unique gene name identifier
        chromosome:  self-evident
        plus_strand: boolean true if the transcript is on the plus strand

    """
    def __init__(self, spec):
        """`spec` is a dict containing the information from a row read from a
        UCSC annotation table.

        Raises an InvalidAnnotationFile error when the spec does not have the
        required fields.

        """
        self._spec = spec
        self._assert_spec_correct()

        self.name = self._spec['name']
        self.chromosome = self._spec['chrom'].lstrip('chr')
        self.gene_id, = (self._spec[field] for field in _GENE_NAME_FIELDS
                         if field in self._spec)
        self.plus_strand = self._spec['strand'] == '+'

    def _assert_spec_correct(self):
        """Raises and InvalidAnnotationFile exception unless all of the
        _REQUIRED_FIELDS and exactly one of the _GENE_NAME_FIELDS are present
        in the _spec.

        """
        if not all(field in self._spec for field in _REQUIRED_FIELDS):
            raise InvalidAnnotationFile(
                "Annotation file is missing required fields: {}".format(
                    [field for field in _REQUIRED_FIELDS
                     if not field in self._spec]))
        gene_names = [field for field in _GENE_NAME_FIELDS
                      if field in self._spec]
        if not len(gene_names) == 1:
            raise InvalidAnnotationFile(
                "Anotation file contains gene id fields: {}. "
                "Expected exactly one of {}".format(
                    gene_names, _GENE_NAME_FIELDS))

    def exons(self):
        """Return the exon positions of a UCSC annotation feature.

        In a UCSC annotation file, the positions of the starts and ends of exons
        are stored as comma-separated strings:

            '20,30,40,'

        Given a dictionary with this data, we return of a list of tuples:

            (exonStart, exonEnd)

        If the 'strand' of the row is '-', the function return the exons in
        reversed order. In this case, the first exon relative the the direction of
        transcription (which is probably what the user means, is the last exon
        along the chromosome reading from left to right along the '+' strand (which
        is how the data are stored in UCSC tables).

        Raises a FormattingError when the `row` does not appear to come from a
        valid UCSC gene table.

        """
        exon_starts = self._spec['exonStarts'].split(',')
        exon_ends = self._spec['exonEnds'].split(',')
        positions = []
        for start, end in zip(exon_starts, exon_ends):
            if start != '' and end != '':
                start, end = int(start), int(end)
                positions.append((start, end))
        if not self.plus_strand:
            positions.reverse()
        return [SequenceRange(self.chromosome, start, end)
                for start, end in positions]

    def coding_exons(self):
        """As in `exons`, but with the UTRs trimmed out.

        """
        cds_start = int(self._spec['cdsStart'])
        cds_end = int(self._spec['cdsEnd'])
        exon_positions = self.exons()
        positions = []

        if not self.plus_strand:
            exon_positions.reverse()
        for exon in exon_positions:
            if exon.end < cds_start:
                pass
            elif exon.start <= cds_start <= exon.end:
                positions.append((cds_start, exon.end))
            elif exon.start <= cds_end <= exon.end:
                positions.append((exon.start, cds_end))
                break
            else:
                positions.append((exon.start, exon.end))
        if not self.plus_strand:
            positions.reverse()
        return [SequenceRange(self.chromosome, start, end)
                for start, end in positions]

    def exon(self, index):
        """Return the exon at the index of the transcript, raising a NoFeatureError
        when the index is out of range.

        """
        try:
            return self.exons()[index-1]
        except IndexError as error:
            raise NoFeature(error)

    def nucleotide_index(self, index):
        """Given a base pair index and a row of a UCSC gene table, return the
        genomic coordinate of the base pair at that index in the transcript.

        'transcript' is a row from a UCSC genome annotation table.

        """
        base_index = self._transcript_index(index)
        return SequenceRange(self.chromosome, base_index, base_index+1)

    def codon_index(self, index):
        """Given a codon index and a row of a UCSC gene table, return the genomic
        coordinate of the second base pair of that codon.

        """
        base_index = self._transcript_index((index * 3) - 2)
        return SequenceRange(self.chromosome, base_index, base_index+3)

    def _transcript_index(self, index):
        """Return the genomic coordinate of the base pair at the index as an
        integer.

        """
        indices = []
        for exon in self.coding_exons():
            nucleotide_range = list(range(exon.start, exon.end))
            if not self.plus_strand:
                nucleotide_range.reverse()
            indices.append(nucleotide_range)
        base_coordinates = itertools.chain(*indices)
        try:
            base_index = next(itertools.islice(
                    base_coordinates,
                    index,
                    None))
        except StopIteration:
            raise OutOfRange(
                "Base {} is outside the range of transcript '{}'".format(
                    index, self.name))
        return base_index


class InvalidAnnotationFile(Exception):
    """Raised when format assumptions about the table used to generate the
    transcript annotations are violated.

    """


class OutOfRange(probe.NonFatalError):
    """Raised when a base index outside the range of a transcript is specified.

    """


class NoFeature(probe.NonFatalError):
    """Raised when the index of the exon is outside the range of the
    transcript.

    """

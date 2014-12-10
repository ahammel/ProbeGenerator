"""Defines the Transcript object, which represents a row in a UCSC table.

"""
import itertools

from probe_generator import probe
from probe_generator.sequence_range import SequenceRange

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

    def __hash__(self):
        return hash(tuple([value for value in sorted(self._spec.values())]))

    def __len__(self):
        """Return the number of coding nucleotides in the transcript.

        """
        return sum(exon.end - exon.start for exon in self.coding_exons())

    def _assert_spec_correct(self):
        """Raises an InvalidAnnotationFile exception unless all of the
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
                "Annotation file contains gene id fields: {}. "
                "Expected exactly one of {}".format(
                    gene_names, _GENE_NAME_FIELDS))

    def exons(self):
        """Return the exon positions of a UCSC annotation feature.

        In a UCSC annotation file, the positions of the starts and ends of exons
        are stored as comma-separated strings:

            '20,30,40,'

        Given a dictionary with this data, we return a list of tuples:

            (exonStart, exonEnd)

        If the 'strand' of the row is '-', the function return the exons in
        reversed order. In this case, the first exon relative to the direction
        of transcription (which is probably what the user means), is the last
        exon along the chromosome reading from left to right along the '+'
        strand (which is how the data are stored in UCSC tables).

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
            elif exon.start <= cds_start <= cds_end <= exon.end:
                positions.append((cds_start, cds_end))
                break
            elif exon.start <= cds_start <= exon.end:
                positions.append((cds_start, exon.end))
            elif cds_start <= exon.start <= exon.end <= cds_end:
                positions.append((exon.start, exon.end))
            elif exon.start <= cds_end <= exon.end:
                positions.append((exon.start, cds_end))
                break
            elif cds_end <= exon.start:
                break
            else:
                assert False, "unreachable: {}/{}".format(self.name, self.gene_id)
        if not self.plus_strand:
            positions.reverse()
        return [SequenceRange(self.chromosome, start, end)
                for start, end in positions]

    def exon(self, index):
        """Given the one-based index of an exon, return a SequenceRange object
        representing the genomic coordinates of that exon.

        Raise a NoFeature error when the exon is out of the bounds of the
        transcript.

        """
        try:
            return self.exons()[index-1]
        except IndexError as error:
            raise NoFeature("{}: {}/{}".format(
                    error, self.gene_id, self.name))

    def nucleotide_index(self, index):
        """Given a 1-based base pair index, return a SequenceRange object
        representing the base pair at that index in the transcript.

        """
        base_index = self._transcript_index(index)
        return SequenceRange(self.chromosome, base_index, base_index+1)

    def codon_index(self, index, reference=None, mutation=None):
        """Given a 1-based codon index, return a SequenceRange object
        representing that codon.

        The `reference` and `mutation` parameters are passed to the
        SequenceRange object.

        """
        base_index = self._transcript_index(index*3)
        if self.plus_strand:
            return SequenceRange(
                    self.chromosome,
                    base_index-2,
                    base_index+1,
                    reference=reference,
                    mutation=mutation)
        else:
            return SequenceRange(
                    self.chromosome,
                    base_index,
                    base_index+3,
                    reference=reference,
                    mutation=mutation)

    def base_index(self, sequence_range):
        """Given a SequenceRange object representing a genomic location within
        the transcript, return the one-based index of nucleotide at the start
        of the sequence-range object.

        Raises an OutOfRange error when the `sequence_range` is not within the
        transcript.

        """
        for i in range(1, len(self)+1):
            nucleotide = self.nucleotide_index(i)
            if sequence_range.start == nucleotide.start:
                return i
        raise OutOfRange

    def transcript_range(self, start, end):
        """Return a list of SequenceRange objects representing the genomic
        location(s) of the transcript from `start` to `end`.

        More than one SequenceRange is returned if the requested range crosses
        exon boundaries.

        The `start` and `end` variables are 1-based left-inclusive,
        right-exclusive.

        """
        ranges = [self.nucleotide_index(i) for i in range(start, end)]
        return SequenceRange.condense(*ranges)

    def _transcript_index(self, index):
        """Given the 1-based index of a nucleotide in the coding sequence,
        return the 0-based genomic index of that nucleotide as an integer.

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
                    index-1,
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

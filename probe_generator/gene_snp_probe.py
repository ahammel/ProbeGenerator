"""Probe for an SNP mutation based on the index of a nucleotide in a coding
sequence.

"""
import re
import itertools

from probe_generator import annotation, probe
from probe_generator.snp_probe import SnpProbe

_STATEMENT_REGEX = re.compile("""
        \s*                # whitespace
        ([a-zA-Z0-9_./-]+) # gene name
        \s*
        :
        \s*
        c\.
        ([0-9]+)           # base number
        \s*
        ([ACGTacgt])       # reference base
        \s*
        >
        \s*
        ([ACGTacgt])       # mutation base
        \s*
        /
        \s*
        ([0-9]+)           # number of base pairs
        \s*
        (--.*|\s*)         # comment
        """, re.VERBOSE)

_STATEMENT_SKELETON = ("{gene}:c.{base}{reference}>{mutation}/{bases}_"
                       "{transcript}_{chromosome}:{index}{comment}")


class GeneSnpProbe(object):
    """Probe for a single nucleotide mutation event at a base pair specified
    relative to the start of a transcript.

    """
    # TODO: This class is somewhat tightly-coupled to the SnpProbe class. If I
    # need to change to format of the spec dictionary for some reason, I will
    # probably have to change both classes. If this starts to become a problem
    # we can use an intermediate ProbeSpecification class.

    def __init__(self, specification):
        self._spec = specification

    def __str__(self):
        return _STATEMENT_SKELETON.format(**self._spec)

    def sequence(self, genome):
        """Return the probe sequence given a genome annotation.

        """
        snp_probe = SnpProbe(self._spec)
        return snp_probe.sequence(genome)

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given a gene SNP probe statement, return all the probes which match
        the specification.

        If more than one probe has identical genomic coordinates, only the
        first is returned.

        """
        if genome_annotation is None:
            genome_annotation = []
        partial_spec = _parse(statement)
        transcripts = annotation.lookup_gene(
            partial_spec["gene"], genome_annotation)
        for transcript in transcripts:
            try:
                index = _relative_index(partial_spec["base"], transcript)
            except OutOfRange as error:
                raise OutOfRange(
                    "{} in statement:\n'{}'".format(error, statement))
            spec = dict(partial_spec,
                        chromosome=transcript["chrom"],
                        transcript=transcript["name"],
                        index=index
                        )
            yield GeneSnpProbe(spec)

def _parse(statement):
    """Return a partial GeneSnpProbe specification given a probe statement.

    Raises an InvalidStatement excpetion when fed an invalid gene snp
    statement.

    """
    match = _STATEMENT_REGEX.match(statement)

    if not match:
        raise probe.InvalidStatement
    (gene,
     base,
     reference,
     mutation,
     bases,
     comment) = match.groups()
    return {"gene": gene,
            "base": int(base),
            "reference": reference,
            "mutation": mutation,
            "bases": int(bases),
            "comment": comment}


def _relative_index(index, transcript):
    """Given a base pair index and a row of a UCSC gene table, return the
    genomic coordinate of the base pair at that index in the transcript.

    """
    strand = transcript["strand"]
    indices = (_base_indices(pair, strand)
               for pair in annotation.exons(transcript))
    base_coorinates = itertools.chain(*indices)
    try:
        base_index = next(itertools.islice(
                base_coorinates,
                index-1, # Convert between zero- and one-indexing
                index))
    except StopIteration:
        raise OutOfRange(
            "\nBase {} is outside the range of transcript '{}'".format(
                index, transcript["name"]))
    return base_index


def _base_indices(exon_range, strand):
    """Given an exon range (int, int) the strand of the exon ("+" or "-"),
    return a generator of the genomic coordinates of the bases in the exon.

    """
    assert strand in "+-"

    p, q = exon_range
    if strand == "+":
        return range(p, q+1)
    else:
        return reversed(range(p, q+1))


class OutOfRange(probe.NonFatalError):
    """Raised when a base index outside the range of a transcript is specified.

    """

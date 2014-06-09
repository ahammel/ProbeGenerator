"""Probe for an SNP mutation based on the index of a nucleotide in a coding
sequence.

"""
import re
import sys

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
        """Return the probe sequence given a genome sequence.

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
        cached_coordinates = set()
        for transcript in transcripts:
            try:
                index = annotation.nucleotide_index(
                    partial_spec["base"], transcript)
            except annotation.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                chromosome = transcript["chrom"].lstrip("chr")
                transcript = transcript["name"]
                spec = dict(partial_spec,
                            chromosome=chromosome,
                            transcript=transcript,
                            index=index)
                if not (chromosome, index) in cached_coordinates:
                    cached_coordinates.add((chromosome, index))
                    yield GeneSnpProbe(spec)


def _parse(statement):
    """Return a partial GeneSnpProbe specification given a probe statement.

    Raises an InvalidStatement exception when fed an invalid gene snp
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

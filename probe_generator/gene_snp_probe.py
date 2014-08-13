"""Probe for an SNP mutation based on the index of a nucleotide in a coding
sequence.

"""
import re
import sys

from probe_generator import annotation, transcript
from probe_generator.variant import TranscriptVariant, GenomeVariant
from probe_generator.probe import AbstractProbe, InvalidStatement

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
        (\[trans\]|)       # transcript-only sequence
        \s*
        /
        \s*
        ([0-9]+)           # number of base pairs
        \s*
        (--.*|\s*)         # comment
        """, re.VERBOSE)


class GeneSnpProbe(AbstractProbe):
    """Probe for a single nucleotide mutation event at a base pair specified
    relative to the start of a transcript.

    """
    _STATEMENT_SKELETON = ("{gene}:c.{base}{reference}>{mutation}"
                           "{transcript_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{index_base}{comment}")

    def __init__(self, *, variant, index, comment):
        self.variant = variant
        self.index = index
        self.comment = comment

    def __str__(self):
        return self._STATEMENT_SKELETON.format(
            gene=self.variant.gene,
            base=self.index,
            reference=self.variant.reference,
            mutation=self.variant.mutation,
            transcript_sequence='[trans]' if self.variant.is_transcript else '',
            bases=self.variant.bases,
            transcript_name=self.variant.transcript_name,
            chromosome=self.variant.chromosome,
            index_base=self.variant.coordinate,
            comment=self.comment)

    def get_ranges(self):
        return self.variant.sequence_ranges()

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given a gene SNP probe statement, return all the probes which match
        the specification.

        If more than one probe has identical genomic coordinates, only the
        first is returned.

        """
        probes = []

        if genome_annotation is None:
            genome_annotation = []

        specification = _parse(statement)

        if specification["transcript_sequence"] == '':
            variant_class = GenomeVariant
        else:
            variant_class = TranscriptVariant

        transcripts = annotation.lookup_gene(
            specification["gene"], genome_annotation)

        cached_coordinates = set()
        for txt in transcripts:
            base = specification["base"]
            try:
                index = txt.nucleotide_index(base)
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                if not (index) in cached_coordinates:
                    cached_coordinates.add(index)
                    variant = variant_class(
                        transcript=txt,
                        reference=specification["reference"],
                        mutation=specification["mutation"],
                        index=index,
                        bases=specification["bases"])
                    probes.append(GeneSnpProbe(
                            variant=variant,
                            index=base,
                            comment=specification["comment"]))
        return probes


def _parse(statement):
    """Return a partial GeneSnpProbe specification given a probe statement.

    Raises an InvalidStatement exception when fed an invalid gene snp
    statement.

    """
    match = _STATEMENT_REGEX.match(statement)

    if not match:
        raise InvalidStatement
    (gene,
     base,
     reference,
     mutation,
     transcript_sequnce,
     bases,
     comment) = match.groups()
    return {"gene":                gene,
            "base":                int(base),
            "reference":           reference,
            "mutation":            mutation,
            "bases":               int(bases),
            "transcript_sequence": transcript_sequnce,
            "comment":             comment,
            }

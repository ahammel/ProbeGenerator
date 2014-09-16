"""Probe for an insertion/deletion event with the coordinate given relative to
a transcript.

"""
import re
import sys

from probe_generator import annotation, transcript
from probe_generator.variant import TranscriptVariant, GenomeVariant
from probe_generator.probe import AbstractProbe, InvalidStatement, NonFatalError

_STATEMENT_REGEX = re.compile(r"""
        \s*                # whitespace
        ([a-zA-Z0-9_./-]+) # gene name
        \s*
        :
        \s*
        c\.
        ([0-9]+)
        \s*
        (del\s*[acgtACGT]+|)
        \s*
        (ins\s*[acgtACGT]+|)
        \s*
        (\[trans\]|)
        \s*
        /
        \s*
        ([0-9]+)
        \s*
        (--.*|\s*)""", re.VERBOSE)


class GeneIndelProbe(AbstractProbe):
    """A probe specifying an insertion or deletion event starting at a base
    pair given relative to the start of a transcript.

    """
    _STATEMENT_SKELETON = ("{gene}:c.{base}{deletion}{insertion}"
                           "{transcript_sequence}/{bases}_"
                           "{transcript_name}_{chromosome}:{index_base}"
                           "{comment}")

    def __init__(self, *, variant, index, comment):
        self.variant = variant
        self.index = index
        self.comment = comment

    def __str__(self):
        reference = self.variant.reference
        mutation = self.variant.mutation

        deletion = 'del'+reference if reference else ''
        insertion = 'ins'+mutation if mutation else ''

        return self._STATEMENT_SKELETON.format(
            gene=self.variant.gene,
            base=self.index,
            deletion=deletion,
            insertion=insertion,
            transcript_sequence='[trans]' if self.variant.is_transcript else '',
            bases=len(self.variant),
            transcript_name=self.variant.transcript_name,
            chromosome=self.variant.chromosome,
            index_base=self.variant.coordinate,
            comment=self.comment)

    def get_ranges(self):
        return self.variant.sequence_ranges()

    @staticmethod
    def explode(statement, genome_annotation=None):
        probes = []

        if genome_annotation is None:
            genome_annotation = []
        specification = _parse(statement)

        if specification["transcript_sequence"] == '':
            variant_class = GenomeVariant
        else:
            variant_class = TranscriptVariant

        reference=specification["reference"]
        mutation=specification["mutation"]

        transcripts = annotation.lookup_gene(
            specification["gene"], genome_annotation)
        cached_coordinates = set()


        for txt in transcripts:
            base = specification["base"]
            try:
                indices = txt.transcript_range(
                    base, base+(len(reference) or 1))
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                if len(indices) != 1:
                    raise DiscontinuousIndelRegion(
                        "In statement {!r}, the requested mutation crosses an "
                        "exon/exon junction. Debug info {!s}".format(
                            statement, indices))
                index, = indices
                variant = variant_class(
                    transcript=txt,
                    index=index,
                    reference=reference,
                    mutation=mutation,
                    length=specification["bases"])
                if not (index, reference, mutation) in cached_coordinates:
                    cached_coordinates.add((index, reference, mutation))
                    probes.append(GeneIndelProbe(
                            variant=variant,
                            index=base,
                            comment=specification["comment"]))
        return probes


def _parse(statement):
    """Return a partial GeneIndelProbe specification given a probe statement.

    Raises an InvalidStatement exception when fed an invalid gene snp
    statement.

    """
    match = _STATEMENT_REGEX.match(statement)

    if not match:
        raise InvalidStatement

    (gene,
     index,
     deletion,
     insertion,
     transcript_sequence,
     bases,
     comment) = match.groups()

    deletion = "".join(deletion.split()) # Remove whitespace
    insertion = "".join(insertion.split())

    return {"gene":                gene,
            "base":                int(index),
            "deletion":            deletion,
            "insertion":           insertion,
            "bases":               int(bases),
            "comment":             comment,
            "transcript_sequence": transcript_sequence,
            "reference":           deletion.lstrip("del"),
            "mutation":            insertion.lstrip("ins")}



class DiscontinuousIndelRegion(NonFatalError):
    """Raised when the transcript region requested by an indel probe crosses an
    exon/exon junction.

    Will be removed if and when we support exon-crossing events.

    """

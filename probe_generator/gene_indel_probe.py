"""Probe for an insertion/deletion event with the coordinate given relative to
a transcript.

"""
import re
import sys

from probe_generator import annotation, transcript
from probe_generator.sequence import SequenceRange
from probe_generator.sequence import reverse_complement
from probe_generator.probe import AbstractProbe, InvalidStatement

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
        /
        \s*
        ([0-9]+)
        \s*
        (--.*|\s*)""", re.VERBOSE)


class GeneIndelProbe(AbstractProbe):
    """A probe specifying an insertion or deletion event starting at a base
    pair given relative to the start of a transcript.

    """
    _STATEMENT_SKELETON = ("{gene}:c.{base}{deletion}{insertion}/{bases}_"
                           "{transcript}_{chromosome}:{index_base}")
    def get_ranges(self):
        chromosome, start, _end, _, _ = self._spec["index"]
        bases = self._spec['bases']

        reference_bases = len(self._spec["reference"])
        mutation_bases = len(self._spec["mutation"])

        left_buffer = bases // 2
        if bases % 2 == 0:
            left_buffer -= 1
        right_buffer = bases - left_buffer - mutation_bases

        return (
            SequenceRange(chromosome,
                          start-left_buffer,
                          start),
            SequenceRange(chromosome,
                          start,
                          start+reference_bases,
                          mutation=True),
            SequenceRange(chromosome,
                          start+reference_bases,
                          start+reference_bases+right_buffer))

    @staticmethod
    def explode(statement, genome_annotation):
        probes = []

        if genome_annotation is None:
            genome_annotation = []
        partial_spec = _parse(statement)
        transcripts = annotation.lookup_gene(
            partial_spec["gene"], genome_annotation)
        cached_coordinates = set()
        for txt in transcripts:
            if txt.plus_strand:
                base = partial_spec["base"]
            else:
                base = partial_spec["base"] - 2
                partial_spec["mutation"] = reverse_complement(
                    partial_spec["mutation"])
                partial_spec["reference"] = reverse_complement(
                    partial_spec["reference"])
            try:
                index = txt.nucleotide_index(base)
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                chromosome = txt.chromosome
                if not (chromosome, index) in cached_coordinates:
                    cached_coordinates.add((chromosome, index))
                    spec = dict(partial_spec,
                                chromosome=chromosome,
                                transcript=txt.name,
                                index=index,
                                index_base=index.start+1)
                    probes.append(GeneIndelProbe(spec))
        return probes

def _parse(statement):
    match = _STATEMENT_REGEX.match(statement)

    if not match:
        raise InvalidStatement

    (gene,
     index,
     deletion,
     insertion,
     bases,
     comment) = match.groups()

    deletion = "".join(deletion.split()) # Remove whitespace
    insertion = "".join(insertion.split())

    return {"gene":       gene,
            "base":       int(index),
            "deletion":   deletion,
            "insertion":  insertion,
            "bases":      int(bases),
            "comment":    comment,
            "reference":  deletion.lstrip("del"),
            "mutation":   insertion.lstrip("ins")}


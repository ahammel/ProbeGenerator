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
                           "{transcript_name}_{chromosome}:{index_base}")
    def get_ranges(self):
        bases = self._spec['bases']

        mutation_bases = len(self._spec["mutation"])

        total_buffer = bases - mutation_bases
        left_buffer = total_buffer // 2
        right_buffer = total_buffer - left_buffer

        if self._spec['transcript_sequence']:
            return self._get_ranges_transcript(left_buffer, right_buffer)
        else:
            return self._get_ranges_genome(left_buffer, right_buffer)

    def _get_ranges_transcript(self, left_buffer, right_buffer):
        """Return the SequenceRange representation of the variant buffered by
        bases taken from the transcript sequence of the gene.

        E.g., if the variant is at the end of an exon on the plus strand, the
        right_buffer bases will be taken from the next exon.

        """
        chromosome, start, _end, _, _ = self._spec["index"]
        reference_bases = len(self._spec['reference'])
        txt = self._spec['transcript']
        base = self._spec['base']

        if not txt.plus_strand:
            left_buffer, right_buffer = right_buffer, left_buffer

        sequence = (
            txt.transcript_range(base-left_buffer, base) +
            [SequenceRange(chromosome,
                           start,
                           start+reference_bases,
                           mutation=True,
                           reverse_complement=not txt.plus_strand)] +
            txt.transcript_range(base+reference_bases,
                                 base+reference_bases+right_buffer))

        if txt.plus_strand:
            return sequence
        else:
            return reversed(sequence)

    def _get_ranges_genome(self, left_buffer, right_buffer):
        """Return the SequenceRagne representation of the variant buffered by
        bases taken from the reference genome seqeunce.

        """
        chromosome, start, _end, _, _ = self._spec["index"]
        reference_bases = len(self._spec['reference'])

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
                                transcript=txt,
                                transcript_name=txt.name,
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


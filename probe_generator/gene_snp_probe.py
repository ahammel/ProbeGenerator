"""Probe for an SNP mutation based on the index of a nucleotide in a coding
sequence.

"""
import re
import sys

from probe_generator import annotation, transcript
from probe_generator.sequence import SequenceRange
from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence import reverse_complement

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
                           "{transcriptome_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{index_base}{comment}")

    def get_ranges(self):
        chromosome, start, end, _, _ = self._spec['index']
        bases = self._spec['bases']

        mutation_bases = len(self._spec["mutation"])

        left_buffer = bases // 2
        if bases % 2 == 0:
            left_buffer -= 1
        right_buffer = bases - left_buffer - mutation_bases

        if self._spec["transcriptome_sequence"]:
            return self._get_ranges_transcript(left_buffer, right_buffer)
        else:
            return self._get_ranges_genome(left_buffer, right_buffer)

    def _get_ranges_transcript(self, left_buffer, right_buffer):
        """Return the SequenceRange representation of the variant buffered by
        bases taken from the transcript sequence of the gene.

        E.g., if the variant is at the end of an exon on the plus strand, the
        right_buffer bases will be taken from the next exon.

        """
        chromosome, start, end, _, _ = self._spec['index']
        txt = self._spec['transcript']
        base, = self._spec["base"],
        return (
            txt.transcript_range(base-left_buffer, base) +
            [SequenceRange(chromosome,
                           start,
                           end,
                           mutation=True,
                           reverse_complement= self._spec['strand'] == '-')] +
            txt.transcript_range(base+1, base+1+right_buffer))

    def _get_ranges_genome(self, left_buffer, right_buffer):
        """Return the SequenceRagne representation of the variant buffered by
        bases taken from the reference genome seqeunce.

        """
        chromosome, start, end, _, _ = self._spec['index']
        return (
            SequenceRange(chromosome,
                          start-left_buffer,
                          start),
            SequenceRange(chromosome,
                          start,
                          end,
                          mutation=True,
                          reverse_complement=self._spec['strand'] == '-'),
            SequenceRange(chromosome,
                          end,
                          end+right_buffer))

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
        partial_spec = _parse(statement)
        transcripts = annotation.lookup_gene(
            partial_spec["gene"], genome_annotation)
        cached_coordinates = set()
        for txt in transcripts:
            base = partial_spec["base"]
            if not txt.plus_strand:
                partial_spec["mutation"] = reverse_complement(
                    partial_spec["mutation"])
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
                                strand='+' if txt.plus_strand else '-',
                                chromosome=chromosome,
                                transcript=txt,
                                transcript_name=txt.name,
                                index=index,
                                index_base=index.start+1)
                    probes.append(GeneSnpProbe(spec))
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
     transcriptome_sequnce,
     bases,
     comment) = match.groups()
    return {"gene":                   gene,
            "base":                   int(base),
            "reference":              reference,
            "mutation":               mutation,
            "bases":                  int(bases),
            "transcriptome_sequence": transcriptome_sequnce,
            "comment":                comment,
            }

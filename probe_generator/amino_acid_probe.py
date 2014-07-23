"""Probe based on a gene name and amino acid change.

"""
import re
import sys
import itertools

from probe_generator import annotation, transcript
from probe_generator.sequence import reverse_complement, SequenceRange
from probe_generator.probe import AbstractProbe, InvalidStatement

_STATEMENT_REGEX = re.compile("""
        \s*                                           # whitespace
        ([A-Za-z0-9_./-]+)                            # gene name
        \s*
        :
        \s*
        ([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*]) # IUPAC amino acid code
        \s*
        ([0-9]+)                                      # codon number
        \s*
        ([ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx*])
        \s*
        (\[trans\]|)                                  # transcript-only-sequence
        \s*
        /
        \s*
        ([0-9]+)                                      # number of base pairs
        \s*
        (--.*|\s*)                                    # comment
        """, re.VERBOSE)

_DNA_CODON_TABLE = {
    'A': ("GCT", "GCC", "GCA", "GCG"),
    'C': ("TGT", "TGC"),
    'D': ("GAT", "GAC"),
    'E': ("GAA", "GAG"),
    'F': ("TTT", "TTC"),
    'G': ("GGT", "GGC", "GGA", "GGG"),
    'H': ("CAT", "CAC"),
    'I': ("ATT", "ATC", "ATA"),
    'K': ("AAA", "AAG"),
    'L': ("CTT", "CTC", "CTA", "CTG", "TTA", "TTG"),
    'M': ("ATG",),
    'N': ("AAT", "AAC"),
    'P': ("CCT", "CCC", "CCA", "CCG"),
    'Q': ("CAA", "CAG"),
    'R': ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    'S': ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    'T': ("ACT", "ACC", "ACA", "ACG"),
    'V': ("GTT", "GTC", "GTA", "GTG"),
    'W': ("TGG",),
    'Y': ("TAT", "TAC"),
    '*': ("TAA", "TAG", "TGA"),
    'X': tuple(''.join(bases) for bases in itertools.product("ACGT", repeat=3)),
    }


class AminoAcidProbe(AbstractProbe):
    """Probe for a single nucleotide polymorphism from an amino acid change at
    a codon.

    """
    _STATEMENT_SKELETON = ("{gene}:{reference_aa}{codon}{mutation_aa}"
                           "({reference_bases}>{mutation_bases})"
                           "{transcript_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{coordinate}{comment}")

    def get_ranges(self):
        chromosome, start, end, _, _ = self._spec['index']
        bases = self._spec['bases']

        mutation_bases = len(self._spec["mutation"])
        left_buffer = bases // 2 - 1
        if bases % 2 == 0:
            left_buffer -= 1
        right_buffer = bases - left_buffer - mutation_bases

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
        txt = self._spec['transcript']
        codon = self._spec['codon']
        codon_start, codon_end = (codon-1)*3+1, (codon*3)+1
        chromosome, start, end, _, _ = self._spec['index']

        return (txt.transcript_range(codon_start-left_buffer, codon_start) +
                [SequenceRange(chromosome,
                               start,
                               end ,
                               mutation=True)] +
                txt.transcript_range(codon_end, codon_end+right_buffer))

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
                          mutation=True),
            SequenceRange(chromosome,
                          end,
                          end+right_buffer))

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given a probe statement and a genome annotation, yield a sequence of
        AminoAcid probes matching the statement.

        If more than one probe has identical genomic coordinates, only the
        first is returned.

        """
        probes = []

        if genome_annotation is None:
            genome_annotation = []
        partial_spec = _parse(statement)
        transcripts = annotation.lookup_gene(
            partial_spec['gene'],
            genome_annotation)
        coordinate_cache = set()
        for txt in transcripts:
            try:
                index = txt.codon_index(partial_spec['codon'])
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                chromosome = txt.chromosome
                if (chromosome, index) not in coordinate_cache:
                    coordinate_cache.add((chromosome, index))
                    reference_codons = _DNA_CODON_TABLE[
                        partial_spec['reference_aa'].upper()]
                    mutation_codons = _DNA_CODON_TABLE[
                        partial_spec['mutation_aa'].upper()]
                    for reference_codon, mutation_codon in itertools.product(
                        reference_codons, mutation_codons):
                        if (mutation_codon not in
                            _DNA_CODON_TABLE[partial_spec["reference_aa"]]):
                            if txt.plus_strand:
                                mutation = mutation_codon
                                reference = reference_codon
                                strand = '+'
                            else:
                                mutation = reverse_complement(mutation_codon)
                                reference = reverse_complement(reference_codon)
                                strand = '-'
                            spec = dict(partial_spec,
                                        index=index,
                                        chromosome=chromosome,
                                        strand=strand,
                                        transcript=txt,
                                        transcript_name=txt.name,
                                        reference=reference,
                                        reference_bases=reference_codon,
                                        mutation_bases=mutation_codon,
                                        mutation=mutation,
                                        coordinate=index.start+1)
                            # 'mutation_bases' is displayed in the probe's string,
                            # while 'mutation' is used internally to determine the
                            # sequence.
                            #
                            # They are different sequences if the transcript is on
                            # the '-' strand.
                            probes.append(AminoAcidProbe(spec))
        return probes


def _parse(statement):
    """Return a partial AminoAcidProbe given a probe statement.

    Raises an InvalidStatement exception when the statement does not match the
    format of an amino acid probe.

    """
    match = _STATEMENT_REGEX.match(statement)
    if not match:
        raise InvalidStatement

    (gene,
     reference_aa,
     codon,
     mutation_aa,
     transcript_sequence,
     bases,
     comment) = match.groups()

    return {"gene":                gene,
            "reference_aa":        reference_aa,
            "codon":               int(codon),
            "mutation_aa":         mutation_aa,
            "bases":               int(bases),
            "transcript_sequence": transcript_sequence,
            "comment":             comment}

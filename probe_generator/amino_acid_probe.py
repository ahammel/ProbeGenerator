"""Probe based on a gene name and amino acid change.

"""
import re
import sys
import itertools

from probe_generator import annotation
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
        ([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*])
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
}


class AminoAcidProbe(AbstractProbe):
    """Probe for a single nucleotide polymorphism from an amino acid change at
    a codon.

    """
    _STATEMENT_SKELETON = ("{gene}:{reference_aa}{codon}{mutation_aa}"
                           "({reference_bases}>{mutation_bases})/"
                           "{bases}_{transcript}_{chromosome}:{index}{comment}")

    def get_ranges(self):
        index = self._spec['index']
        bases = self._spec['bases']
        chromosome = self._spec['chromosome']

        left_buffer = bases // 2
        if bases % 2 == 0:
            left_buffer -= 1
        right_buffer = bases - left_buffer

        return (
            SequenceRange(chromosome,
                          index-left_buffer,
                          index-1),
            SequenceRange(chromosome,
                          index-1,
                          index+2,
                          mutation=True),
            SequenceRange(chromosome,
                          index+2,
                          index+right_buffer))

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given a probe statement and a genome annotation, yield a sequence of
        AminoAcid probes matching the statement.

        If more than one probe has identical genomic coordinates, only the
        first is returned.

        """
        if genome_annotation is None:
            genome_annotation = []
        partial_spec = _parse(statement)
        transcripts = annotation.lookup_gene(
            partial_spec['gene'],
            genome_annotation)
        coordinate_cache = set()
        for transcript in transcripts:
            try:
                index = annotation.codon_index(partial_spec['codon'], transcript)
            except annotation.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                chromosome = transcript['chrom'].lstrip('chr')
                if (chromosome, index) not in coordinate_cache:
                    coordinate_cache.add((chromosome, index))
                    reference_codons = _DNA_CODON_TABLE[
                        partial_spec['reference_aa'].upper()]
                    mutation_codons = _DNA_CODON_TABLE[
                        partial_spec['mutation_aa'].upper()]
                    for reference_codon, mutation_codon in itertools.product(
                        reference_codons, mutation_codons):
                        if transcript['strand'] == '+':
                            mutation = mutation_codon
                            reference = reference_codon
                        else:
                            mutation = reverse_complement(mutation_codon)
                            reference = reverse_complement(reference_codon)
                        spec = dict(partial_spec,
                                    index=index,
                                    chromosome=chromosome,
                                    strand=transcript['strand'],
                                    transcript=transcript['name'],
                                    reference=reference,
                                    reference_bases=reference_codon,
                                    mutation_bases=mutation_codon,
                                    mutation=mutation)
                        # 'mutation_bases' is displayed in the probe's string,
                        # while 'mutation' is used internally to determine the
                        # sequence.
                        #
                        # They are different sequences if the transcript is on
                        # the '-' strand.
                        yield AminoAcidProbe(spec)


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
     bases,
     comment) = match.groups()

    return {"gene":          gene,
            "reference_aa":  reference_aa,
            "codon":         int(codon),
            "mutation_aa":   mutation_aa,
            "bases":         int(bases),
            "comment":       comment}

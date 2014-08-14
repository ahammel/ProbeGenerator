"""Probe based on a gene name and amino acid change.

"""
import re
import sys
import itertools

from probe_generator import annotation, transcript
from probe_generator.variant import TranscriptVariant, GenomeVariant
from probe_generator.probe import AbstractProbe, InvalidStatement

_STATEMENT_REGEX = re.compile(r"""
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

_AMINO_ACID_TABLE = {
    codon: amino_acid
    for amino_acid, codons in _DNA_CODON_TABLE.items() if amino_acid != 'X'
    for codon in codons
    }


class AminoAcidProbe(AbstractProbe):
    """Probe for a single nucleotide polymorphism from an amino acid change at
    a codon.

    """
    _STATEMENT_SKELETON = ("{gene}:{reference_aa}{codon}{mutation_aa}"
                           "({reference}>{mutation})"
                           "{transcript_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{coordinate}{comment}")

    def __init__(self, *, variant, index, comment):
        self.variant = variant
        self.index = index
        self.comment = comment

    def __str__(self):
        return self._STATEMENT_SKELETON.format(
            gene=self.variant.gene,
            reference_aa=_AMINO_ACID_TABLE[self.variant.reference],
            codon=self.index,
            mutation_aa=_AMINO_ACID_TABLE[self.variant.mutation],
            reference=self.variant.reference,
            mutation=self.variant.mutation,
            transcript_sequence='[trans]' if self.variant.is_transcript else '',
            bases=len(self.variant),
            transcript_name=self.variant.transcript_name,
            chromosome=self.variant.chromosome,
            coordinate=self.variant.coordinate,
            comment=self.comment)

    def get_ranges(self):
        return self.variant.sequence_ranges()

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

        specification = _parse(statement)
        transcripts = annotation.lookup_gene(
            specification['gene'],
            genome_annotation)

        if specification["transcript_sequence"] == '':
            variant_class = GenomeVariant
        else:
            variant_class = TranscriptVariant

        reference_aa = specification['reference_aa'].upper()
        mutation_aa = specification['mutation_aa'].upper()
        reference_codons = _DNA_CODON_TABLE[reference_aa]
        mutation_codons = _DNA_CODON_TABLE[mutation_aa]

        index = specification["index"]

        coordinate_cache = set()
        for txt, reference, mutation in itertools.product(
            transcripts, reference_codons, mutation_codons):
            try:
                sequence_range_index = txt.codon_index(specification['index'])
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                variant = variant_class(
                    transcript=txt,
                    index=sequence_range_index,
                    reference=reference,
                    mutation=mutation,
                    length=specification["bases"])
                if _AMINO_ACID_TABLE[mutation] != reference_aa:
                    coordinate_cache.add((variant, index))
                    probe = AminoAcidProbe(
                        variant=variant,
                        index=index,
                        comment=specification["comment"])
                    probes.append(probe)
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
     index,
     mutation_aa,
     transcript_sequence,
     bases,
     comment) = match.groups()

    return {"gene":                gene,
            "reference_aa":        reference_aa,
            "index":               int(index),
            "mutation_aa":         mutation_aa,
            "bases":               int(bases),
            "transcript_sequence": transcript_sequence,
            "comment":             comment}

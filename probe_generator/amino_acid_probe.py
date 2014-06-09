"""Probe based on a gene name and amino acid change.

"""
import re
import sys

from probe_generator import annotation, probe, reference
from probe_generator.sequence import reverse_complement

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

_STATEMENT_SKELETON = ("{gene}:{reference}{codon}{mutation}({mutation_bases})"
                       "/{bases}_{transcript}_{chromosome}:{index}{comment}")

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


class AminoAcidProbe(object):
    """Probe for a single nucleotide polymorphism from an amino acid change at
    a codon.

    """
    def __init__(self, specification):
        self._spec = specification

    def __str__(self):
        return _STATEMENT_SKELETON.format(**self._spec)

    def sequence(self, genome):
        """Return the probe sequence given a genome sequence.

        """
        if self._spec["strand"] == '+':
            index = self._spec["index"] + 1
        else:
            index = self._spec["index"] - 1
        # The index must be adjusted so that the codon starts in the centre of
        # the probe.
        #
        # TODO: It would probably be better to have the second base of the
        # codon cetnred in the probe when possible.
        start, end = annotation.get_bases(
            self._spec['bases'],
            index)
        raw_bases = reference.bases(
            genome,
            self._spec["chromosome"],
            start,
            end)
        return self._mutate(raw_bases)

    def _mutate(self, bases):
        codon_index = (self._spec["bases"] // 2)
        spec_codons = _DNA_CODON_TABLE[self._spec["reference"]]
        if self._spec['strand'] == '+':
            reference_codon = bases[codon_index-2:codon_index+1].upper()
            mutation_bases = self._spec["mutation_bases"]
        else:
            reference_codon = reverse_complement(
                bases[codon_index-2:codon_index+1].upper())
            mutation_bases = reverse_complement(
                self._spec["mutation_bases"])
        if reference_codon not in spec_codons:
            raise CodonMismatch(
                "Reference sequence {!r} does not match requested amino "
                "acid {!r} in probe {}".format(
                    reference_codon,
                    self._spec["reference"],
                    self))
        return (bases[:codon_index-2] +
                mutation_bases        +
                bases[codon_index+1:])

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
                    codons = _DNA_CODON_TABLE[partial_spec['mutation'].upper()]
                    for codon in codons:
                        spec = dict(partial_spec,
                                    index=index,
                                    chromosome=chromosome,
                                    strand=transcript['strand'],
                                    transcript=transcript['name'],
                                    mutation_bases=codon)
                        yield AminoAcidProbe(spec)


def _parse(statement):
    """Return a partial AminoAcidProbe given a probe statement.

    Raises an InvalidStatement exception when the statement does not match the
    format of an amino acid probe.

    """
    match = _STATEMENT_REGEX.match(statement)
    if not match:
        raise probe.InvalidStatement

    (gene,
     reference_aa,
     codon,
     mutation_aa,
     bases,
     comment) = match.groups()

    return {"gene":       gene,
            "reference":  reference_aa,
            "codon":      int(codon),
            "mutation":   mutation_aa,
            "bases":      int(bases),
            "comment":    comment}


class CodonMismatch(probe.NonFatalError):
    """Raised when the codon at the genome reference does not match the codon
    specified in the statement.

    """

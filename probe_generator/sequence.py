"""Nucleotide sequence manipulation utility functions.

"""
import itertools

_COMPLEMENT = str.maketrans('acgtACGT', 'tgcaTGCA')

DNA_CODON_TABLE = {
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

AMINO_ACID_TABLE = {
    codon: amino_acid
    for amino_acid, codons in DNA_CODON_TABLE.items() if amino_acid != 'X'
    for codon in codons
    }


def complement(string):
    """Return the complement of a string of nucleotides.

    """
    return string.translate(_COMPLEMENT)


def reverse_complement(string):
    """Return the reverse-complement of a string of nucleotides.

    """
    return ''.join(reversed(complement(string)))


def translate(sequence):
    """Translate a DNA sequence to a peptide.

    """
    assert len(sequence) % 3 == 0
    return ''.join(AMINO_ACID_TABLE[sequence.upper()[i: i+3]]
                   for i in range(0, len(sequence), 3))


def reverse_translate(peptide):
    """Yield all the DNA sequences that translate to the original
    peptide.

    """
    sequences = (DNA_CODON_TABLE[amino_acid] for amino_acid in peptide.upper())
    for seq in itertools.product(*sequences):
        yield ''.join(seq)

"""Nucleotide seqeunce manipulation library for ProbeGenerator.

Provides complement() and reverse_complemement() functions.

"""

_COMPLEMENT = str.maketrans('acgtACGT', 'tgcaTGCA')


def complement(string):
    """Return the complement of a string of nucleotides.

    """
    return string.translate(_COMPLEMENT)


def reverse_complement(string):
    """Return the reverse-complement of a string of nucleotides.

    """
    return ''.join(reversed(complement(string)))

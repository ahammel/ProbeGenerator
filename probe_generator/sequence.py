"""Nucleotide sequence manipulation library for ProbeGenerator.

Provides complement() and reverse_complemement() functions and the
SequenceRange object.

"""
from collections import namedtuple

_COMPLEMENT = str.maketrans('acgtACGT', 'tgcaTGCA')


class SequenceRange(namedtuple("SequenceRange",
                               ["chromosome",
                                "start",
                                "end",
                                "reverse_complement",
                                "mutant"])):
    """Data object for specifying a range of base pairs to be extracted from
    the genome.

    The 'start' and 'end' fields are 0-indexed, left-inclusive, right-exclusive
    (like slices in Python).

    The optional 'mutant' flag is used to mark a range which will be replaced
    with a different sequence by a Probe object.

    """
    # I've subclassed the namedtuple rather than using it raw so I can use
    # keyword-only args and default values.

    __slots__ = () # Performance hack

    def __new__(self, chromosome, start, end, *,
                # Keyword only arguments:
                reverse_complement=False, mutant=False):
        return super().__new__(self,
                               chromosome,
                               start,
                               end,
                               reverse_complement,
                               mutant)


def complement(string):
    """Return the complement of a string of nucleotides.

    """
    return string.translate(_COMPLEMENT)


def reverse_complement(string):
    """Return the reverse-complement of a string of nucleotides.

    """
    return ''.join(reversed(complement(string)))

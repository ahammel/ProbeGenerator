"""Provides the SequenceRange object.  

"""
from probe_generator import sequence

from collections import namedtuple


class SequenceRange(namedtuple("SequenceRange",
                               ["chromosome",
                                "start",
                                "end",
                                "reverse_complement",
                                "reference",
                                "mutation"])):
    """Data object for specifying a range of base pairs to be extracted from
    the genome.

    The 'start' and 'end' fields are 0-indexed, left-inclusive, right-exclusive
    (like slices in Python).

    The optional 'mutation' flag is used to mark a range which will be replaced
    with a different sequence by a Probe object.

    """
    __slots__ = ()
    # This is a performance hack which reduces the memory footprint of the
    # object.

    def __new__(cls, chromosome, start, end, *,
                reverse_complement=False, reference=None, mutation=None):
        """As in the standard `namedtuple` __new__ method, but
        `reverse_complement` and `mutation` are keyword-only arguments with
        default values.

        If reverse_complement is True and the mutation is not None, the
        mutation is reverse-complemented.

        """
        if reverse_complement and mutation is not None:
            the_mutation = sequence.reverse_complement(mutation)
        else:
            the_mutation = mutation
        return super().__new__(cls,
                               chromosome,
                               start,
                               end,
                               reverse_complement,
                               reference,
                               the_mutation)

    def concat(self, other):
        """Return a new SequenceRange object representing the combined genomic
        region of the two SequenceRanges.

        Raise a ValueError if the two SequenceRanges are not adjacent.

        """
        if not self.adjacent(other):
            raise ValueError(
                "Cannot concatenate non-adjacent "
                "SequenceRange objects {} and {}".format(self, other))
        if self.end == other.start:
            return SequenceRange(self.chromosome, self.start, other.end,
                                 mutation=self.mutation,
                                 reverse_complement=self.reverse_complement)
        elif self.start == other.end:
            return SequenceRange(self.chromosome, other.start, self.end,
                                 mutation=self.mutation,
                                 reverse_complement=self.reverse_complement)
        else:
            assert False, "unreachable"

    def adjacent(self, other):
        """Return True if the SequenceRange object ends where the other starts
        or vice-versa.

        The chromosome, mutation, and reverse_complement status must be equal
        for both SequenceRange objects.

        """
        return (((self.chromosome         == other.chromosome)          and
                 (self.mutation           == other.mutation)            and
                 (self.reverse_complement == other.reverse_complement)) and
                ((self.start              == other.end)                 or
                 (self.end                == other.start)))

    def condense(self, *others):
        """Given an iterable of SequenceRange objects, recursively concatenate
        ranges which are adjacent to `self` until no two consecutive ranges in
        the list are adjacent.

        """
        ranges = []
        chunk = self
        for other in others:
            if chunk.adjacent(other):
                chunk = chunk.concat(other)
            else:
                ranges.append(chunk)
                chunk = other
        ranges.append(chunk)
        return ranges

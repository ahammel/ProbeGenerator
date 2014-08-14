"""Objects representing substitution variants (indels, SNPs, multiple
substitutions, etc.)

"""
from probe_generator.sequence import SequenceRange

class AbstractVariant(object):
    """Super-class for Variant objects.

    The *Variant collection of classes represent arbitrary-length sequence
    substitutions. Common uses include one-for-one substitutions (SNPs),
    three-for-three (codon substitutions) and n-for-n (indels).

    The super-class enforces the __init__ method and provides __len__ (which
    just returns the length provided at init time).

    `length` is the sum of the length of the variant and the surrounding buffer
    sequence.

    """
    def __init__(self, *, transcript, index, reference, mutation, length):
        self.transcript = transcript
        self.index = index
        self.reference = reference
        self.mutation = mutation
        self._length = length

        self.gene = self.transcript.gene_id
        self.transcript_name = self.transcript.name
        self.chromosome = self.transcript.chromosome
        self.coordinate = self.index.start + 1

    def __len__(self):
        return self._length


class TranscriptVariant(AbstractVariant):
    """A substitution variant using buffer sequence from the surrounding
    transcript.

    """
    is_transcript = True

    def sequence_ranges(self):
        chromosome, start, _end, _, _ = self.index
        reference_length = len(self.reference)

        reference_length = len(self.reference)
        mutation_length = len(self.mutation)

        total_buffer = self.length - mutation_length
        left_buffer = total_buffer // 2
        right_buffer = total_buffer - left_buffer

        base = self.transcript.base_index(self.index)

        if not self.transcript.plus_strand:
            left_buffer, right_buffer = right_buffer, left_buffer

        sequence = (
            self.transcript.transcript_range(base-left_buffer, base) +
            [SequenceRange(chromosome,
                           start,
                           start+reference_length,
                           mutation=self.mutation,
                           reverse_complement=not self.transcript.plus_strand)] +
            self.transcript.transcript_range(base+reference_length,
                                             base+reference_length+right_buffer)
            )

        if self.transcript.plus_strand:
            return sequence
        else:
            return reversed(sequence)


class GenomeVariant(AbstractVariant):
    """A substitution variant using buffer sequence from the surrounding genome
    sequence.

    """
    is_transcript = False

    def sequence_ranges(self):
        chromosome, start, _end, _, _ = self.index

        reference_length = len(self.reference)
        mutation_length = len(self.mutation)

        total_buffer = self.length - mutation_length
        left_buffer = total_buffer // 2
        right_buffer = total_buffer - left_buffer

        return (
            SequenceRange(chromosome,
                          start-left_buffer,
                          start),
            SequenceRange(chromosome,
                          start,
                          start+reference_length,
                          mutation=self.mutation,
                          reverse_complement=not self.transcript.plus_strand),
            SequenceRange(chromosome,
                          start+reference_length,
                          start+reference_length+right_buffer))

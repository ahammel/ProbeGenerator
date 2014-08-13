"""Objects representing substitution variants (indels, SNPs, multiple
substitutions, etc.)

"""
from probe_generator.sequence import SequenceRange

class AbstractVariant(object):
    def __init__(self, *, transcript, index, reference, mutation, bases):
        self.transcript = transcript
        self.index = index
        self.reference = reference
        self.mutation = mutation
        self.bases = bases

        self.gene = self.transcript.gene_id
        self.transcript_name = self.transcript.name
        self.chromosome = self.transcript.chromosome
        self.coordinate = self.index.start + 1


class TranscriptVariant(AbstractVariant):
    """A substitution variant using buffer sequence from the surrounding
    transcript.

    """
    is_transcript = True

    def sequence_ranges(self):
        chromosome, start, _end, _, _ = self.index
        reference_bases = len(self.reference)

        reference_bases = len(self.reference)
        mutation_bases = len(self.mutation)

        total_buffer = self.bases - mutation_bases
        left_buffer = total_buffer // 2
        right_buffer = total_buffer - left_buffer

        base = self.transcript.base_index(self.index)

        if not self.transcript.plus_strand:
            left_buffer, right_buffer = right_buffer, left_buffer

        sequence = (
            self.transcript.transcript_range(base-left_buffer, base) +
            [SequenceRange(chromosome,
                           start,
                           start+reference_bases,
                           mutation=self.mutation,
                           reverse_complement=not self.transcript.plus_strand)] +
            self.transcript.transcript_range(base+reference_bases,
                                 base+reference_bases+right_buffer))

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

        reference_bases = len(self.reference)
        mutation_bases = len(self.mutation)

        total_buffer = self.bases - mutation_bases
        left_buffer = total_buffer // 2
        right_buffer = total_buffer - left_buffer

        return (
            SequenceRange(chromosome,
                          start-left_buffer,
                          start),
            SequenceRange(chromosome,
                          start,
                          start+reference_bases,
                          mutation=self.mutation,
                          reverse_complement=not self.transcript.plus_strand),
            SequenceRange(chromosome,
                          start+reference_bases,
                          start+reference_bases+right_buffer))

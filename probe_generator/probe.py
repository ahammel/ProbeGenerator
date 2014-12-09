"""Shared objects for probe classes.

"""
from abc import ABCMeta, abstractmethod

from probe_generator import reference
from probe_generator.exceptions import NonFatalError


class AbstractProbe(object, metaclass=ABCMeta):
    """Super-class for Probe objects.

    Subclasses provide the _STATEMENT_SKELETON property, an 'explode' static
    method, and a 'get_ranges' method. The '__init__', '__str__', and
    'sequence' methods are mixed-in.

    """
    variant = NotImplemented # provided by children

    def sequence(self, genome):
        """Return the sequence of the probe given a reference genome object
        using the SequenceRange objects returned by the get_ranges method.

        Raises a MissingChromosome exception (non-fatal) when the chromosome is
        not present in the reference genome.

        Raises a ReferenceMismatch (non-fatal) when the reference sequence
        specified by the probe does not match the sequence taken from the
        reference genome.

        Raises a NonContainedRange error (fatal) when the range requested falls
        outside the chromosome.

        """
        ranges = self.get_ranges()
        return ''.join(self._get_sequence(seq_range, genome)
                       for seq_range in ranges)

    def _get_sequence(self, seq_range, genome):
        """Return a string of base pairs given a SequenceRange object and a
        reference genome.

        If the SequenceRange is a mutation, check that the reference base
        specified by the probe matches the reference sequence and return the
        mutation bases specified by the probe.

        """
        bases = reference.bases(seq_range, genome)
        if seq_range.reference is not None:
            self._assert_reference_matches(seq_range.reference, bases)
        if seq_range.mutation is not None:
            return seq_range.mutation
        else:
            return bases

    def _assert_reference_matches(self, expected, actual):
        """If the 'bases' string is not the same as the reference bases
        specified by the probe, raise a ReferenceMismatch exception.

        Raises a NotImplementedError if the probe doesn't support sequence
        variants (i.e., if the probe is not for sequence variant mutations,
        this method should never be called).

        """
        if not expected.lower() == actual.lower():
            raise ReferenceMismatch(
                "In probe '{!s}': expected sequence {!r} but got {!r}" .format(
                    self, expected, actual))

    @abstractmethod
    def get_ranges(self):
        """Return an iterable of SequenceRange objects representing the bases
        to be collected.

        """

    @staticmethod
    @abstractmethod
    def explode(statement, genome_annotation=None):
        """Return a list of probes from a statement and, optionally, a genome
        annotation.

        """


class ReferenceMismatch(NonFatalError):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """


class InvalidStatement(NonFatalError):
    """Raised when a probe statement cannot be parsed.

    Note that this is a fatal error.

    """

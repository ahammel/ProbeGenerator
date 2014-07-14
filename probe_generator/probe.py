"""Shared objects for probe classes.

"""
import abc

from probe_generator import reference
from probe_generator.exceptions import NonFatalError


class AbstractProbe(metaclass=abc.ABCMeta):
    """Super-class for Probe objects.

    Subclasses provide the _STATEMENT_SKELETON property, an 'explode' static
    method, and a 'get_ranges' method. The '__init__', '__str__', and
    'sequence' methods are mixed-in.

    """
    def __init__(self, specification, *rows):
        self._spec = specification
        self._rows = rows

    def __str__(self):
        return self._STATEMENT_SKELETON.format(**self._spec)

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
        if seq_range.mutation:
            self._assert_reference_matches(bases)
            return self._spec['mutation']
        else:
            return bases

    def _assert_reference_matches(self, bases):
        """If the 'bases' string is not the same as the reference bases
        specified by the probe, raise a ReferenceMismatch exception.

        """
        if not self._spec['reference'].lower() == bases.lower():
            raise ReferenceMismatch(
                    "Reference sequence {!r} does not match requested mutation "
                    "{!r} => {!r}".format(bases,
                                     self._spec["reference"],
                                     self._spec["mutation"]))

    @abc.abstractproperty
    def _STATEMENT_SKELETON(self):
        """The template of the __str__ property of the probe.

        The internal _spec property of an instance is passed to
        _STATEMENT_SKELETON.format to produce the probe's string.

        """

    @abc.abstractmethod
    def get_ranges(self):
        """Return an iterable of SequenceRange objects representing the bases
        to be collected.

        """

    @abc.abstractstaticmethod
    def explode(statement, genome_annotation=None):
        """Return a list of probes from a statement and, optionally, a genome
        annotation.

        """


class ReferenceMismatch(NonFatalError):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """


class InvalidStatement(Exception):
    """Raised when a probe statement cannot be parsed.

    Note that this is a fatal error.

    """


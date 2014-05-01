"""Parse single-nucleotide polymorphism events from human-readable statements.

"""
import re

from probe_generator import reference, sequence
from probe_generator.probe import InvalidStatement, NonFatalError

_SNP_REGEX = re.compile(r"""
        \s*
        ([a-zA-Z0-9.]+) # chromosome
        \s*
        :               # colon separator
        \s*
        (\d+)           # base pair index
        \s*
        ([acgtACGT])    # reference base
        \s*
        >               # arrow separator
        \s*
        ([acgtACGT])    # mutant base
        \s*
        /               # solidus separator
        \s*
        (\d+)           # bases
        \s*
        """, re.VERBOSE)

_SNP_STATEMENT_SKELETON = "{chromosome}:{index}_{reference}>{mutation}/{bases}"


class SnpProbe(object):
    """A probe for a single-nucleotide polymorphism event.

    The statement is in the following form:

        chromosome:index reference>mutant / bases

    For example: '1:100c>g/50' would be a statement for a 50-base pair probe
    with a mutation from a C to a G at the 100th base pair of the first
    chromosome.

    """
    def __init__(self, specification):
        self._spec = specification

    def __str__(self):
        return _SNP_STATEMENT_SKELETON.format(**self._spec)

    def sequence(self, genome):
        """Return the sequence of the probe.

        The mutant at the index (probe_length // 2).

        """
        start, end = _get_bases(self._spec)
        raw_bases = reference.bases(
                genome,
                self._spec["chromosome"],
                start,
                end)
        return self._mutate(raw_bases)

    @staticmethod
    def from_statement(statement):
        spec = _parse(statement)
        return SnpProbe(spec)

    def _mutate(self, bases):
        """Return the base pair sequence with the reference base of the spec
        replaced with the mutation base.

        Automatically reverse-complements the sequence if necessary.

        """
        mutation_index = (self._spec["bases"] // 2) - 1
        genome_ref_base = bases[mutation_index].lower()
        spec_ref_base = self._spec["reference"].lower()
        if genome_ref_base == spec_ref_base:
            return (bases[:mutation_index] +
                    self._spec["mutation"] +
                    bases[mutation_index+1:])
        elif genome_ref_base == sequence.complement(spec_ref_base):
            return (bases[:mutation_index]                      +
                    sequence.complement(self._spec["mutation"]) +
                    bases[mutation_index+1:])
        else:
            raise ReferenceMismatch(
                    "In probe {!s}: "
                    "Reference base {!r} does not match requested mutation "
                    "'{}>{}'".format(
                        self,
                        bases[mutation_index],
                        self._spec["reference"],
                        self._spec["mutation"]))


def _parse(statement):
    match = _SNP_REGEX.match(statement)

    if not match:
        raise InvalidStatement(
                "could not parse snp statement {!r}".format(
                    statement))

    chromosome, index, reference, mutation, bases = match.groups()
    return {"chromosome": chromosome,
            "index":      int(index),
            "reference":  reference,
            "mutation":   mutation,
            "bases":      int(bases)}


def _get_bases(spec):
    """Return the start and end indices from a SNP probe spec.

    """
    bases, index = spec["bases"], spec["index"]
    buffer = bases // 2
    return (index - buffer + 1), (index + buffer)


class ReferenceMismatch(NonFatalError):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """

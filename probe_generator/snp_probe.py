"""Parse single-nucleotide polymorphism events from human-readable statements.

"""
import re

from probe_generator import reference, sequence
from probe_generator.probe import InvalidStatement

_SNP_REGEX = re.compile(r"""
        \s*
        ([a-z0-9.]+) # chromosome
        \s*
        :            # colon separator
        \s*
        (\d+)        # base pair index
        \s*
        ([acgt])     # reference base
        \s*
        >            # arrow separator
        \s*
        ([acgt])     # mutant base
        \s*
        /            # solidus separator
        \s*
        (\d+)        # bases
        \s*
        """, re.VERBOSE | re.IGNORECASE)

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

        The mutant base is placed at the centre of the probe. If the number of
        bases in the probe is even, the location of the mutation in the probe
        is detemined by integer division (i.e., the mutation will be at the
        second base of a 4-base pair probe.

        """
        start, end = _get_bases(self._spec)
        raw_bases = reference.bases(
                genome,
                self._spec["chromosome"],
                start,
                end)
        return _mutate(raw_bases, self._spec)

    @staticmethod
    def from_statement(statement):
        spec = _parse(statement)
        return SnpProbe(spec)


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


def _mutate(bases, spec):
    """Return the base pair sequence with the reference base of the spec
    replaced with the mutation base.

    Automatically reverse-complements the sequence if necessary.

    Raises a ReferenceError if the reference base of the probe does
    not match the reference genome.

    """
    mutation_index = (spec["bases"] // 2) - 1
    genome_ref_base = bases[mutation_index].lower()
    spec_ref_base = spec["reference"].lower()
    if genome_ref_base == spec_ref_base:
        return (bases[:mutation_index] +
                spec["mutation"]       +
                bases[mutation_index+1:])
    elif genome_ref_base == sequence.complement(spec_ref_base):
        return (bases[:mutation_index]                +
                sequence.complement(spec["mutation"]) +
                bases[mutation_index+1:])
    else:
        raise ReferenceError(
                "Reference base {!r} does not match requested mutation "
                "'{}>{}'".format(
                    bases[mutation_index],
                    spec["reference"],
                    spec["mutation"]))


class ReferenceError(Exception):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """

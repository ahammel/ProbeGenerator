"""Parse single-nucleotide polymorphism events from human-readable statements.

"""
### TODO: make the rest of the *_statement modules look like this
import re

from probe_generator import reference

_SNP_REGEX = re.compile("""
        \s*
        ([a-z0-9.]+) # chromosome
        \s*
        :            # colon separator
        \s*
        (\d+)        # base pair index
        \s*
        ([acgt])     # reference base
        \s*
        >            # separator
        \s*
        ([acgt])     # mutant base
        \s*
        /            # solidus separator
        \s*
        (\d+)        # bases
        \s*
        """, re.VERBOSE | re.IGNORECASE)

_SNP_STATEMENT_SKELETON = "{chromosome}:{index}_{reference}>{mutation}/{bases}"

_COMPLEMENT = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        }


class SnpProbe(object):
    """A probe for a single-nucleotide polymorphism event.

    The statement is in the following form:

        chromosome:index reference>mutant / bases

    For example: '1:100c>g/50' would be a statement for a 50-base pair probe
    with a mutation from a C to a G at the 100th base pair of the first
    chromosome.

    """
    def __init__(self, statement, genome):
        self._spec = _parse(statement)
        self._genome = genome

    def __str__(self):
        return _SNP_STATEMENT_SKELETON.format(**self._spec)

    def sequence(self):
        """Return the sequence of the probe.

        The mutant base is placed at the centre of the probe. If the number of
        bases in the probe is even, the location of the mutation in the probe
        is detemined by integer division (i.e., the mutation will be at the
        second base of a 4-base pair probe.

        """
        start, end = _get_bases(self._spec)
        raw_bases = reference.bases(
                self._genome,
                self._spec["chromosome"],
                start,
                end)
        return _mutate(raw_bases, self._spec)


def _parse(statement):
    match = _SNP_REGEX.match(statement)

    if not match:
        raise InvalidFormat(
                "could not parse snp statement {!r}".format(
                    statement))

    chromosome, index, reference, mutation, bases = match.groups()
    return {"chromosome": chromosome,
            "index":      int(index),
            "reference":  reference,
            "mutation":   mutation,
            "bases":      int(bases)}


def _get_bases(spec):
    """Return the start and end indecies from a SNP probe spec.

    """
    bases, index = spec["bases"], spec["index"]
    buffer = bases // 2
    return (index - buffer + 1), (index + buffer)


def _mutate(bases, spec):
    """Return the base pair sequence with the reference base of the spec
    replaced with the mutation base.

    Automatically reverse-complements the sequence if necessary.

    """
    mutation_index = (spec["bases"] // 2) - 1
    if bases[mutation_index].lower() == spec["reference"].lower():
        return (bases[:mutation_index] +
                spec["mutation"]       +
                bases[mutation_index+1:])
    elif bases[mutation_index].lower() == _COMPLEMENT[spec["reference"]].lower():
        return (bases[:mutation_index]        +
                _COMPLEMENT[spec["mutation"]] +
                bases[mutation_index+1:])
    else:
        raise ReferenceError(
                "Reference base {!r} does not match requested mutation "
                "'{}>{}'".format(
                    bases[mutation_index],
                    spec["reference"],
                    spec["mutation"]))


class InvalidFormat(Exception):
    """Raised when a snp statement cannot be parsed.

    """


class ReferenceError(Exception):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """

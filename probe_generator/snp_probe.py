"""Parse single-nucleotide polymorphism events from human-readable statements.

"""
import re

from probe_generator import reference, sequence, annotation
from probe_generator.probe import InvalidStatement, NonFatalError

_STATEMENT_REGEX = re.compile(r"""
        \s*
        ([a-zA-Z0-9.]+) # chromosome
        \s*
        :               # colon separator
        \s*
        (\d+)           # base pair index
        \s*
        ([acgtACGT*])   # reference base
        \s*
        >               # arrow separator
        \s*
        ([acgtACGT*])   # mutant base
        \s*
        /               # solidus separator
        \s*
        (\d+)           # bases
        \s*
        (--.*|\s*)      # comment
        """, re.VERBOSE)

_STATEMENT_SKELETON = ("{chromosome}:{index}_"
                       "{reference}>{mutation}/{bases}{comment}")


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
        return _STATEMENT_SKELETON.format(**self._spec)

    def sequence(self, genome):
        """Return the sequence of the probe.

        The mutation is at the index (probe_length // 2).

        """
        start, end = annotation.get_bases(
            self._spec['bases'],
            self._spec['index'])
        raw_bases = reference.bases(
                genome,
                self._spec["chromosome"],
                start,
                end)
        return self._mutate(raw_bases)

    @staticmethod
    def explode(statement):
        """Yield probe statements with globbed reference and mutation
        bases filled in.

        """
        partial_spec = _parse(statement)
        specs = _expand(partial_spec)
        return [SnpProbe(spec) for spec in specs]

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
                    "Reference base {!r} does not match requested mutation "
                    "'{}>{}'".format(
                        bases[mutation_index],
                        self._spec["reference"],
                        self._spec["mutation"]))


def _parse(statement):
    """Return a partial SNP probe specification given a probe statement.

    """
    match = _STATEMENT_REGEX.match(statement)

    if not match:
        raise InvalidStatement(
            "could not parse snp statement {!r}".format(
                    statement))

    chromosome, index, reference_base, mutation, bases, comment = match.groups()
    return {"chromosome": chromosome,
            "index":      int(index),
            "reference":  reference_base,
            "mutation":   mutation,
            "bases":      int(bases),
            "comment":    comment}


def _expand(partial_spec):
    """Given a possibly globbed SNP probe specification, fill in all
    possible combinations of reference and mutation bases.

    """
    spec_reference = partial_spec['reference']
    spec_mutation = partial_spec['mutation']
    if spec_reference == spec_mutation == '*':
        ref_bases = "AC"
        # When both the reference and the mutant bases are globbed,
        # it's we only need one purine and one pyrimidine in the
        # reference in order to specify all six possible
        # mutations. The automatic reverse-complementing behaviour
        # takes care of T and G.
    elif spec_reference == '*':
        ref_bases = "ACGT"
    else:
        ref_bases = spec_reference

    mutant_bases = 'ACGT' if spec_mutation == '*' else spec_mutation

    for ref_base in ref_bases:
        for mutant_base in mutant_bases:
            if ref_base.upper() != mutant_base.upper():
                yield dict(partial_spec,
                           reference=ref_base,
                           mutation=mutant_base)


class ReferenceMismatch(NonFatalError):
    """Raised when the reference base of the genome does not match the
    reference base of the spec.

    """

"""Parse single-nucleotide polymorphism events from human-readable statements.

"""
import re

from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence import SequenceRange

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



class SnpProbe(AbstractProbe):
    """A probe for a single-nucleotide polymorphism event.

    The statement is in the following form:

        chromosome:index reference>mutant / bases

    For example: '1:100c>g/50' would be a statement for a 50-base pair probe
    with a mutation from a C to a G at the 100th base pair of the first
    chromosome.

    """
    _STATEMENT_SKELETON = ("{chromosome}:{index}_"
                           "{reference}>{mutation}/{bases}{comment}")

    def get_ranges(self):
        bases = self._spec['bases']
        chromosome = self._spec['chromosome']
        index = self._spec['index'] - 1 # Convert from 0- to 1-based indexing
        left_buffer = bases // 2 - 1
        right_buffer = bases - left_buffer
        return (
            SequenceRange(chromosome,
                          index-left_buffer,
                          index),
            SequenceRange(chromosome,
                          index,
                          index+1,
                          mutation=True),
            SequenceRange(chromosome,
                          index+1,
                          index+right_buffer))

    @classmethod
    def explode(cls, statement):
        """Yield probe statements with globbed reference and mutation
        bases filled in.

        """
        partial_spec = cls._parse(statement)
        specs = _expand(partial_spec)
        return [SnpProbe(spec) for spec in specs]


    @staticmethod
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

    ref_bases    = 'ACGT' if spec_reference == '*' else spec_reference
    mutant_bases = 'ACGT' if spec_mutation  == '*' else spec_mutation

    for ref_base in ref_bases:
        for mutant_base in mutant_bases:
            if ref_base.upper() != mutant_base.upper():
                yield dict(partial_spec,
                           reference=ref_base,
                           mutation=mutant_base)

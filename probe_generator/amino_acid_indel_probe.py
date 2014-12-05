"""Probe based on a blah blah blah FIXME.

"""
import re
import sys
import itertools

from probe_generator import annotation, transcript
from probe_generator.variant import TranscriptVariant, GenomeVariant
from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence import translate, reverse_translate

_STATEMENT_REGEX = re.compile(r"""
        \s*                                                     # whitespace
        ([A-Za-z0-9_./-]+)                                      # gene name
        \s*
        :
        \s*
        ([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*])           # IUPAC amino acid code
        \s*
        ([0-9]+)                                                # codon number
        \s*
        -
        \s*
        ([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*])
        \s*
        ([0-9]+)
        \s*
        (del|)
        \s*
        (ins\s*[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*]+|)
        \s*
        (\[trans\]|)                                            # transcript-only-sequence
        \s*
        /
        \s*
        ([0-9]+)                                                # number of base pairs
        \s*
        (--.*|\s*)                                              # comment
        """, re.VERBOSE)

class AminoAcidIndelProbe(AbstractProbe):
    """Probe for a single nucleotide polymorphism from an amino acid change at
    a codon.

    """
    _STATEMENT_SKELETON = ("{gene}:{codon}.{deletion}{insertion}"
                           "({deletion_nucleotides}{insertion_nucleotides})"
                           "{transcript_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{coordinate}{comment}")

    def __init__(self, *, variant, index, comment):
        self.variant = variant
        self.index = index
        self.comment = comment

    def __str__(self):
        reference_aa=translate(self.variant.reference)
        mutation_aa=translate(self.variant.mutation)
        return self._STATEMENT_SKELETON.format(
            gene=self.variant.gene,
            codon=self.index,
            insertion=_ins(mutation_aa),
            deletion=_del(reference_aa),
            insertion_nucleotides=_ins(self.variant.mutation),
            deletion_nucleotides=_del(self.variant.reference),
            transcript_sequence='[trans]' if self.variant.is_transcript else '',
            bases=len(self.variant),
            transcript_name=self.variant.transcript_name,
            chromosome=self.variant.chromosome,
            coordinate=self.variant.coordinate,
            comment=self.comment)

    def get_ranges(self):
        return self.variant.sequence_ranges()

    @staticmethod
    def explode(statement, genome_annotation=None):
        """Given a probe statement and a genome annotation, yield a sequence of
        AminoAcid probes matching the statement.

        If more than one probe has identical genomic coordinates, only the
        first is returned.

        """
        probes = []

        if genome_annotation is None:
            genome_annotation = []

        specification = _parse(statement)
        transcripts = annotation.lookup_gene(
            specification['gene'],
            genome_annotation)

        reference_peptide = _xify(
                specification["left_aa"],
                specification["right_aa"],
                (specification["right_codon"]-specification["left_codon"])-1)

        if specification["deletion"]:
            mutation_peptide = specification["insertion"]
        else:
            mutation_peptide = (
                    specification["left_aa"] +
                    specification["insertion"] +
                    specification["right_aa"])

        reference_sequences = reverse_translate(reference_peptide)
        mutation_sequences = reverse_translate(mutation_peptide)

        if specification["transcript_sequence"] == '':
            variant_class = GenomeVariant
        else:
            variant_class = TranscriptVariant

        index = specification["left_codon"]

        coordinate_cache = set()
        for txt, reference, mutation in itertools.product(
            transcripts, reference_sequences, mutation_sequences):
            try:
                sequence_range_index = txt.codon_index(
                        specification["left_codon"])
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                variant = variant_class(
                    transcript=txt,
                    index=sequence_range_index,
                    reference=reference,
                    mutation=mutation,
                    length=specification["bases"])
                if (variant.index, reference, mutation) not in coordinate_cache:
                    coordinate_cache.add((variant.index, reference, mutation))
                    probe = AminoAcidIndelProbe(
                        variant=variant,
                        index=index,
                        comment=specification["comment"])
                    probes.append(probe)
        return probes


def _parse(statement):
    """Return a partial AminoAcidProbe given a probe statement.

    Raises an InvalidStatement exception when the statement does not match the
    format of an amino acid probe.

    """
    match = _STATEMENT_REGEX.match(statement)
    if not match:
        raise InvalidStatement

    (gene,
     left_aa,
     left_codon,
     right_aa,
     right_codon,
     deletion,
     insertion,
     transcript_sequence,
     bases,
     comment) = match.groups()

    if deletion == '' and insertion == '':
        raise InvalidStatement

    return {"gene" :                gene,
            "left_aa" :             left_aa,
            "left_codon" :          int(left_codon),
            "right_aa" :            right_aa,
            "right_codon" :         int(right_codon),
            "deletion" :            bool(deletion),
            "insertion" :           re.sub(r"^ins +", "", insertion),
            "bases" :               int(bases),
            "transcript_sequence" : transcript_sequence,
            "comment" :             comment}


def _xify(head, tail, n):
    """Return the string generated by joining the head to the tail with n "x"s
    in between.

    Pronounced "ex-if-fye"

    """
    return head + ("X"*n) + tail


def _ins(string):
    return string if string == '' else "ins"+string

def _del(string):
    return string if string == '' else "del"+string

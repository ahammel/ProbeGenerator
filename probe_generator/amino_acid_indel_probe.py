"""Probe based on the insertion or deletion of a specified range of codons.

"""
import re
import sys
import itertools

from probe_generator import annotation, transcript
from probe_generator.variant import TranscriptVariant, GenomeVariant
from probe_generator.probe import AbstractProbe, InvalidStatement
from probe_generator.sequence import translate, reverse_translate
from probe_generator.sequence_range import SequenceRange

_STATEMENT_REGEX = re.compile(r"""
        \s*                                                     # whitespace
        ([A-Za-z0-9_./-]+)                                      # gene name
        \s*
        :
        \s*
        (del|)
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
        (ins\s*[ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy*]+|)
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
    _STATEMENT_SKELETON = ("{gene}:{deletion}"
                           "{left_aa}{left_index}-{right_aa}{right_index}"
                           "{insertion}{insertion_nucleotides}"
                           "{transcript_sequence}/{bases}_{transcript_name}_"
                           "{chromosome}:{coordinate}{comment}")

    def __init__(self, *, variant, index, comment):
        self.variant = variant
        self.index = index
        self.comment = comment

    def __str__(self):
        is_deletion = self.variant.index.start != self.variant.index.end
        left_index, right_index = self.index
        left_check, right_check = self.variant.check_sequences
        mutation_aa=translate(self.variant.mutation)
        return self._STATEMENT_SKELETON.format(
            gene=self.variant.gene,
            deletion='del' if is_deletion else '',
            left_aa=translate(left_check.reference),
            left_index=left_index,
            right_aa=translate(right_check.reference),
            right_index=right_index,
            insertion=_ins(mutation_aa),
            insertion_nucleotides=_bracket(self.variant.mutation),
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

        left_codon_seqs = reverse_translate(specification["left_aa"])
        right_codon_seqs = reverse_translate(specification["right_aa"])
        insertion_seqs = reverse_translate(specification["insertion"])

        if specification["transcript_sequence"] == '':
            variant_class = GenomeVariant
        else:
            variant_class = TranscriptVariant

        coordinate_cache = set()
        for txt, left_seq, right_seq, insertion in itertools.product(
            transcripts, left_codon_seqs, right_codon_seqs, insertion_seqs):
            try:
                left_index = txt.codon_index(
                        specification["left_codon"],
                        reference=left_seq,
                        mutation='')
                right_index = txt.codon_index(
                        specification["right_codon"],
                        reference=right_seq,
                        mutation='')
            except transcript.OutOfRange as error:
                print("{} in statement: {!r}".format(error, statement),
                      file=sys.stderr)
            else:
                if specification["deletion"] == '':
                    sequence_range_index = SequenceRange.between(
                            left_index, right_index)
                else:
                    sequence_range_index = SequenceRange.span(
                            left_index, right_index)
                variant = variant_class(
                    transcript=txt,
                    index=sequence_range_index,
                    mutation=insertion,
                    reference=None,
                    length=specification["bases"],
                    checks=[left_index, right_index]
                    )
                variant_spec = (variant.index, left_index, right_index, insertion)
                if variant_spec not in coordinate_cache:
                    coordinate_cache.add(variant_spec)
                    probe = AminoAcidIndelProbe(
                        variant=variant,
                        index=[specification["left_codon"],
                               specification["right_codon"]],
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
     deletion,
     left_aa,
     left_codon,
     right_aa,
     right_codon,
     insertion,
     transcript_sequence,
     bases,
     comment) = match.groups()

    if deletion == '' and insertion == '':
        raise InvalidStatement

    if deletion == '' and int(right_codon) - int(left_codon) != 1:
        raise InvalidStatement(
                "Amino acids must be adjacent when specifying an insertion")

    return {"gene" :                gene,
            "left_aa" :             left_aa,
            "left_codon" :          int(left_codon),
            "right_aa" :            right_aa,
            "right_codon" :         int(right_codon),
            "deletion" :            deletion,
            "insertion" :           re.sub(r"^ins +", "", insertion),
            "bases" :               int(bases),
            "transcript_sequence" : transcript_sequence,
            "comment" :             comment}


def _ins(string):
    return string if string == '' else "ins"+string
def _bracket(string):
    return string if string == '' else "({})".format(string)

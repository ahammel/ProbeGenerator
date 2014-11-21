"""Unified parser for all probe statements.

"""
from simpleparse.parser import Parser

DECLARATION = r"""

root             := statement

>statement<      := variant_stmnt / fusion_stmnt

>fusion_stmnt<   := coord_stmnt / exon_stmnt
>variant_stmnt<  := snp_stmnt / gene_snp_stmnt / amino_acid_stmnt / gene_indel_stmnt

coord_stmnt      := breakpoint, slash, breakpoint, comment?
exon_stmnt       := exon, w*, separator, w*, exon, comment?

snp_stmnt        := chromosome, colon, index, w*, nuc_mutation, slash, length, comment?
gene_snp_stmnt   := gene, colon, nuc_index, w*, nuc_mutation, trans?, slash, length, comment?
amino_acid_stmnt := gene, colon, aa_mutation, trans?, slash, length, comment?
gene_indel_stmnt := gene, colon, nuc_index, w*, indel, trans?, slash, length, comment?

breakpoint       := chromosome, colon, index, w*, direction, length
exon             := gene, exon_marker, opt_index, close_brkt, direction, w*, length

chromosome       := 'X' / 'Y' / ('2', [0-3]) / ('1', [0-9]) / [0-9]
gene             := [A-Za-z0-9_./-]+

nuc_mutation     := reference_nuc, short_arrow, mutation_nuc
aa_mutation      := reference_aa, w*, index, w*, mutation_aa

indel            := (deletion,w*,insertion) / deletion / insertion
deletion         := del_marker, w*, reference_seq
insertion        := ins_marker, w*, mutation_seq

reference_nuc    := [ACGTacgt]
mutation_nuc     := [ACGTacgt*]

reference_seq    := [ACGTacgt]+
mutation_seq     := [ACGTacgt]+

reference_aa     := [ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*]
mutation_aa      := [ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy*] # Includes 'X'

index            := [0-9]+
nuc_index        := index_marker, index
opt_index        := '*' / [0-9]+
length           := [0-9]+

direction        := [*+-]

separator        := '/' / '->'

comment          := w*, '--', non_newline*
trans            := w*,'[trans]',w*

<index_marker>   := w*, 'c.', w*

<exon_marker>    := w*, '#', w*, 'exon', w*, '[', w*
<close_brkt>     := w*, ']', w*

<del_marker>     := 'del.'
<ins_marker>     := 'ins.'

<colon>          := w*, ':', w*
<slash>          := w*, '/', w*
<short_arrow>    := w*, '>', w*

<non_newline>    := [ \t!-~]
<w>              := [ \t]

"""

STATEMENT_PARSER = Parser(DECLARATION)


def _flatten(parse_tree):
    for name, start, end, subtree in parse_tree:
        if subtree is not None and subtree != []:
            for s_name, s_start, s_end in _flatten(subtree):
                yield s_name, s_start, s_end
        else:
            yield name, start, end


def _extract_string(parse_element, string):
    _, start, end = parse_element
    return string[start:end]


def tokenize_statement(statement):
    success, parse_tree, _ = STATEMENT_PARSER.parse(statement)
    if success != 1:
        raise ParseError("{!r} could not be parsed".format(statement))
    return [(element[0], _extract_string(element, statement))
            for element in _flatten(parse_tree)]


class ParseError(ValueError):
    """Raised on invalid probe language statements.

    """

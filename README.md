probe-generator: make short sequences to probe for fusion events

`probe-generator` is a tool to make short sequences of base pairs representing
fusion events. These probes can be used to screen high-throughput sequencing
libraries for evidence of the events represented by the probe.

# Probe language

The input to `probe-generator` consists of statements in *probe language*.
Probe language is used to specify the genomic location of the fusion which will
be detected by the probe.

Statements in probe lanugage are in the form:

    "<gene>#<feature>[<number>] <side><bases> / <gene>#<feature>[<number>] <side><bases>"


    <gene>:    the name of the gene of interest. Acceptable characters are
               alphanumerics plus '_-/.' Case is not significant.

    <feature>: the name of the feature of interest ('exon', 'intron', etc.), or
               '*' Case is not significant.

    <number>:  the cardinality of the feature of intrest (1 for the first exon,
               etc.). Must be a digit or '*'.

    <side>:    Whether to return a sequence at the start or end of the feature.
               Acceptable characters are '+-*'.

    <bases>:   The length of probe sequence to return for this feature. Must be
               a digit or '*'.

 Any of "<feature>", "<number>", "<side>", or "<bases>" can be replaced with the
glob character ("*"), to indicate that any value is acceptable. In the "<bases>"
field, the interpretation is that the entire feature is desired.

Whitespace is not significant.

*EXONS ARE CURRENTLY THE ONLY FEATURE WHICH IS SUPPORTED*

The exact genomic location of the breakpoints of the fusion event can also be
specified directly using the coordinate-statement format:

    "chr:breakpoint(+|-)bases/chr:breakpoint(+|-)bases"

No globbing is allowed in coordinate statements. Whitespace is ignored.

## Ambiguity

Probe statements—unlike coordinate statements—do not necessarily specify a
unique location in the genome. When a probe statement could refer to any of
several genomic locations, `probe-generator` returns probes for every possible
location.

A probe statement can be ambiguous for three reasons:

    1. Globbing
    2. Alternative transcripts
    3. Non-unique gene names

Consider the following probe statement:

    FOO#exon[3] -25 / BAR#exon[*] +25

Imagine that FOO is alternatively spliced, so that there are two different
exons that could possibly be called the third. Furthermore, we will assume that
the symbol BAR identifies two different genes, each with two exons. In this
case, eight different probes will be generated:

    > FOO exon 3a / BARa exon 1
    > FOO exon 3a / BARa exon 2
    > FOO exon 3a / BARb exon 1
    > FOO exon 3a / BARb exon 2
    > FOO exon 3b / BARa exon 1
    > FOO exon 3b / BARa exon 2
    > FOO exon 3b / BARb exon 1
    > FOO exon 3b / BARb exon 2

where _3a_ and _3b_ are the two possible third exons of FOO and _BARa_ and
_BARb_ are the two genes called 'BAR'.

## Examples

To specifiy a probe covering the last 20 bases of the first exon of the gene
ABC and the first 30 bases of the third exon of DEF, you would pass the
following probe statement:

    "ABC#exon[1] -20 / DEF#exon[3] +30"

The same fusion, but with *any* exon of DEF:

    "ABC#exon[1] -20 / DEF#exon[*] +30"

Any fusion between introns of FOO and BAR with at least 40 bases covered:

    "FOO#intron[*] *20 / BAR#intron[*] *20"

Any fusion between any two features of SPAM and EGGS, with the entirity of both
features covered:

    "SPAM#*[*] ** / EGGS#*[*] **"

A probe for a fusion event between the 100th base pair of chromosome 1 and the
200th base pair of chromsome Y, with 25 bases on either side:

    "1:100-25/Y:200+25"

# Usage

    probe-generator --statement STMT  --genome GENOME --annotation FILE...
    probe-generator --coordinate COORD  --genome GENOME

    Options:
        -c COORD --coordinate=COORD     a file containing coordinate statements
        -s STMT --statement=STMT        a file containing fusion statements
        -g GENOME --genome=GENOME       the Ensembl reference genome
                                        (FASTA format)
        -a FILE --annotation=FILE       a genome annotation file in UCSC format


Currently, the RefSeq Genes and UCSC Genes annotation files are supported. More
than one annotation file can be specified for a single run:

    $ probe-generator -s statements.txt -g genome.fa \
                      -a refseq_genes.txt            \
                      -a ucsc_genes.txt

The resulting probes are printed to stdout in FASTA format. The titles of the
probes are the probe statements, followed by the unique identifiers of the rows
in the annotation file which were used, if applicable.

Annotations can be downloaded from [the UCSC table browser][ucsc_tables].


[ucsc_tables]: http://genome.ucsc.edu/cgi-bin/hgTables

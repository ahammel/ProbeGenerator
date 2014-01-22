probe-generator: make short sequences of virtual fusion events

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
               alphanumerics plus '_-/.'. Case is not significant.

    <feature>: the name of the feature of interest ('exon', 'intron', etc.), or
               '*'. Case is not significant.

    <number>:  the cardinality of the feature of interest (1 for the first exon,
               etc.). Must be a digit or '*'.

    <side>:    Whether to return a sequence at the start or end of the feature.
               Acceptable characters are '+-*'.

    <bases>:   The length of probe sequence to return for this feature. Must be
               a digit or '*'.

 Any of "<feature>", "<number>", "<side>", or "<bases>" can be replaced with the
glob character ("*"), to indicate that any value is acceptable. In the "<bases>"
field, the interpretation is that the entire feature is desired.

Whitespace is not significant.

## Examples

To specify a probe covering the last 20 bases of the first exon of the gene
ABC and the first 30 bases of the third exon of DEF, you would pass the
following probe statement:

    "ABC#exon[1] -20 / DEF#exon[3] +30"

The same fusion, but with *any* exon of DEF:

    "ABC#exon[1] -20 / DEF#exon[*] +30"

Any fusion between introns of FOO and BAR with exactly 40 bases covered:

    "FOO#intron[*] *20 / BAR#intron[*] *20"

Any fusion between any two features of SPAM and EGGS, with the entirety of both
features covered:

    "SPAM#*[*] ** / EGGS#*[*] **"

# Usage

[comment]: TODO

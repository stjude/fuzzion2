# fuzzion2

There are seven programs in the Fuzzion2 suite: `fuzzion2`, `fuzzort`,
`fuzzum`, `fuzzion2html`, `fuzzall`, `fuzzhop`, and `kmerank`.

The main program is `fuzzion2` which searches paired-end RNA or DNA to
identify read pairs that match fusion patterns.  The matching read pairs
are sorted by `fuzzort`; summarized by `fuzzum`; and converted to HTML by
`fuzzion2html` for viewing using a browser.  The `fuzzall` program
aggregates summaries produced by `fuzzum`.  The `fuzzhop` program reports
matching read pairs that may be artifacts due to index hopping.

The `fuzzion2` program requires as input a k-mer rank table.
Download the file named `fuzzion2_hg38_k15.krt` from
[`https://doi.org/10.5281/zenodo.6122447`](https://doi.org/10.5281/zenodo.6122447).
It holds a 4-GB 15-mer rank table that was constructed from the GRCh38 human
reference genome.  Use this file only when searching human RNA or DNA.
The `kmerank` program is provided to construct k-mer rank tables for other
species.  Depending on the size of the table under construction, you may
need to run `kmerank` with as much as 21 GB of memory.

*When `fuzzion2` reads a 4-GB k-mer rank table into memory, be sure to run
it with at least 5 GB of memory.*  Each of the other programs in the Fuzzion2
suite requires less than 1 GB of memory.

## Build

```
$ git clone https://github.com/stjude/fuzzion2.git
$ cd fuzzion2
$ make
```

This builds executable files for all seven programs of the Fuzzion2 suite and
puts them in `build/bin`.

#### Dependencies

These programs are written in C++ and compiled using [g++] version 6 or later.

[HTSlib] version 1.10.2 or later is used to read BAM files.
Note: Set `CPATH` and `LIBRARY_PATH` before running `make` and set
`LD_LIBRARY_PATH` before running `fuzzion2` using these commands:

```
$ HTSLIB=HTSlib-installation-directory
$ export CPATH=$CPATH:$HTSLIB/include
$ export LIBRARY_PATH=$LIBRARY_PATH:$HTSLIB/lib
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HTSLIB/lib
```

[gunzip] is used to decompress gzipped FASTQ data.

[g++]: https://gcc.gnu.org/
[HTSlib]: https://github.com/samtools/htslib
[gunzip]: https://www.gnu.org/software/gzip/

## Usage

```
Usage: fuzzion2 OPTION ... [filename ...] > hits

These options are required:
  -pattern=filename   name of pattern input file
  -rank=filename      name of binary  input file containing the k-mer rank table

Specify -fastq1 and -fastq2, or -ifastq or -ubam, or list filenames on command line
  -fastq1=filename    name of FASTQ Read 1 input file
  -fastq2=filename    name of FASTQ Read 2 input file
  -ifastq=filename    name of interleaved FASTQ input file (may be /dev/stdin)
  -ubam=filename      name of unaligned Bam input file

The following are optional:
   N is a numeric value, e.g., -threads=4
  -maxins=N     maximum insert size in bases. . . . . . . . . . . . default 500
  -maxrank=N    maximum rank percentile of minimizers . . . . . . . default 99.9
  -maxtrim=N    maximum bases second read aligned ahead of first. . default 5
  -minbases=N   minimum percentile of matching bases. . . . . . . . default 90.0
  -minmins=N    minimum number of matching minimizers . . . . . . . default 1
  -minov=N      minimum overlap in number of bases. . . . . . . . . default 5
  -show=N       show best only (1) or all patterns (0) that match . default 1
  -single=N     show single-read (1) or just read-pair (0) matches. default 0
  -threads=N    number of threads . . . . . . . . . . . . . . . . . default 8
  -w=N          window length in number of bases. . . . . . . . . . default 10
```

The `-pattern` option is required and specifies the name of a text file containing
RNA or DNA patterns.  Look in the `patterns` directory for pattern files provided
with this distribution.

The `-rank` option is also required and specifies the name of a binary file
containing a k-mer rank table.  See above for where to download this file.

`fuzzion2` expects read pairs as input; single-read formats are unsupported.
Each read pair is examined to see if it matches any of the patterns in the pattern
file.  Read pairs are obtained from:

1. a pair of FASTQ files identified by the `-fastq1` and `-fastq2` options; or
1. an interleaved FASTQ file named by the `-ifastq` option; or
1. an unaligned Bam file named by the `-ubam` option; or
1. one or more files listed on the command line.

`fuzzion2` expects that the mates of a read pair are adjacent in interleaved FASTQ
and unaligned Bam files.  If a pair of FASTQ files is specified, the mates are
separated; "Read 1" mates are in one file and corresponding "Read 2" mates are
in another file.  If the name of a FASTQ file ends with ".gz", `fuzzion2` assumes
that the file is gzipped and uses `gunzip` to decompress it.

If `-fastq1`, `-fastq2`, `-ifastq`, and `-ubam` options are omitted, `fuzzion2`
looks for file names on the command line.  The named files can be any combination
of the above file types and can be listed in any order.  `fuzzion2` automatically
recognizes and pairs up corresponding "Read 1" and "Read 2" FASTQ files.

Each read pair that matches a pattern is a "hit" and is written along with the
pattern on three lines to the standard output stream.  In the following example,
`BCR-ABL1` is the name of the pattern.  The second column of the first line shows
the substring of the pattern sequence that matches the read pair.  The second and
third lines show the read name and entire sequence of each mate.

```
pattern BCR-ABL1                 CATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAATAA]GGAAGA[AGCCCTTCAGCGGCCAGTAGCATCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGC
read EXAMPLE:1105:12909:66982/2  CATCCGTGGAGCTGCAGATGCTGACCAACTCGTGTGTGAAACTCCAGACTGTCCACAGCATTCCGCTGACCATCAATAAGGAAGAAGCCCTTCAGCGGCCA
read EXAMPLE:1105:12909:66982/1                                                                                                               TCTGACTTTGAGCCTCAGGGTCTGAGTGAAGCCGCTCGTTGGAACTCCAAGGAAAACCTTCTCGCTGGACCCAGTGAAAATGACCCCAACCTTTTCGTTGC
```

A final line is written to the standard output stream showing the total number
of read pairs processed by the program.

#### Additional Options

You will normally not need to specify any of the options described below, but we mention
them here for completeness.

If there are similar patterns, a read pair may match more than one of them.  By default,
only the best match is reported.  Specify `-show=0` to see all of them.

By default, both mates of a read pair must match a pattern.  If there are pattern
sequences too short to match both mates, set `-single=1` to see also matches of one mate
to patterns.

If there are patterns consisting of very common k-mers, some matches may be missed.
In this case, the value of the `-maxrank` option should be increased.  By default, this
option is set to 99.9, which means the 0.1% most common k-mers in the reference genome are
ignored.  Increasing the value to 100.0 processes all k-mers, but the program may run
slowly.

The `-w` option specifies the window length; reducing this value increases the number
of minimizers representing each sequence.  The `-minmins` option specifies the minimum
number of minimizers shared by a read sequence and pattern sequence in a candidate match.

An alignment of a read pair to a pattern sequence is reported as a hit only if the
alignment has a reasonable insert size (not greater than the value of the `-maxins`
option); the read pair overlaps each side of the pattern by at least the value of the
`-minov` option; the percentage of bases in agreement on each side must be at least the
value of the `-minbases` option; and if the second read of the pair aligns ahead of
the first read, it is by no more than the value of the `-maxtrim` option.  Normally,
the first read will align ahead of the second read, but the latter option is provided
to accommodate imprecise adapter trimming.

#### Postprocessing

When multithreading is used (i.e., the value of the `-threads` option is greater
than 1), the order of the hits is indeterminate and a simple `diff` cannot be
used to compare output files.  It is therefore recommended to run the `fuzzort`
program to sort the hits.  This program sorts the hits by pattern name so that
the hits are in a determinate order (and `diff` can be used to compare files)
and the hits of a pattern are together.

```
Usage: fuzzort < fuzzion2_hits > sorted_hits
```

Run `fuzzion2html` to produce an HTML file that provides an attractive display
of hits when opened in a browser such as Google Chrome or Microsoft Edge.
SNPs, indels, and sequencing errors are highlighted in the display.  An example
can be seen [here](https://htmlpreview.github.io/?https://github.com/stjude/fuzzion2/blob/master/test/SJBALL020765_D1_by_pattern.html).

```
Usage: fuzzion2html OPTION ... < fuzzion2_hits > html

The following are optional:
  -group=string   comma-separated list of column headings, default is no grouping
  -strong=N       minimum overlap of a strong match in #bases, default is 15
  -title=string   string to include in the title of the HTML page
```

The input to `fuzzort`, `fuzzion2html`, and `fuzzum` may be the output from a
single run of `fuzzion2` or the concatenation of outputs from multiple runs.

`fuzzum` produces a tab-delimited summary of hits, and these summaries may be
aggregated using the `fuzzall` program.

```
Usage: fuzzum OPTION ... < fuzzion2_hits > hit_summary

This option is required:
  -id=string      identifies the sample

The following are optional:
  -group=string   comma-separated list of column headings, default is no grouping
  -strong=N       minimum overlap of a strong match in #bases, default is 15

Usage: fuzzall OPTION fuzzum_filename ... > pattern_summary

The following is optional:
  -dataset=name   name associated with this dataset
```

These summaries indicate the number of distinct read pairs matching each pattern,
and of those the number of "strong" versus "weak" matches.  A match is considered
to be strong if the alignment of the read pair to the pattern overlaps each side of
the pattern by at least N bases, where the value of N is given by the `-strong` option.
Otherwise, a match is regarded as weak due to insufficient overlap.  The "strong" matches
are divided into two types: "strong+" if at least one read is junction spanning, and
"strong-" if neither read spans the junction.

In the leftmost column of a display of hits produced by `fuzzion2html`, each hit is
labeled as either "weak," "strong-", or "strong+", or as "dup" if the read pair is a
duplicate of another hit.  Furthermore, each read is marked as "+" if it qualifies as
junction spanning and "-" if it does not.  As used here, "+" and "-" do not indicate
orientation.

In `fuzzall` output, each sample ID is followed by two numbers in parentheses,
e.g., (24/22), indicating the number of distinct matches (24) and "strong+" matches (22).

It is possible to summarize the hits by pattern "group."  The `-group` option to
`fuzzion2html` and `fuzzum` specifies a comma-separated list of annotation column
headings in the pattern file.  The first heading in the list identifies the column by
which hits will be grouped in the summaries; hits of patterns having the same value
in this column are grouped together.  The other headings in the list identify "group"
annotation columns.

It is possible that a read pair matching a pattern was assigned to the wrong sample
during the sequencing process.  This phenomenon is known as "index hopping."  Given the
hits from two or more samples that were sequenced together, the `fuzzhop` program reports
the numbers of hits of a pattern that came from the same flow cell and lane but were
assigned to different samples.  Each input file contains the hits from a single sample.

```
Usage: fuzzhop fuzzion2_filename1 fuzzion2_filename2 ... > possible_index_hops
``` 

#### Example Run

The `test` directory contains some files you can use to run a simple test:

```
fuzzion2 -pattern=example_patterns.txt -rank=fuzzion2_hg38_k15.krt \
   -fastq1=example_input1.fq -fastq2=example_input2.fq > my_output.txt

fuzzort < my_output.txt > my_sorted_output.txt

fuzzion2html -title="Fuzzion2 Example" < my_output.txt > my_output.htm

fuzzum -id=example < my_output.txt > my_output_summary.txt
```

## Patterns

Pattern files describe the nucleotide-level breakpoints of sequences of
interest, e.g., fusions, ITD (internal tandem duplication) boundaries,
or other targets.

#### File Format

A pattern file is formatted as tab-delimited text with a header line.  Two columns
must be present in the file:

* "pattern" contains the pattern identifier.  In our pattern set
  this is the gene pairing with a numbered suffix, e.g., BCR-ABL1-01 (note
  that these identifiers are not stable between releases).
* "sequence" contains the sequence spanning the breakpoint.  A single
  pair of brackets is used in each sequence to indicate the boundaries
  of the breakpoint.  Two types of brackets may be
  used, square brackets ("]" and "[") for fusions, and curly brackets ("}" and "{")
  for ITD boundaries.  The right or closing bracket appears first, marking the end
  of the sequence upstream of the breakpoint (e.g., the first gene fusion partner),
  followed by the left or opening bracket indicating the start of the sequence
  downstream of the breakpoint (e.g., the second fusion gene partner). Sequence may
  optionally appear between brackets, indicating either interstitial sequence or a
  region of microhomology (i.e., a portion of the sequence that is ambiguous between
  the two sides of the breakpoint).  If curly brackets are used, `fuzzion2` will
  require at least one read to span the breakpoint.  Flanking sequence of 400-500 nt
  on either side of the breakpoint is recommended, or whatever length is appropriate
  for your sequencing's insert size.

Additional columns may also be added to the pattern file for any other desired
information or annotations.  Below is an example pattern sequence for a BCR-ABL1
fusion:

```
AGGGCGCCTTCCATGGAGACGCAGA][AGCCCTTCAGCGGCCAGTAGCATCT
```

This is a just a very short excerpt of the pattern sequence around the
breakpoint for illustrative purposes.  Square brackets appear in the
pattern, indicating a fusion event.  The sequence to the left of the
"]" is from BCR, the sequence to the right of the "[" is from ABL1.

#### Sources

Pattern sets distributed with Fuzzion2 can be found in the "patterns"
subdirectory of this repository.  These sets were generated from fusion and ITD
data from various pediatric cancer projects and collaborations at St. Jude, such
as PCGP and NCI TARGET, as well as from clinical sequencing.  Patterns were also
generated from fusions described in the COSMIC database.  Pattern sets are works
in progress.

#### Programs for Creating Pattern Files

Programs to generate pattern files from either fusion/ITD contig sequences or genomic
breakpoints are available here: [`https://github.com/stjude/fuzzion2_patgen/`](https://github.com/stjude/fuzzion2_patgen/).

## COPYRIGHT

Copyright 2025 St. Jude Children's Research Hospital

## LICENSE

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

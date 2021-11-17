# fuzzion2

There are six programs in the Fuzzion2 suite: `fuzzion2`, `fuzzort`,
`fuzzum`, `fuzzion2html`, `fuzzall`, and `kmerank`.

The main program is `fuzzion2` which searches paired-end RNA or DNA to
identify read pairs that match fusion patterns.  The matching read pairs
are sorted by `fuzzort`; summarized by `fuzzum`; and converted to HTML by
`fuzzion2html` for viewing using a browser.  The `fuzzall` program
aggregates summaries produced by `fuzzum`.

The `fuzzion2` program requires as input a k-mer rank table.  Download
this file from `http://ftp.stjude.org/pub/software/fuzzion2_hg38_k15.krt`.
It holds a 4-GB 15-mer rank table that was constructed from the GRCh38 human
reference genome.  Use this file only when searching human RNA or DNA.
The `kmerank` program is provided to construct k-mer rank tables for other
species.  Depending on the size of the table under construction, you may
need to run `kmerank` with as much as 21 GB of memory.

*When `fuzzion2` reads a 4-GB k-mer rank table into memory, be sure to run
it with at least 5 GB of memory.*  Each of the other programs in the Fuzzion2
suite requires less than 1 GB of memory.

## Quick Start

### Developer's Build

```
$ git clone https://github.com/stjude/fuzzion2.git
$ cd fuzzion2
$ make
```

This builds executable files for all six programs of the Fuzzion2 suite and
puts them in `build/bin`.

#### Dependencies

These programs are written in C++ and compiled using [g++] version 6 or later.

[HTSlib] version 1.10.2 or later is used to read BAM files.
Set `CPATH` and `LIBRARY_PATH` before running `make` and set
`LD_LIBRARY_PATH` before running `fuzzion2`.

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
Usage: fuzzion2 OPTION ... [ubam_filename ...] > hits

These options are required:
  -pattern=filename   name of pattern input file
  -rank=filename      name of binary  input file containing the k-mer rank table

The following are optional:
  -fastq1=filename    name of FASTQ Read 1 input file
  -fastq2=filename    name of FASTQ Read 2 input file
  -ifastq=filename    name of interleaved FASTQ input file (may be /dev/stdin)

The following are optional:
   N is a numeric value, e.g., -threads=4
  -maxins=N     maximum insert size in bases. . . . . . . . . . . . default 500
  -maxrank=N    maximum rank percentile of minimizers . . . . . . . default 95.0
  -maxtrim=N    maximum bases second read aligned ahead of first. . default 5
  -minbases=N   minimum percentile of matching bases. . . . . . . . default 90.0
  -minmins=N    minimum number of matching minimizers . . . . . . . default 3
  -minov=N      minimum overlap in number of bases. . . . . . . . . default 5
  -show=N       show best only (1) or all patterns (0) that match . default 1
  -single=N     show single-read (1) or just read-pair (0) matches. default 0
  -threads=N    number of threads . . . . . . . . . . . . . . . . . default 8
  -w=N          window length in number of bases. . . . . . . . . . default 5
```

The `-pattern` option is required and specifies the name of a text file containing
RNA or DNA patterns.  Look in the `patterns` directory for pattern files provided
with this distribution.

The `-rank` option is also required and specifies the name of a binary file
containing a k-mer rank table.  See the note above.

`fuzzion2` examines each read pair to see if it matches any of the patterns in the
pattern file.  Read pairs are obtained from:

1. a single interleaved FASTQ file named by the `-ifastq` option;
1. a pair of FASTQ files identified by the `-fastq1` and `-fastq2` options; or
1. one or more unaligned Bam files specified on the command line.

`fuzzion2` expects that the mates of a read pair are adjacent in interleaved FASTQ
files and unaligned Bam files.  If a pair of FASTQ files is specified, the mates
are separated; "Read 1" mates are in the first file and corresponding "Read 2"
mates are in the second file.  If the name of a FASTQ file ends with ".gz",
`fuzzion2` assumes that the file is gzipped and uses `gunzip` to decompress it.

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

When multithreading is used (i.e., the value of the `-threads` option is greater
than 1), the order of the hits is indeterminate and a simple `diff` cannot be
used to compare output files.  It is therefore recommended to run the `fuzzort`
program to sort the hits.  This program sorts the hits by pattern name so that
the hits are in a determinate order (and `diff` can be used to compare files)
and the hits of a pattern are grouped together.

```
Usage: fuzzort < fuzzion2_hits > sorted_hits
```

If the read pairs for a sample are stored in multiple pairs of FASTQ files, run
`fuzzion2` once for each pair of FASTQ files and concatenate the resulting output
files.  Then run `fuzzort` on the concatenated file to get a single sorted file
of hits.

Run `fuzzion2html` to produce an HTML file that provides an attractive display
of hits when opened in a browser such as Google Chrome or Microsoft Edge.
SNPs, indels, and sequencing errors are highlighted in the display.

```
Usage: fuzzion2html OPTION ... < fuzzion2_hits > html

The following is optional:
  -title=string   string to include in the title of the HTML page
```

`fuzzum` produces a tab-delimited summary of hits, and these summaries may be
aggregated using the `fuzzall` program.

```
Usage: fuzzum OPTION ... < fuzzion2_hits > hit_summary

This option is required:
  -id=string   identifies the sample

Usage: fuzzall fuzzum_filename ... > pattern_summary
```

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

#### File format

A pattern file is formatted as tab-delimited text with a header line.  Two columns must be present in the file:

* "pattern" contains the pattern identifier.  In our pattern set
  this is the gene pairing with a numbered suffix, e.g. BCR-ABL1-01 (note
  that these identifiers are not stable between releases).
* "sequence" contains the sequence spanning the breakpoint.  A single
  pair of brackets is used in each sequence to indicate the boundaries
  of the breakpoint.  Two types of brackets may be
  used, square brackets ("]" and "[") for fusions, and curly brackets ("}" and "{") for ITD
  boundaries.  The right or closing bracket appears first, marking the end of the sequence upstream of the breakpoint (e.g. the first gene fusion partner), followed by the left or opening bracket indicating the start of the sequence downstream of the breakpoint (e.g. the second fusion gene partner). Sequence may optionally appear between brackets,
  indicating either interstitial sequence or a region of microhomology 
  (i.e., a portion of the sequence that is ambiguous between the two 
  sides of the breakpoint).  If curly brackets are used, fuzzion2 will
  require at least one read to span the breakpoint.  Flanking sequence
  of 400-500 nt on either side of the breakpoint is recommended, or
  whatever length is appropriate for your sequencing's insert size.

Additional columns may also be added to the pattern file for any other desired information or annotations.  Below is an example pattern sequence for a BCR-ABL1 fusion:

```
AGGGCGCCTTCCATGGAGACGCAGA][AGCCCTTCAGCGGCCAGTAGCATCT
```

This is a just a very short excerpt of the pattern sequence around the
breakpoint for illustrative purposes.  Square brackets appear in the
pattern, indicating a fusion event.  The sequence to the left of the
"]" is from BCR, the sequence to the right of the "[" is from ABL1.


#### Sources

The pattern set distributed with fuzzion2 can be found in the "patterns"
subdirectory of this repo.  This set was generated from fusion and ITD
data from various pediatric cancer projects and collaborations
at St. Jude, such as PCGP and NCI TARGET, as well as from clinical
sequencing.  Patterns were also generated from fusions described in the COSMIC database.  The pattern set is a work in progress.

#### Programs for creating pattern files

Coming soon.

## COPYRIGHT

Copyright 2021 St. Jude Children's Research Hospital

## LICENSE

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# assessPrimers
A tool to aid PCR primer design and evaluation. Developed in association with the [In Situ Lab](https://insitulabs.org/) initiative.

The primer assessment tool requires several inputs including: 1) list of nucleotide reference sequences, 2) list of forward and reverse primers, respectively (degenerate bases are allowed), 3) a reference protein sequence (optional). The inputs are used to create a non-redundant multiple sequence alignment of all references sequences to each other as well as to each primer pair. From this alignment, the following statistics are printed to stdout:

Primer number: The number of unique primers calculated after converting all degenerate bases to their non-degenerate equivalents.

Entropy: Cumulative entropy score for each length of k-mer along the alignment (lower entropy scores reflect more conserved sequences)

Start coordinate position

Number mismatches: For each primer provided, histogram of the number of mismatches for each reference sequence

Additional outputs

Sequence statistics text file (tab-delimited)
Provides the number of mismatches detected for each sequence primer combination

K-mer statistics text file (tab-delimited)
For each k-mer along the alignment: coordinate start position, cumulative entropy score, number of gaps, number of non-degenerate primers to match all reference sequences, degenerate consensus sequence, primer name (default name is ‘blank’, unless the k-mer matches one of the input primers)

Multiple sequence alignment files


Notes: 1) temporary and in-between outputs are saved to the same location as the final outputs. 2) Any number of forward and reverse primers can be analyzed, and the current version of this software requires that they be paired.

Stepwise Analysis pipeline:
Reverse-complement reverse primers
Filter reference sequences with N content > filtN parameter (default = 0.05)
Make blast database from reference sequence file
If reference protein sequence has been provided remove references that do not match
Collapse reference sequences with identity > idcutoff parameter (default = 0.9).
Create multiple sequence alignment file
Generate primer and sequence statistics


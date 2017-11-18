## indelope: find indels and SVs too small for structural variant callers and too large for GATK

`indelope` was started with the goal of increasing the diagnostic rate in exomes. To do this it must be:

+ fast : ~2.5 CPU-minutes per exome (25% slower than `samtools view -c`)
+ easy-to-use : goes from BAM to VCF in a single command.
+ novel : it does local assembly and then [aligns](https://github.com/lh3/ksw2) assembled contigs
  to the genome to determine the event, and then does k-mer counting (not alignment) to genotype
  without k-mer tables.
+ accurate : because of the genotyping method, we know that called variants are not present
  in the reference.

These features will help ensure that it is actually used (fast, easy-to-use) and that it finds new
and valid variation.

As of November 2017, `indelope` is working -- it finds large indels that are clearly valid by visual
inspection that are missed by GATK/freebayes/lumpy. I am in the process of evaluating whether this
truly increases the diagnostic rate and determining ways to increase the sensitivity.

`indelope` also works on whole genomes, but, for now, that is not the target use-case.

## how it works

`indelope` sweeps over a single bam and finds regions that are likely to harbor indels--reads that have
more than 1 cigar event and split-reads (work on split reads is in progress). As it finds these it increments
a counter for the genomic position of the event. Upon finding a gap in coverage, it goes back, finds any
previous position with sufficient evidence (this is a parameter) of an event, gathers reads that have been
aligned across that genomic position (and unaligned reads from that region) and does assembly on those reads.
It then aligns the assembled contigs to the genome using [ksw2](https://github.com/lh3/ksw2) and uses the CIGAR
to determine the event as it's represented in the VCF. Any event will result in a novel k-mer not present in
the reference genome; `indelope` gets the k-mer of the reference genome at the event and the novel k-mer of
the alternate event. It again iterates through the reads that were originally aligned to the event and counts
reference and alternate k-mers. Those counts are used for genotyping. Note that this reduces reference bias
because we are aligning a contig (often >400 bases) to the genome and never re-aligning the actual reads.

As `indelope` sweeps across the genome, it keeps the reads for each chunk in memory. A chunk bound is defined
by a gap in coverage; this occurs frequently enough that the memory use is negligible. Once a new chunk is reached,
all events from the previous chunk are called and then those reads are discarded. This method, along with the
assembly method make `indelope` extremely fast--given 2 BAM decompression threads, it can call variants in an
exome in ~ 1 minute (2.5 CPU-minutes).

## assembly

A read (contig) slides along another read (contig) to find the offset with the most matches. At each offset, if
more than $n mismatches are found, the next offset is attempted. This is efficient enough that a random read to
a random (non-matching) contig of length $N will incur ~ 1.25 * $N equality (char vs. char) tests.

Within each contig `indelope` tracks the number of reads supporting each base. Given a sufficient number of
reads supporting a contig, it can account for sequencing errors with a simple voting scheme. That is: if contig a,
position x has more than 7 supporting reads and contig b has fewer than 3 supporting reads (and we know that
otherwise, `a` and `b` have no mismatches), we can vote to set the mismatch in `b` to the apparent match in `a`.
This allows us to first create contigs allowing no mismatches within the reads and then to combine and extend contigs
using this voting method.

## installation and usage

get a binary from [here](https://github.com/brentp/indelope/releases)
and make sure that libhts.so is in your `LD_LIBRARY_PATH`

then run `./indelope -h` for usage. recommended is:

```
indelope --min-event-len 5 --min-reads 5 $fasta $bam > $vcf
```

## to do

+ somatic mode / filter mode. allow filtering on a set of k-mers from a parental genome (parent for 
  mendelian context or normal for somatic variation).

+ use SA tag. (and possibly discordant reads)

## see also

+ [svaba](https://github.com/walaj/svaba) does local assembly, but then genotypes by alignment to those
  assemblies. It is slower than `indelope` but it is an extremely useful tool and has a series of
  careful and insightful analyses in its paper. (highly recommend!!)

+ [rufus](https://github.com/jandrewrfarrell/RUFUS) does k-mer based variant detection; Andrew described
  to me the RUFUS assembly method that inspired the one now used in `indelope`.

+ [lancet](https://github.com/nygenome/lancet), [scalpel](http://scalpel.sourceforge.net/),
  [mate-clever](https://academic.oup.com/bioinformatics/article/29/24/3143/194997),  and [prosic2](https://github.com/prosic/prosic2) are all
  great tools that are similar in spirit that are worth checking out (the latter 2 are focused on somatic variation).


## notes and TODO

# need a better way to combine contigs

sometimes, can have 2 contigs, each of length ~ 80 and they overlap for 60 bases but cutoff is
e.g. 65. Need a way to recover this as it happens a lot in low-coverage scenarios. maybe it can
first combine, then trim (currently, it's trim, combine).
This should also allow more permissive overlaps if the correction list is empty.

# contigs

`min_overlap` in contig::best_match should be a float between 0 and 1 that will make sure that at least
that portion of the shortest contig overlaps the other.

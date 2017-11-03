+ iterate through bam and note locations of
  - discordant (by ISIZE)
  - split reads (SA tag)
  - > 1 cigar op
+ extract all primary reads from union of locations into
  fastq

+ assemble fastq with fermil into contigs
+ align contig to genomic position to get coordinates
+ align all reads from those coords to both contig and reference.
  - get ref, alt counts.

+ genotype based on ref and alt
+ for non hom-ref translate back to ref coords (TODO)


## TODO:
1. first merge reads looking at variant quality.
2. merge and vote but keep separate anything where support for either is >= 5
3. TODO: first trim reads before putting them into contig?
4. expand region by keeping normal reads at bounds.

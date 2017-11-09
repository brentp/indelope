## indelope: find medium sized indels with local realignment

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
0. set contig default to not insert if corrections is not empty. then in contigs.merge
   that can be used.

1. first merge reads looking at variant quality.
2. merge and vote but keep separate anything where support for either is >= 5

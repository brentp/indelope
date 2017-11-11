import osproc
import contig
import sequtils
import os
import times
import random
import math
import algorithm
import strutils
import docopt
import hts
import ./ksw2/ksw2
import kmer
import genotyper

type
  event = tuple[start: int, stop: int, len: int]
  roi = tuple[start: int, stop: int, reads:seq[Record]]

proc trim(sequence: var string, base_qualities: seq[uint8], min_quality:int=12) =
  var a = 0
  while a < base_qualities.high and base_qualities[a] < uint8(min_quality):
    a += 1

  if a == base_qualities.high:
    sequence.set_len(0)
    return

  var b = base_qualities.high
  while b > a and base_qualities[b] < uint8(min_quality):
    b -= 1

  if a != 0 or b != base_qualities.high:
    sequence = sequence[a..b]

proc skippable(r: Record, allow_unmapped:bool=false): bool {.inline.} =
  if r.chrom == "hs37d5": return true
  if r.chrom.startswith("GL"): return true
  var f = r.flag
  if f.dup or f.qcfail: return true
  if not allow_unmapped and f.unmapped: return true
  if f.supplementary or f.secondary: return true
  return false

type
  Variant* = ref object of RootObj
    chrom*: string
    start*: int
    reference*: string
    alternate*: string
    genotype*: Genotype
    ref_kmer*: string
    alt_kmer*: string
    AD*: array[2, int]

proc info(v:Variant): string =
   return ("DP=" & $(v.AD[0] + v.AD[1]) &
           ";AD=" & $v.AD[0] & "," & $v.AD[1] &
           ";ref_kmer=" & $v.ref_kmer &
           ";alt_kmer=" & $v.alt_kmer)

const header = """##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=DP,Number=1,Type=Integer,Description="supporting k-mer depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="supporting k-mer depth">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=ref_kmer,Number=1,Type=String,Description="reference kmer used for genotyping">
##INFO=<ID=alt_kmer,Number=1,Type=String,Description="alternate kmer used for genotyping">
$1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$2"""

proc `$`*(v:Variant): string =
  return format("$chrom\t$pos\t$id\t$ref\t$alt\t$qual\tPASS\t$info\tGT:GQ:GL\t$gt" %
                ["chrom", v.chrom, "pos", $v.start, "id", ".",
                 "ref", v.reference, "alt", v.alternate,
                 "qual", formatFloat(v.genotype.qual, precision=2, format=ffDecimal), "info", v.info(), "gt", $(v.genotype)])

proc same(a:Variant, b:Variant): bool =
  if a == nil or b == nil: return false
  return (a.start == b.start) and (a.chrom == b.chrom) and (a.reference == b.reference) and (a.alternate == b.alternate)

iterator callsemble(r: roi, fai: Fai, ez: Ez, min_ctg_len:int=74, min_reads:int=4, min_event_len:int=4, K:int=29): Variant =

  var contigs = new_seq[Contig]()
  var read_seq = ""
  # assemble the alt contig
  for read in r.reads:
    if read.skippable(allow_unmapped=true): continue
    contigs.insert(read.sequence(read_seq), read.start)
  var ucontigs = new_seq_of_cap[Contig](contigs.len)

  for c in contigs:
    if c.nreads < min_reads or c.len < min_ctg_len: continue
    ucontigs.add(c)
  contigs = ucontigs

  var sequence = ""
  var ss = ""
  var chrom = r.reads[0].chrom

  for ctg in contigs:

    var max_stop = ctg.start
    for read in r.reads:
      max_stop = max(max_stop, read.stop)

    var width = int((K + 1) / 2 - 1)
    var reference = fai.get(chrom, ctg.start, max_stop + width)
    ctg.sequence.align_to(reference, ez)

    var qlocs = toSeq(ez.query_locations())
    var ii = 0

    for loc in ez.target_locations(ctg.start):
      if loc.len.int < min_event_len: continue
      var rstart = max(0, loc.start - ctg.start - width)


      # BUG: this can fail if we are at end of reference e.g. MT
      var ref_kmer = reference[rstart..<(rstart + K)]
      var qloc = qlocs[ii]
      ii += 1

      var qstart = max(qloc.start - width, 0)
      if qstart + K >= ctg.len:
        qstart = (ctg.len - K - 1)
      var alt_kmer = ctg.sequence[qstart..<(qstart + K)]

      if len(ref_kmer) != len(alt_kmer) or len(ref_kmer) != K:
        stderr.write_line("bug!!! kmer lengths should be equal, ref:" & $ref_kmer.len & " alt:" & $alt_kmer.len)
        stderr.write_line("bug!!! qloc start:" & $(loc.start - ctg.start - width) & " ctg len:" & $ctg.len & " genomic position:" & $chrom & ":" & $loc.start)
        #quit()

      var refe = ref_kmer.mincode()
      var alte = alt_kmer.mincode()
      var alt_support = 0
      var ref_support = 0

      var both_found = false
      for read in r.reads:
        var ref_found = false
        var alt_found = false
        for e in read.sequence(ss).slide(K):
            if e == refe: ref_support += 1; ref_found = true
            if e == alte: alt_support += 1; alt_found = true
        if ref_found and alt_found:
            both_found = true
            break

      if alt_support < min_reads or both_found: continue

      var gt = genotype(ref_support, alt_support, 1e-4)
      if gt.GT == GT.HOM_REF: continue
      var v = Variant(chrom: chrom, start: loc.start, genotype: gt, ref_kmer: ref_kmer,
                  alt_kmer: alt_kmer, AD: [ref_support, alt_support])

      if loc.event_type == Deletion:
        v.reference = fai.get(chrom, loc.start - 1, loc.stop - 1)
        v.alternate = v.reference[0..<1]
      else:
        v.alternate = ctg.sequence[(qloc.start - 1)..<(qloc.stop)]
        v.reference = v.alternate[0..<1]
      yield v

iterator event_locations(r: Record): event {.inline.} =
  ## the genomic start-end of the location of the event
  ## TODO: SA tags
  var off: int
  for c in r.cigar:
    var cons = c.consumes.reference
    if c.op != CigarOp.match:
      if cons:
        yield (r.start + off, r.start + off + c.len, c.len)
      else:
        yield (r.start + off, r.start + off + 1, c.len)
    if cons:
      off += c.len

proc overlaps(r: Record, start:int, stop:int): bool {.inline.} =
  if r.start > stop: return false
  if r.stop < start: return false
  return true

iterator gen_roi_internal(evidence: seq[uint8], cache:seq[Record], min_evidence:uint8, min_reads:int, max_reads:int, cache_start:int, cache_end:int): roi {.inline.} =
  ## given the counts in evidence and a region to look in (cache_start .. cache_end)
  ## yield regions where the values in evidence are >= min_evidence along with reads from that region.
  ## futher filter to those regoins that have > min_reads that overlap the putative event 
  ## TODO: filter out high-coverage regions with, e.g. 4 supporting reads out of 1000.
  var in_roi = false
  var roi_start = 0
  var roi_end = 0

  for i in cache_start..<cache_end:
    var ev = evidence[i]
    if ev >= min_evidence:
      if not in_roi:
        in_roi = true
        roi_start = i
      roi_end = i
      continue

    # ending a region, yield the roi
    if in_roi:
      var reads = new_seq_of_cap[Record](16)
      for r in cache:
        if r.overlaps(roi_start, roi_end):
          reads.add(r)
          if len(reads) > max_reads: break
        if r.start > roi_end: break
      if len(reads) >= min_reads and len(reads) <= max_reads:
        yield (roi_start, roi_end, reads)
      in_roi = false

  if in_roi:
    var reads = new_seq_of_cap[Record](16)
    for r in cache:
      if r.overlaps(roi_start, roi_end):
        reads.add(r)
        if len(reads) > max_reads: break
      if r.start > roi_end: break

    if len(reads) >= min_reads and len(reads) <= max_reads:
      yield (roi_start, roi_end, reads)

iterator gen_roi(b:Bam, t:Target, min_event_support:uint8=4, min_read_coverage:int=4, max_read_coverage:int=200): roi =
  # we iterate over the bam an increment an evidence counter in an genomic array for positions
  # that appear to have an event (any non match)
  # whenever we have a gap in coverage where start > last_end, we check the evidence
  # array for any regions with >= min_event_support that indicate an event.
  # we yield the bounds of the event and the reads that overlapped it.

  var evidence = new_seq[uint8](t.length + 1)
  var cache = new_seq_of_cap[Record](100000)
  # we use last_start to make sure we dont keep iterating over the same chunk of evidence.
  var last_start = 0

  for r in b.querys(t.name):
    if r.skippable: continue
    if cache.len > 0 and r.start > cache[cache.high].stop:
      for roi in gen_roi_internal(evidence, cache, min_event_support, min_read_coverage, max_read_coverage, last_start, r.start):
        yield roi
      # reset
      last_start = r.start
      cache.set_len(0)

    cache.add(r.copy())
    for e in r.event_locations:
      for i in e.start..<e.stop:
        evidence[i] += 1
        # reset after avoid overflow
        if evidence[i] == 0:
          evidence[i] = 255
  for roi in gen_roi_internal(evidence, cache, min_event_support, min_read_coverage, max_read_coverage, last_start, evidence.len):
    yield roi


## make proper vcf header.
proc make_contig_header(t:Target): string =
  return "##contig=<ID=$1,length=$2>" % [t.name, $t.length]

proc contig_header(b: Bam): string =
  return join(map(b.hdr.targets, make_contig_header), "\n")

when isMainModule:
  let version = "indelope 0.0.1"
  let doc = format("""
  $version

  Usage: indelope [options] <reference> <BAM-or-CRAM>

Arguments:

  <reference>     reference fasta file.
  <BAM-or-CRAM>   call variants in this file.

Options:

  -m --min-reads <INT>        minimum number of reads to send for alignment [default: 3]
  -c --min-contig-len <INT>   minimum contig length to send for alignment [default: 73]
  -e --min-event-len <INT>    minimum size of indel to report [default: 4]
  -t --threads <INT>          number of cram/bam decompression threads [default: 1]
  -h --help                   show help

  """ % ["version", version])

  let
     args = docopt(doc, version = version)
     min_reads = parse_int($args["--min-reads"])
     min_ctg_len = parse_int($args["--min-contig-len"])
     min_event_len = parse_int($args["--min-event-len"])
     bam_path = $args["<BAM-or-CRAM>"]
     fai = open_fai($args["<reference>"])

  var b:Bam
  open(b, bam_path, index=true, threads=parse_int($args["--threads"]))
  var targets = b.hdr.targets

  var ez = new_ez()
  var last_var: Variant
  echo header % [b.contig_header, "sample"]
  for target in targets:
    for r in gen_roi(b, target, min_read_coverage=min_reads, min_event_support=min_reads.uint8):
      for v in callsemble(r, fai, ez, min_ctg_len=min_ctg_len, min_reads=min_reads, min_event_len=min_event_len):
        if v.same(last_var): continue
        echo v
        last_var = v

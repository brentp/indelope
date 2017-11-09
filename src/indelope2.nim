import osproc
import contig
import binaryheap
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

proc skippable(r: Record): bool {.inline.} =
  if r.chrom == "hs37d5": return true
  if r.chrom.startswith("GL"): return true
  var f = r.flag
  if f.dup or f.qcfail or f.unmapped: return true
  if f.supplementary or f.secondary: return true
  return false


proc assemble(r: roi, fai: Fai, ez: Ez, bp_overlap:int=4) =
  var contigs = new_seq[Contig]()
  var sequence = ""
  var read_seq = ""
  var ss = ""
  # assemble the alt contig
  for read in r.reads:
    if read.skippable: continue
    contigs.insert(read.sequence(sequence), read.start)
  var k = 0
  var before = len(contigs)
  contigs = filter(contigs.combine(), proc(x: Contig): bool = x.nreads >= 4 and x.len > 75)
  if len(contigs) == 0: return

  echo "#########################################################"
  echo r.reads[0].chrom, ":", r.start, "-", r.stop, " n contigs:", len(contigs), " was:", before
  
  for ctg in contigs:
    var max_stop = ctg.start
    for read in r.reads:
      max_stop = max(max_stop, read.stop)

    var reference = fai.get(r.reads[0].chrom, ctg.start, max_stop)
    ctg.sequence.align_to(reference, ez)

    if ez.max_event_length < 3: continue
    var locs = toSeq(ez.event_locations(ctg.start))
    if locs.len > 4: continue

    echo "\n\n", ez.cigar_string(sequence),  " nreads:", ctg.nreads, " len c reads:", len(r.reads)
    echo "locs:", locs
    echo ez.draw(ctg.sequence, reference)
    echo ctg.support

    var ref_support = 0
    var alt_support = 0

    for read in r.reads:
      if read.stop < (locs[locs.high].stop + bp_overlap): continue
      if read.start > (locs[locs.high].start - bp_overlap): continue

      var offset = max(0, read.start - ctg.start)
      read.sequence(read_seq).align_to(ctg.sequence[offset..<ctg.len], ez)
      echo "\nread:", read.start, "..", read.stop
      echo "contig score:",  ez.score, " cigar:", ez.cigar_string(ss), " full cigar:", ez.cigar_string(ss, full=true)
      var cevents = toSeq(ez.event_locations(ctg.start + offset))
      echo "events:", cevents
      echo ez.draw(read_seq, ctg.sequence[offset..<ctg.len])
      var ctg_score = ez.score

      read_seq.align_to(reference[offset..<reference.len], ez)
      echo "ref score::", ez.score, " ", "cigar:", ez.cigar_string(ss), " full cigar:", ez.cigar_string(ss, full=true)
      echo "events:", toSeq(ez.event_locations(ctg.start + offset))
      echo ez.draw(read_seq, reference[offset..<reference.len])
      var ref_score = ez.score
      if ctg_score >= 50 and ctg_score > 3 + ref_score:
        alt_support += 1
      elif ref_score >= 50 and ref_score > 3 + ctg_score:
        ref_support += 1
      elif ref_score >= 50 and len(cevents) == 0:
        ref_support += 1
       
      
    echo ref_support, " ", alt_support
    quit()



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

iterator gen_roi_internal(evidence: seq[uint8], cache:seq[Record], min_evidence:uint8, min_reads:int, cache_start:int, cache_end:int): roi {.inline.} =
  ## given the counts in evidence and a region to look in (cache_start .. cache_end)
  ## yield regions where the values in evidence are >= min_evidence along with reads from that region.
  ## futher filter to those regoins that have > min_reads that overlap the putative event 
  ## TODO: filter out high-coverage regions with, e.g. 4 supporting reads out of 1000.
  var in_roi = false
  var roi_start = 0
  var roi_end = 0

  for i in cache_start..cache_end:
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
        if r.start > roi_end: break
      if len(reads) >= min_reads:
        yield (roi_start, roi_end, reads)
      in_roi = false

  if in_roi:
    var reads = new_seq_of_cap[Record](16)
    for r in cache:
      if r.overlaps(roi_start, roi_end):
        reads.add(r)
      if r.start > roi_end: break

    if len(reads) >= min_reads:
      yield (roi_start, roi_end, reads)

iterator gen_roi(b:Bam, t:Target, min_event_support:uint8=4, min_read_coverage:int=4): roi =
  # we iterate over the bam an increment an evidence counter in an genomic array for positions
  # that appear to have an event (any non match)
  # whenever we have a gap in coverage where start > last_end, we check the evidence
  # array for any regions with >= min_event_support that indicate an event.
  # we yield the bounds of the event and the reads that overlapped it.

  var evidence = new_seq[uint8](t.length)
  var cache = new_seq_of_cap[Record](100000)
  # we use last_start to make sure we dont keep iterating over the same chunk of evidence.
  var last_start = 0
  
  for r in b.querys(t.name):
    if r.skippable: continue
    if cache.len > 0 and r.start > cache[cache.high].stop:
      for roi in gen_roi_internal(evidence, cache, min_event_support, min_read_coverage, last_start, r.start):
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
  for roi in gen_roi_internal(evidence, cache, min_event_support, min_read_coverage, last_start, evidence.len):
    yield roi

when isMainModule:

  var b:Bam
  #open(b, "x.bam", index=true)
  open(b, "../mosdepth/GT04008021.bam", index=true)
  var targets = b.hdr.targets
  var fai = open_fai("/data/human/g1k_v37_decoy.fa")
  #[
  for r in b:
    if r.skippable: continue
    echo r.chrom, " ", r.start, " ", r.cigar, " ", toSeq(r.event_locations)
  ]#
  var ez = new_ez()
  for r in gen_roi(b, targets[0]):
    echo r.start, " ", r.stop, " ", r.reads.len
    assemble(r, fai, ez)

  #[
  let version = "indelope 0.0.1"
  let doc = format("""
  $version

  Usage: indelope [options] <prefix> <reference> <BAM-or-CRAM>

Arguments:

	<prefix>      output prefix.
  <BAM-or-CRAM> call variants in this file.

Options:
  
  -m --min-reads <INT>      minimum number of reads to send for alignment [default: 3]
  -c --min-contig-len <INT>      minimum contig length to send for alignment [default: 90]
  -h --help                       show help

  """ % ["version", version])

  let
     args = docopt(doc, version = version)
     min_reads = parse_int($args["--min-reads"])
     min_ctg_length = parse_int($args["--min-contig-len"])
     cram_path = $args["<BAM-or-CRAM>"]
     fai = open_fai($args["<reference>"])

  var ez = new_ez()
  for contig in assembler(cram_path, $args["<prefix>"], min_reads=min_reads, min_ctg_length=min_ctg_length):
    contig.align(fai, ez)

  #align($args["<prefix>"], $args["<reference>"])
  let bam = $args["<prefix>"] & ".bam"
  echo $bam
  echo "now find indels"
  ]#


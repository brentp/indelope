import osproc
import assemble
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
  cached_seq = tuple[start: int, stop: int, sequence: string]

  InfoContig* = ref object of RootObj
    ## InfoContig wraps a contig with additional info.
    ctg*: Contig
    chrom: string
    qreads*: Heap[cached_seq]


proc region(c: InfoContig): string =
  return c.chrom & ":" & $c.ctg.start & "-" & $c.stop & "(" & $c.qreads.size & ")"

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


proc dump_contigs(fq: File, ctgs: var Contigs, chrom: string, start: int, stop: int) =
  ctgs = ctgs.combine(p_overlap=0.5)
  for ctg in ctgs:
    ctg.trim()
  ctgs = ctgs.combine(p_overlap=0.5)
  var fq_line = new_string_of_cap(500)
  for nc, ctg in ctgs:
    fq.write(ctg.fastq(fq_line, name=chrom & ":" & $(start + 1) & "-" & $stop & "#" & $(nc + 1) & "_of_" & $ctgs.len))

proc align(prefix: string, reference: string) =
   var ret = execCmd("bash -c 'set -eo pipefail; bwa mem -k 21 -L 3 $1 $2.fastq | samtools sort -o $3.bam && samtools index $4.bam'" % [reference, prefix, prefix, prefix])
   if ret != 0:
     quit(ret)
  


proc skippable(r: Record): bool {.inline.} =
  if r.chrom == "hs37d5": return true
  if r.chrom.startswith("GL"): return true
  var f = r.flag
  if f.dup or f.qcfail or f.unmapped: return true
  if f.supplementary or f.secondary: return true
  return false

proc align(ctg:InfoContig, fai:Fai, ez:Ez) =
  var target = fai.get(ctg.chrom, ctg.ctg.start, ctg.stop)
  ctg.ctg.sequence.align_to(target, ez)
  var s = ""
  echo ctg.region, " ", ez.cigar_string(s), " ", ez.max_event_length

  #[
  var cigar_str = ""
  if aln.score > 0:
    #echo "contig to reference"
    var refCounts = 0
    var altCounts = 0
    var totCounts = 0
    for r in ctg.qreads:
      if r.stop < ctg.start: continue
      if r.start > ctg.stop: continue
      var q = r.sequence
      totCounts += 1
      var alnAlt = q.align_to(ctg.ctg.sequence, cfg)
      var alnRef = q.align_to(target, cfg)
      var altScore = alnAlt.score
      var refScore = alnRef.score
 
      #echo "reference:", refScore
      ##echo alnRef
      #echo "alternate:", altScore
      #echo alnAlt
      if altScore == refScore: continue
      if altScore > refScore:
        if float64(altScore) > 0.90 * float64(len(q)):
          altCounts += 1
      else:
        if float64(refScore) > 0.98 * float64(len(q)):
          refCounts += 1

    if true and float64(refCounts + altCounts) >= 0.5 * float64(totCounts):
      echo "#######################"
      echo ctg.chrom, ":", ctg.start + aln.start, " " & aln.cigar(cigar_str)
      echo $aln
      echo aln.score
      echo refCounts, " ", altCounts, " of total reads:", totCounts
    if ctg.start + aln.start == 227636811811:
      echo target
      echo ctg.ctg.sequence
      quit(1)
    ]#

iterator assembler(path: string, prefix: string, min_reads:int=5, min_ctg_length:int=150): InfoContig =
  var b: Bam
  open(b, path, index=true)
  var fq: File
  var fname = prefix
  fname = fname.strip(chars={'.'}) & ".fastq"

  if not open(fq, fname, fmWrite):
    quit(2)

  var sequence = ""
  var bqs = new_seq[uint8]()
  var ctgs : seq[Contig]
  var region : string
  var targets: seq[string]
  
  #region = "3:120971087-120971287"
  if region == nil:
    targets = map(b.hdr.targets, proc(a: Target): string = return a.name)
  else:
    targets = @[region]

  for target in targets:
    ctgs = new_seq[Contig]()
    var
      contig_start = -1000
      contig_stop = -1000
      in_region = false

    var qreads = new_seq_of_cap[cached_seq](1000)
 
    for r in b.querys(target):
      if r.skippable: continue
      discard r.sequence(sequence)
      discard r.base_qualities(bqs)
      var tmp = sequence
      tmp.trim(bqs)
      if tmp.len < 50: continue
      var is_informative = r.aux("SA") != nil or r.cigar.len > 1

      for i, contig in contigs:
        if r.start > contig.start + contig.sequence.length

      # too far, so start a new region
      if r.start > contig_stop and in_region:
        #fq.dump_contigs(ctgs, target, contig_start, contig_stop)
        for contig in ctgs:
          if contig.len < min_ctg_length: continue
          if contig.nreads < min_reads: continue
          yield InfoContig(ctg: contig, chrom:target, stop: contig_stop, qreads: qreads)
        #align_to_contigs(ctgs, target, qreads, contig_start, contig_end)

        # make a new set of contigs
        ctgs.set_len(0)
        in_region = false

      if is_informative:
        if not in_region:
          contig_start = r.start
          in_region = true

        discard ctgs.insert(tmp, r.start)
        contig_stop = max(contig_stop, r.stop)

      # add all reads that are long enough to the queue
      qreads.push((r.start, r.stop, tmp))

      while qreads.size > 0 and qreads.peek.stop < contig_start:
        discard qreads.pop()

    fq.dump_contigs(ctgs, target, contig_start, contig_stop)
    for contig in ctgs:
      if contig.len < min_ctg_length: continue
      if contig.nreads < min_reads: continue
      yield InfoContig(ctg: contig, chrom:target, stop:contig_stop, qreads: qreads)

  close(fq)
  #return fname


when isMainModule:
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

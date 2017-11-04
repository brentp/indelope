import osproc
import assemble
import queues
import sequtils
import os
import times
import random
import math
import algorithm
import strutils
import docopt
import hts

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

type cached_seq = tuple[start: int, stop: int, sequence: string]

proc dump_contigs(fq: File, ctgs: var Contigs, chrom: string, start: int, stop: int, min_reads:int, min_ctg_length:int) =
  ctgs = ctgs.combine(p_overlap=0.5)
  for ctg in ctgs:
    ctg.trim()
  ctgs = ctgs.combine(p_overlap=0.5)
  var fq_line = new_string_of_cap(500)
  for nc, ctg in ctgs:
    if ctg.len < min_ctg_length: break
    if ctg.nreads < min_reads: continue
    ctg.trim()
    if ctg.len < min_ctg_length: break
    fq.write(ctg.fastq(fq_line, name=chrom & ":" & $(start + 1) & "-" & $stop & "_" & $nc))

proc align(prefix: string, reference: string) =
   var ret = execCmd("bash -c 'set -eo pipefail; bwa mem -k 21 -L 3 $1 $2.fastq | samtools sort -o $3.bam && samtools index $4.bam'" % [reference, prefix, prefix, prefix])
   if ret != 0:
     quit(ret)
  

proc assembler(path: string, prefix: string, min_reads:int=5, min_ctg_length:int=150): string =
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
  
  #region = "2:40882229-40882499"
  if region == nil:
    targets = map(b.hdr.targets, proc(a: Target): string = return a.name)
  else:
    targets = @[region]

  for target in targets:
    ctgs = new_seq[Contig]()
    var q = new_seq_of_cap[cached_seq](20)
    var last_start = -1000
    var last_stop = -1000
    for r in b.querys(target):
      if r.chrom == "hs37d5": continue
      if r.chrom.startswith("GL"): continue
      if r.chrom.startswith("GL"): continue
      var f = r.flag
      if f.dup or f.qcfail or f.unmapped: continue
      if f.supplementary or f.secondary: continue

      discard r.sequence(sequence)
      discard r.base_qualities(bqs)
      var tmp = sequence
      tmp.trim(bqs)
      if tmp.len < 50: continue

      var is_informative = r.aux("SA") != nil or r.cigar.len > 1

      # too far, so start a new region
      if (r.start - last_stop) > 50:
        while q.len > 0:
          var cs = q.pop()
          # try to make the contigs a bit larger, here on the right side
          if cs.start > last_stop: break
          discard ctgs.insert(cs.sequence)
        fq.dump_contigs(ctgs, target, last_start, last_stop, min_reads=min_reads, min_ctg_length=min_ctg_length)

        # make a new set of contigs
        ctgs = new_seq[Contig]()
        if is_informative:
          last_start = r.start
          # extend the contig on the left size
          while q.len > 0:
            var cs = q.pop()
            if cs.stop < r.start:
              continue
            discard ctgs.insert(cs.sequence)

      if is_informative:
        discard ctgs.insert(tmp)
        last_stop = max(last_stop, r.stop)

      else: # even perfect reads can be used to extend contigs.

        while q.len > 30: discard q.pop
        q.add((r.start, r.stop, tmp))

    fq.dump_contigs(ctgs, target, last_start, last_stop, min_reads=min_reads, min_ctg_length=min_ctg_length)

  close(fq)
  return fname


when isMainModule:
  let version = "indelope 0.0.1"
  let doc = format("""
  $version

  Usage: indelope [options] <prefix> <BAM-or-CRAM>

Arguments:

	<prefix>      output prefix.
  <BAM-or-CRAM> call variants in this file.

Options:
  
  -r --reference <STR>      path to bwa indexed reference
  -m --min-reads <INT>      minimum number of reads to send for alignment [default: 3]
  -c --min-contig-len <INT>      minimum contig length to send for alignment [default: 120]
  -h --help                       show help

  """ % ["version", version])

  let
     args = docopt(doc, version = version)
     min_reads = parse_int($args["--min-reads"])
     min_ctg_length = parse_int($args["--min-contig-len"])
     cram_path = $args["<BAM-or-CRAM>"]

  let fastq = assembler(cram_path, $args["<prefix>"], min_reads=min_reads, min_ctg_length=min_ctg_length)
  align($args["<prefix>"], $args["--reference"])
  let bam = $args["<prefix>"] & ".bam"
  echo $bam
  echo "now find indels"

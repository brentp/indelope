import assemble
import queues
import sequtils
import os
import times
import random
import math
import algorithm
import strutils

import hts

iterator random_regions(b:Bam, n_per:int=7000, n:int=5000): Record =
  randomize(42)
  var length:int = 10000000
  var otargets = b.hdr.targets
  var targets = new_seq[Target]()
  var lens = new_seq[int]()
  var sum: int
  for i, t in otargets:
    if int(t.length) < length: continue
    targets.add(t)
    lens.add( sum + int(t.length))
    sum += int(t.length)

  for k in 0..<n:
    var off = random(sum)
    var idx = lowerBound(lens, off)
    var tgt = targets[idx]
    var toff = off - lens[idx]
    var nr: int
    for record in b.queryi(tgt.tid, toff, toff+length):
      if nr > n_per: break
      var f = record.flag
      if record.b.core.mtid != record.b.core.tid: continue
      if f.dup or f.unmapped or f.qcfail: continue
      if record.b.core.pos > record.b.core.mpos: continue
      if not f.properpair: continue
      yield record
      nr += 1

type insert_size_stats = tuple[mean:float64, std:float64]

proc mean_std(arr: seq[float64]): insert_size_stats =
  var m, std: float64
  var L = float64(arr.len)
  for a in arr:
    m += a / L

  for a in arr:
    std += pow(a - m, 2) / L

  return (m, sqrt(std))

proc insert_size(r: Record): float64 =
  return float64(r.b.core.mpos - r.stop)

proc get_insert_size(b: Bam): insert_size_stats =
  var insert_sizes = new_seq_of_cap[float64](10000)
  for record in b.random_regions(n=3000):
    insert_sizes.add(record.insert_size())
  return mean_std(insert_sizes)


proc write_to(r:Record, fh:File, sequence:var string, baseqs:var seq[uint8]) =
  fh.write_line("@" & r.qname & " " & r.chrom & ":" & intToStr(r.start+1) & "-" & intToStr(r.stop))
  fh.write_line(r.sequence(sequence))
  fh.write_line("+")
  discard r.base_qualities(baseqs)
  for i, v in baseqs:
    baseqs[i] = v + 33
  fh.write_line(cast[string](baseqs))

type se = tuple[start:int, stop:int]

iterator to_intervals(arr: seq[bool]): se =
  var gt0: bool = false
  var start:int = -1
  for i, b in arr:
    if b != gt0:
      if gt0 and start >= 0:
        yield (start, i)
      else:
        start = i
      gt0 = b
  if gt0:
    yield (start, arr.len)


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

proc assemble(tid:int, arr:seq[bool], b:Bam) =
  var counter: int
  var sequence = ""
  var bqs = new_seq[uint8]()
  var chrom = b.hdr.targets[tid].name
  for iv in arr.to_intervals():
    var ctgs = new_seq[Contig]()
    var k: int
    var t = cpuTime()
    for r in b.queryi(tid, iv.start, iv.stop):
      if r.cigar.len == 1 and not r.flag.unmapped: continue

      discard r.sequence(sequence)
      #echo r.qname
      #echo r.base_qualities(bqs)
      var tmp = sequence
      discard ctgs.insert(tmp)
      k += 1
    # TODO: for each contig, extend the flanks with good reads to improve alignment.
    if k < 4: continue
    ctgs = ctgs.combine(p_overlap=0.5)
    for ctg in ctgs:
      ctg.trim()
    ctgs = ctgs.combine(p_overlap=0.5)
    if ctgs[0].len < 125: continue

    var fq = ""
    for nc, ctg in ctgs:
      ctg.trim()
      if ctg.len < 125: break
      write(stdout, ctg.fastq(fq, name=chrom & ":" & $(iv.start) & "-" & $(iv.stop) & "_" & $nc))

const pad:int = 100

type cached_seq = tuple[start: int, stop: int, sequence: string]

proc dump_contigs(ctgs: var Contigs, chrom: string, start: int, stop: int) =
  ctgs = ctgs.combine(p_overlap=0.5)
  for ctg in ctgs:
    ctg.trim()
  ctgs = ctgs.combine(p_overlap=0.5)
  var fq = ""
  for nc, ctg in ctgs:
    if ctg.len < 125: break
    ctg.trim()
    if ctg.len < 125: break
    write(stdout, ctg.fastq(fq, name=chrom & ":" & $(start + 1) & "-" & $stop & "_" & $nc))

proc assembler(path: string, opath: string) =
  var b: Bam
  open(b, path, index=true)

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
        dump_contigs(ctgs, target, last_start, last_stop)

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



    dump_contigs(ctgs, target, last_start, last_stop)




proc find_reads_of_interest(path:string, opath:string) =
  var b1:Bam
  open(b1, path, index=true)
  if b1.idx == nil:
    quit(2)

  var fastq: File
  if not open(fastq, opath & ".fastq", fmWrite):
    quit(2)

  #var isize = get_insert_size(b1)
  #echo isize
  var b:Bam
  open(b, path, index=true)
  #var imax = isize.mean + 7 * isize.std
  var last_tid = -1
  var targets = b.hdr.targets
  var arr: seq[bool]

  for r in b:
    if r.b.core.tid != last_tid:
      if last_tid != -1:
        var k = 0
        assemble(last_tid, arr, b1)
      last_tid = r.b.core.tid
      arr = new_seq[bool](targets[last_tid].length)
    if r.chrom == "hs37d5": continue
    if r.chrom.startswith("GL"): continue
    var f = r.flag
    if f.dup or f.qcfail or f.unmapped: continue
    if f.supplementary or f.secondary: continue
    if r.aux("SA") != nil or r.cigar.len > 1: # or abs(r.insert_size) > imax:
      for p in max(0, r.start-pad)..<min(r.stop+pad, arr.len):
        arr[p] = true
      #r.write_to(fastq, sequence, bqs)

  #for p in (69200-pad)..(69900+pad):
  #  arr[p] = true
  assemble(last_tid, arr, b1)

when isMainModule:

  #find_reads_of_interest(commandLineParams()[0], "xx")
  assembler(commandLineParams()[0], "xx")
  

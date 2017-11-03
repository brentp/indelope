import assemble
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
  fh.write_line("@" & r.qname & " " & r.chrom & ":" & intToStr(r.start) & "-" & intToStr(r.stop))
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

proc assemble(tid:int, arr:seq[bool], b:Bam) =
  var sequence = ""
  var nn = 0
  for iv in arr.to_intervals():
    var ctgs = new_seq[Contig]()
    var k: int
    var t = cpuTime()
    for r in b.queryi(tid, iv.start, iv.stop):
      if r.cigar.len == 1 and not r.flag.unmapped: continue
      discard r.sequence(sequence)
      var tmp = sequence
      discard ctgs.insert(tmp)
      k += 1
    # TODO: for each contig, extend the flanks with good reads to improve alignment.
    if k < 4: continue
    var cmb = ctgs.combine(p_overlap=0.5)
    if cmb[0].len < 125: continue
    #echo ">time:", cpuTime() - t
    #echo ">", tid, " ", iv, " alignments:", k, " ctgs:", len(ctgs), " combined:", len(cmb), " longest:", cmb[0].len, " nreads:", cmb[0].nreads
    for ctg in cmb:
      if ctg.len < 125: break
      echo "@seq" & $nn & ":" & b.hdr.targets[tid].name & ":" & $iv.start & "-" & $iv.stop
      echo ctg.sequence
      echo "+"
      echo join(cast[string](map(ctg.base_count, proc(i:uint32): char = return cast[char](min(90, int(i+33))) )))

      nn += 1
    if nn > 10:
      quit(0)
    #  echo ctg.base_count
    #  quit(0)
    #  assert ctg.sequence.len == ctg.base_count.len

const pad:int = 100

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
        echo targets[last_tid].name, ":", k
      last_tid = r.b.core.tid
      arr = new_seq[bool](targets[last_tid].length)
    if r.chrom == "hs37d5": continue
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

  echo(commandLineParams())
  find_reads_of_interest(commandLineParams()[0], "xx")
  

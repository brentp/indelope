import osproc
import contig
import sequtils
import strutils
import sets
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
    filter*: string
    qual*: float64
    reference*: string
    alternate*: string
    genotype*: Genotype
    ref_kmer*: string
    alt_kmer*: string
    info_str*: string
    AD*: array[2, int]

proc info(v:Variant): string =
  result = ("DP=" & $(v.AD[0] + v.AD[1]) &
           ";AD=" & $v.AD[0] & "," & $v.AD[1] &
           ";ref_kmer=" & $v.ref_kmer &
           ";alt_kmer=" & $v.alt_kmer)
  if v.info_str.len != 0:
    result &= ";" & v.info_str

proc info_add(v:Variant, kv:string) =
  ## add a new field to the info
  if v.info_str == nil or v.info_str.len == 0:
    v.info_str = kv
  else:
    v.info_str &= ';' & kv

const header = """##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=DP,Number=1,Type=Integer,Description="supporting k-mer depth">
##INFO=<ID=BS,Number=1,Type=Integer,Description="number of times there was support for both ref and alt k-mer in a single read">
##INFO=<ID=MF,Number=1,Type=Integer,Description="minimum matching bases around this event when BS > 0. Higher gives more confidence">
##INFO=<ID=CF,Number=1,Type=Integer,Description="minimum flank of the event from either end of the contig. higher is better.">
##INFO=<ID=NC,Number=1,Type=Integer,Description="number of contigs at the site of this variant.">
##INFO=<ID=CC,Number=1,Type=String,Description="contig cigar from alignment to reference">
##INFO=<ID=LO,Number=0,Type=Flag,Description="low-offset: the event occurred near at the start of the contig so we may not have the full variant">
##INFO=<ID=AKE,Number=1,Type=Float,Description="mean alt-kmer distance from end of read">
##INFO=<ID=RKE,Number=1,Type=Float,Description="mean ref-kmer distance from end of read">
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
  if v.filter == nil: v.filter = ""
  var filter = v.filter
  if filter == "": filter = "PASS"
  return format("$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\tGT:GQ:GL\t$gt" %
                ["chrom", v.chrom, "pos", $v.start, "id", ".",
                 "ref", v.reference, "alt", v.alternate,
                 "qual", formatFloat(v.qual, precision=2, format=ffDecimal),
                 "filter", filter, "info", v.info(), "gt", $(v.genotype)])

proc same(a:Variant, b:Variant): bool =
  if a == nil or b == nil: return false
  return (a.start == b.start) and (a.chrom == b.chrom) and (a.reference == b.reference) and (a.alternate == b.alternate)

const init_len = 100000000
proc get_min_flank(e:ksw2.event, ez:Ez): int =
  result = init_len
  var found_event = false
  for c in ez.cigar:
    if c.op == 0:
      if found_event:
        result = min(c.length.int, result)
      else:
        result = c.length.int
      if found_event: return
    elif c.op != 0 and int(c.op) - 1 == int(e.event_type) and c.length.int == e.len.int:
      if init_len == result: result = 0
      found_event = true
  return 0

proc mean(a:seq[int]): float64 =
    result = 0'f64
    for v in a:
        result += v.float64
    return result / float64(a.len)

iterator callsemble(r: roi, fai: Fai, ez: Ez, min_ctg_len:int=74, min_reads:int=4, min_event_len:int=4, K:int=27): Variant =

  var contigs = new_seq[Contig]()
  var read_seq = ""
  var base_q = new_seq[uint8](300)
  # assemble the alt contig
  for read in r.reads:
    if read.skippable(allow_unmapped=true): continue
    discard read.sequence(read_seq)
    discard read.base_qualities(base_q)
    trim(read_seq, base_q)
    contigs.insert(read_seq, read.start, min_overlap=int(0.88 * float64(read_seq.len)))

  var n_contigs = len(contigs)
  when not defined(release):
    echo "n contigs:", n_contigs
  contigs = contigs.combine()
  var n_contigs_a = len(contigs)
  var sequence = ""
  var chrom = r.reads[0].chrom

  for ctg in contigs:
    if n_contigs > 20: continue
    #echo "ctg:", ctg.nreads, " len:", ctg.len
    if ctg.nreads < min_reads or ctg.len < min_ctg_len: continue

    var max_stop = ctg.start
    for read in r.reads:
      max_stop = max(max_stop, read.stop)

    var width = int((K + 1) / 2 - 1)
    # TODO: for an SV, for the alignment, we need both ends.
    var reference = fai.get(chrom, ctg.start, max_stop + K + 1)
    ctg.sequence.align_to(reference, ez)
    when not defined(release):
      echo ez.cigar_string(read_seq)

    #[
    discard ez.cigar_string(read_seq, full=true)
    echo "roi:", r.start, "-", r.stop
    echo "ctg range:", ctg.start, " ", ctg.start + ctg.len
    echo "full:", read_seq
    echo "part:", ez.cigar_string(read_seq, full=true)
    echo ez.draw(ctg.sequence, reference)
    ]#

    var qlocs = toSeq(ez.query_locations())
    if len(qlocs) == 0 or len(qlocs) > 3 or (qlocs[0].start != 0 and len(qlocs) > 2): continue
    ## TODO: revisit this.
    if len(qlocs) == 1 and qlocs[0].start == 0 and qlocs[0].event_type == Deletion and qlocs[0].len < 10: continue
    var ii = -1

    #echo "qlocs:", qlocs
    #echo "tlocs:", toSeq(ez.target_locations(ctg.start))

    for tloc in ez.target_locations(ctg.start):
      ii += 1
      if tloc.len.int < min_event_len: continue
      var tstart = max(0, tloc.start - ctg.start - width)
      if tstart + K > reference.len:
        tstart = reference.len - K

      var ref_kmer = reference[tstart..<(tstart + K)]
      var qloc = qlocs[ii]

      var low_offset = qloc.start < 5 or qloc.stop + 6 > ctg.len
      var offset = min(qloc.start, ctg.len - qloc.stop - 1)
      var qstart = max(qloc.start - width, 0)
      if qstart + K > ctg.len:
        qstart = ctg.len - K

      var alt_kmer = ctg.sequence[qstart..<(qstart + K)]
      # this can happen given a long homopolymer in the insert
      # so we put as much of the insertion into the kmer as we can an hope that makes it unique
      # TODO: this part has to be dynamic and iterative
      # 1. maximize counts
      # 2. minimize bs
      # 3. vary offset and k-mer length.
      #[
      echo qlocs
      for i in ez.cigar:
          echo i
      ]#
      if alt_kmer == ref_kmer:
        qstart = max(qloc.start - 3, 0)
        if qstart + K > ctg.sequence.len:
          var qend = min(qloc.stop + 4, ctg.len)
          alt_kmer = ctg.sequence[qend - K..<(qend)]
        else:
          alt_kmer = ctg.sequence[qstart..<(qstart + K)]
      if ref_kmer == alt_kmer and (qloc.start == 0 or alt_kmer.toSet.len == 1): continue
      # simple repeats are hard so require at least 2 unique bases.
      if ref_kmer.toSet.len < 3: continue

      if ref_kmer == alt_kmer:
        stderr.write_line("bug!!! ref and alt kmers are same!! chrom:" & chrom & " " & $qloc & " alt:" & $tloc)
        stderr.write_line(">>>>>> contig length:" & $ctg.len & " ctg start:" & $ctg.start & " ctg nreads:" & $ctg.nreads )
        stderr.write_line(">>>>>> kmer:" & ref_kmer)
        var d = ez.draw(ctg.sequence, reference).split("\n")
        stderr.write_line ">>>>>> " & d[0]
        stderr.write_line ">>>>>> " & d[1]
        continue

      if len(ref_kmer) != len(alt_kmer) or len(ref_kmer) != K:
        stderr.write_line("bug!!! kmer lengths should be equal, ref:" & $ref_kmer.len & " alt:" & $alt_kmer.len)
        stderr.write_line("bug!!! " & $qloc & " alt:" & $tloc)
        stderr.write_line("bug!!! tloc start:" & $(tloc.start - ctg.start - width) & " ctg len:" & $ctg.len & " genomic position:" & $chrom & ":" & $tloc.start)
        quit(2)

      var refe = ref_kmer.mincode()
      var alte = alt_kmer.mincode()
      var alt_support = 0
      var ref_support = 0
      var adists = new_seq_of_cap[int](len(r.reads))
      var rdists = new_seq_of_cap[int](len(r.reads))

      var both_found = 0
      for read in r.reads:
        #var e2 = new_ez()
        #discard read.sequence(read_seq)
        #read_seq.align_to(reference, e2)
        #echo toSeq(items(read.cigar)), ".. ", e2.cigar_string(read_seq)
        var ref_found = false
        var alt_found = false
        for d, e in read.sequence(read_seq).dists(K):
            if not ref_found and e == refe:
                ref_support += 1; ref_found = true
                rdists.add(d)

            if not alt_found and e == alte:
                alt_support += 1; alt_found = true
                adists.add(d)
        if ref_found and alt_found:
            both_found += 1

      if alt_support < min_reads: continue
      if float64(alt_support) / float64(alt_support + ref_support) < 0.1: continue

      var gt = genotype(ref_support, alt_support, 1e-3)
      if gt.GT == GT.HOM_REF: continue
      var v = Variant(chrom: chrom, start: tloc.start, genotype: gt, ref_kmer: ref_kmer,
                  qual: gt.qual, alt_kmer: alt_kmer, AD: [ref_support, alt_support])
      # this line works extremeley well to remove false positives with no loss of sensitivity.
      if offset == 0 and both_found >= int(0.75 * float64(min(ref_support, alt_support))): continue

      if low_offset:
        v.info_add("LO")
        v.qual /= 2'f64
      if both_found > 0:
        v.info_add("BS=" & $both_found)
        v.qual /= 1.5
      else:
        v.qual *= 2
      v.info_add("CC=" & ez.cigar_string(read_seq))
      var min_flank = qloc.get_min_flank(ez)
      # if we have a big event and a small flank, bail
      if (min_flank - 1) < max(tloc.stop - tloc.start, qloc.stop - qloc.start): continue
      v.info_add("MF=" & $min_flank)
      v.info_add("CF=" & $offset)
      v.info_add("NC=" & $n_contigs_a)
      if offset == 0:
          v.qual /= 4'f64
      v.info_add("AKE=" &  formatFloat(mean(adists), precision=2, format=ffDecimal))
      v.info_add("RKE=" &  formatFloat(mean(rdists), precision=2, format=ffDecimal))
      if mean(adists) < 5: continue
      if tloc.event_type == Deletion:
        v.reference = fai.get(chrom, tloc.start - 1, tloc.stop - 1)
        v.alternate = v.reference[0..<1]
      else:
        v.reference = fai.get(chrom, tloc.start - 1, tloc.start - 1)
        v.alternate = ctg.sequence[(qloc.start-1)..<qloc.stop]
        v.start = tloc.start
        var vset = v.alternate[1..v.alternate.high].toSet
        if vset.len == 1 and alt_kmer[alt_kmer.high-10..alt_kmer.high].toSet.len == 1 and 
            ref_kmer[ref_kmer.high-10..ref_kmer.high].toSet.len == 1:
          continue
      yield v

iterator event_locations(r: Record, max_tags:int=1): event {.inline.} =
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
  #var tags = r.aux("SA")
  #if tags != nil and tags.count(';') <= max_tags:
  #  tags = tags[0..<tags.high] # strip(';')
  #  for tag in tags:


proc overlaps(r: Record, start:int, stop:int): bool {.inline.} =
  if r.start > stop: return false
  if r.stop < start: return false
  return true

proc single_roi(b:Bam, region:string): roi =
  var reads = new_seq_of_cap[Record](16)
  for r in b.querys(region):
    reads.add(r.copy())
  var se = region.split(":")[1].split("-")
  return (parse_int(se[0]), parse_int(se[1]), reads)

iterator gen_roi_internal(evidence: seq[uint8], cache:seq[Record], min_evidence:uint8, min_reads:int, max_reads:int, cache_start:int, cache_end:int): roi {.inline.} =
  ## given the counts in evidence and a region to look in (cache_start .. cache_end)
  ## yield regions where the values in evidence are >= min_evidence along with reads from that region.
  ## further filter to those regions that have > min_reads that overlap the putative event
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

# store the max stop we've seen.
type cache_t = tuple[records: seq[Record], stop:int]

proc add(c:var cache_t, r:Record) =
  c.stop = max(c.stop, r.stop)
  c.records.add(r)

proc clear(c:var cache_t) =
  c.records.set_len(0)
  c.stop = 0

proc len(c:cache_t): int {.inline.} =
    return c.records.len

iterator gen_roi(b:Bam, t:Target, min_event_support:uint8=4, min_read_coverage:int=4, max_read_coverage:int=200): roi =
  # we iterate over the bam an increment an evidence counter in an genomic array for positions
  # that appear to have an event (any non match)
  # whenever we have a gap in coverage where start > last_end, we check the evidence
  # array for any regions with >= min_event_support that indicate an event.
  # we yield the bounds of the event and the reads that overlapped it.

  var evidence = new_seq[uint8](t.length + 1)
  var cache: cache_t = (new_seq_of_cap[Record](100000), 0)
  # we use last_start to make sure we dont keep iterating over the same chunk of evidence.
  var last_start = 0

  for r in b.querys(t.name):
    #echo "r.start:", r.start
    if cache.len > 0 and r.start > cache.stop:
      for roi in gen_roi_internal(evidence, cache.records, min_event_support, min_read_coverage, max_read_coverage, last_start, r.start):
        yield roi
      # reset
      last_start = int(r.start)
      cache.clear()

    if r.skippable: continue
    cache.add(r.copy())
    for e in r.event_locations:
      for i in e.start..<e.stop:
        evidence[i] += 1
        # reset after avoid overflow
        if evidence[i] == 0:
          evidence[i] = 255
  for roi in gen_roi_internal(evidence, cache.records, min_event_support, min_read_coverage, max_read_coverage, last_start, evidence.len):
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
  var b:Bam
  var ez = new_ez()

  if  len(commandLineParams()) > 0 and commandLineParams()[0] == "single-site":
    # single-site region fai bam
    open(b, commandLineParams()[3], index=true)
    var fai = open_fai(commandLineParams()[2])
    var r = single_roi(b, commandLineParams()[1])
    echo "got ", len(r.reads), " reads"
    for v in callsemble(r, fai, ez, min_event_len=4):
      echo v
    quit(0)
  let
     args = docopt(doc, version = version)
     min_reads = parse_int($args["--min-reads"])
     min_ctg_len = parse_int($args["--min-contig-len"])
     min_event_len = parse_int($args["--min-event-len"])
     bam_path = $args["<BAM-or-CRAM>"]
     fai = open_fai($args["<reference>"])

  open(b, bam_path, index=true, threads=parse_int($args["--threads"]))
  var targets = b.hdr.targets

  var last_var: Variant
  echo header % [b.contig_header, "sample"]
  for target in targets:
    for r in gen_roi(b, target, min_read_coverage=min_reads, min_event_support=max(3, min_reads-2).uint8):
      for v in callsemble(r, fai, ez, min_ctg_len=min_ctg_len, min_reads=min_reads, min_event_len=min_event_len):
        if v.same(last_var): continue
        echo v
        last_var = v

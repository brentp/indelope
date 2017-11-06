import os
import strutils
import sequtils
import hts

discard """
2:40882328-40882721_0#nreads=6#len=130	0	2	40882429	60	70M60S	*	0	0	TTTTCTCAGCATATGTGTGTGTATGCATGTATAATGCATGTACTATATGCAGGTATATTTAAAATTGTTTAGAATTAAGTTGCCATGGCAGTTTTAGAAGAGTTTGACTGATTCTATTGACACACATAAT	#################################%%%%%%%%%%%%%%%%%%%%%%%''''''''''''''&&&&$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######################	NM:i:0	MD:Z:70	AS:i:70	XS:i:0	SA:Z:2,40882559,+,67S63M,60,0;
"""

type
  Info* = ref object of RootObj
    chrom*: string
    start*: int
    stop*:int
    nreads*: int
    nof*:string

  Variant* = ref object of RootObj
    chrom*: string
    start*: int
    stop*: int
    reference*: string
    alternate*: string
    quality*: uint8
    INFO: string
    AD*: array[2, int]
    GT*: string
    GQ*: int

proc info(r: Record): Info =
  ## information was encoded in the read name. extract it.
  var rnrl = r.qname.split("#", 4)
  var region = rnrl[0]
  var cse = region.split(":", 1)
  var se = cse[1].split("-", 1)
  var nof = rnrl[1]
  var nr = rnrl[2].split("=")[1]
  return Info(chrom: cse[0], start: parseInt(se[0]) - 1, stop: parseInt(se[1]), nreads: parseInt(nr), nof:nof)

proc `$`*(i: Info): string =
  return "Info(" & i.chrom & ":" & $(i.start + 1) & "-" & $(i.stop) & " nreads:" & $i.nreads & ")"

proc likely_variant(r: Record, i:Info): bool =
  if r.chrom != i.chrom: return false # TODO: might mis translocations.

  #if r.start > (i.stop + 100): return false
  #if i.start > (r.stop + 100): return false
  if r.qual < 5: return false
  if i.nreads < 3: return false

  var sa = r.aux("SA")
  if r.cigar.len == 1 and sa == nil: return false
  if r.cigar.len > 5: return false

  var nm = r.aux("NM")
  var insertion_bases = 0
  for o in r.cigar:
      if o.op in {CigarOp.insert, CigarOp.deletion}:
        insertion_bases += o.len 
  # don't allow more than 2 mismatches. insertions are counted as mismatches.
  if nm != nil and nm.integer() > (insertion_bases + 2): return false

  var min_len = 200
  for op in r.cigar:
    if op.len < min_len:
      min_len = op.len
  if min_len < 3: return false
  var s = "nil"
  if sa != nil:
    s = sa.tostring()

  return true

iterator as_sa_variant(r: Record, sa: string, loc:Info, bqs: seq[uint8]): Variant =
  # 2:40882328-40882721_0#nreads=6#len=130  0 2 40882429  60  70M60S  * 0 0 TTTTCTCAGCATATGTGTGTGTATGCATGTATAATGCATGTACTATATGCAGGTATATTTAAAATTGTTTAGAATTAAGTTGCCATGGCAGTTTTAGAAGAGTTTGACTGATTCTATTGACACACATAAT  #################################%%%%%%%%%%%%%%%%%%%%%%%''''''''''''''&&&&$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######################  NM:i:0  MD:Z:70 AS:i:70 XS:i:0  SA:Z:2,40882559,+,67S63M,60,0;
  var sas = sa.strip(chars={';'}).split(';')
  # using while loop so we can break as no StopIteration in nim
  while sas.len == 1: # only do stuff with a single split
    var v = Variant(chrom:r.chrom, start: r.stop)
    var sat = sas[0].split(',')
    if sat[0] != r.chrom:
      break

    if parse_int(sat[4]) < 10:
      break

    var svstop = sat[1]
    v.stop = parse_int(svstop)
    v.reference = "N"

    if v.stop > v.start:
      v.alternate = "<DEL>"
    else:
      echo "FIX ME ", r
      break

    var alt_left = float64(bqs[r.stop - r.start])
    var alt_right = float64(bqs[r.stop - r.start + 1])
    v.AD[0] = 0 # TODO
    v.AD[1] = loc.nreads
    v.AD[1] = int(0.51 + ((alt_left + alt_right) / 2.0))
    yield v
    break

proc add_ref_alt(v:var Variant, vstart:int, vstop:int, sequence: string, o:Op, fa:Fai) =
  if not o.consumes.reference:
    v.alternate = sequence[(vstart - 1)..<(vstart + o.len)]
    #v.reference = fa.get(r.chrom, v.start, v.start)
    v.reference = v.alternate[0..<1]
  else:
    # subtract 1 because we need an extra anchor base for the start
    v.reference = fa.get(v.chrom, v.start - 1, v.stop - 1)
    v.alternate = v.reference[0..<1]

proc createVariant(vstart:int, vstop:int, r:Record, bqs: seq[uint8], loc:Info, o:Op, fa:Fai): Variant =
  var sequence: string = ""
  var v = Variant(chrom: r.chrom, start: vstart + r.start, stop: vstop + r.start)
  discard r.sequence(sequence)
  v.add_ref_alt(vstart, vstop, sequence, o, fa)
  v.quality = uint8(r.qual)
  v.AD[0] = 0 # TODO:
  v.AD[1] = loc.nreads
  v.INFO = "nof=" & $loc.nof & ";"
  return v

iterator as_variant*(r: Record, loc: Info, bqs: seq[uint8], fa: Fai): Variant =
  var sa = r.aux("SA")
  if sa != nil:
    discard sa.tostring
    #for res in r.as_sa_variant(sa.tostring, loc, bqs):
    #  yield res
    
  var vstart = 0
  for o in r.cigar:
    if o.op in {CigarOp.insert, CigarOp.deletion}:
      var vstop = vstart
      if o.consumes.reference:
        vstop += o.len

      yield createVariant(vstart, vstop, r, bqs, loc, o, fa)
    if o.consumes.reference:
        vstart += o.len

var header = """##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=n_of,Number=1,Type=String,Description="indicates how many contigs were created for this region and which was ##the alternate">
$1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""
#

proc make_contig_header(t:Target): string =
  return "##contig=<ID=$1,length=$2>" % [t.name, $t.length]

proc contig_header(b: Bam): string =
  return join(map(b.hdr.targets, make_contig_header), "\n")

when isMainModule:

  var b:Bam
  var bqs = new_seq[uint8]()
  open(b, commandLineParams()[0], index=true)
  var fai = open_fai(commandLineParams()[1])
  echo(header % [b.contig_header])
  #.query("16", 29088058, 29088060):
  for rec in b:
    if rec.flag.secondary: continue
    if rec.flag.supplementary: continue
    var loc = rec.info
    if not rec.likely_variant(loc): continue
    discard rec.base_qualities(bqs)

    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
    for v in rec.as_variant(loc, bqs, fai):
      echo(v.chrom, "\t", $(v.start), "\t",
        "ID\t",
        v.reference, "\t",
        v.alternate, "\t",
        $(v.quality), "\t",
        "PASS\t",
        v.INFO,
        "END=" & $(v.stop),
        ";SVLEN=" & $(v.stop - v.start),
        ";AD=" & $v.AD[0] & "," & $v.AD[1])

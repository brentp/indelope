import strutils
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

  Variant* = ref object of RootObj
    chrom*: string
    pos*: int
    stop*: int
    reference*: string
    alternate*: string
    quality*: string
    AD*: array[2, int]
    GT*: string
    GQ*: int

proc info(r: Record): Info =
  var rnrl = r.qname.split("#", 4)
  var region = rnrl[0]
  var cse = region.split(":", 1)
  var se = cse[1].split("-", 1)
  var nr = rnrl[1].split("=")[1]
  return Info(chrom: cse[0], start: parseInt(se[0]) - 1, stop: parseInt(se[1]), nreads: parseInt(nr))

proc `$`*(i: Info): string =
  return "Info(" & i.chrom & ":" & $(i.start + 1) & "-" & $(i.stop) & " nreads:" & $i.nreads & ")"

proc likely_variant(r: Record, i:Info): bool =
  if r.chrom != i.chrom: return false # TODO: might mis translocations.

  if r.start > (i.stop + 100): return false
  if i.start > (r.stop + 100): return false
  if r.qual < 10: return false
  if i.nreads < 3: return false

  var sa = r.aux("SA")
  if r.cigar.len == 1 and sa == nil: return false

  var nm = r.aux("NM")
  if nm != nil and nm.integer() > 2: return false

  var min_len = 200
  for op in r.cigar:
    if op.len < min_len:
      min_len = op.len
  if min_len < 4: return false
  var s = "nil"
  if sa != nil:
    s = sa.tostring()

  return true

iterator as_sa_variant(r: Record, sa: string, loc:Info, bqs: seq[uint8]): Variant =
  # 2:40882328-40882721_0#nreads=6#len=130  0 2 40882429  60  70M60S  * 0 0 TTTTCTCAGCATATGTGTGTGTATGCATGTATAATGCATGTACTATATGCAGGTATATTTAAAATTGTTTAGAATTAAGTTGCCATGGCAGTTTTAGAAGAGTTTGACTGATTCTATTGACACACATAAT  #################################%%%%%%%%%%%%%%%%%%%%%%%''''''''''''''&&&&$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#######################  NM:i:0  MD:Z:70 AS:i:70 XS:i:0  SA:Z:2,40882559,+,67S63M,60,0;
  var sas = sa.strip(chars={';'}).split(';')
  # using while loop so we can break as no StopIteration in nim
  while sas.len == 1: # only do stuff with a single split
    var v = Variant(chrom:r.chrom, pos: r.stop)
    var sat = sas[0].split(',')
    if sat[0] != r.chrom:
      break

    if parse_int(sat[4]) < 10:
      break

    var svstop = sat[1]
    v.stop = parse_int(svstop)
    v.reference = "N"

    if v.stop > v.pos:
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


iterator as_variant*(r: Record, loc: Info, bqs: seq[uint8]): Variant =
  var sa = r.aux("SA")
  if sa != nil:
    for res in r.as_sa_variant(sa.tostring, loc, bqs):
      yield res
  else: # TODO: sa reads might also have other stuff
    echo r.qname, " ", $r.cigar, " ", r.chrom, ":", r.start

when isMainModule:

  var b:Bam
  var bqs = new_seq[uint8]()
  open(b, "t.bam", index=true)
  #.query("16", 29088058, 29088060):
  for rec in b:
    if rec.flag.secondary: continue
    if rec.flag.supplementary: continue
    var loc = rec.info
    if not rec.likely_variant(loc): continue
    discard rec.base_qualities(bqs)

    for v in rec.as_variant(loc, bqs):
      echo "YAY:", v.chrom, " ", v.pos, " ", v.stop, " ", v.AD[1]

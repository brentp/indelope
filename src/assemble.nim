import strutils
import sequtils
import algorithm
import sets

const complement = ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']

type 
  Contig* = ref object of RootObj
    ## a contig is a collection of sequences that have been
    ## merged. the base_count indicates the number of reads supporting
    ## the contig at each base.
    sequence*: string
    base_count*: seq[uint32]
    mismatch_count*: seq[uint32]
    nreads*: int

  Match* = ref object of RootObj
    matches: int
    added: bool
    offset: int
    mm: int
    sub: string
    contig: Contig

proc revcomp(c: var Contig) =
  c.base_count.reverse
  c.sequence.reverse
  c.mismatch_count.reverse
  for i, ch in c.sequence:
    c.sequence[i] = complement[ch.int]

var DEBUG* = false

proc match_sort(a, b: Match): int =
  if a.matches == b.matches:
    return b.mm - a.mm
  return a.matches - b.matches

type Contigs* = seq[Contig]
type Matches* = seq[Match]

proc len*(c:Contig): int {.inline.} =
  ## the length of the sequence of the contig
  return c.sequence.len

proc high*(c:Contig): int {.inline.} =
  return c.sequence.len - 1

proc trim*(c:Contig, min_support:int=2) =
  ## trim bases in a contig that do not have at least `min_support`
  var a = 0
  while a < c.high and c.base_count[a] < uint32(min_support):
    a += 1

  if a >= c.high:
    c.sequence.set_len(0)
    c.base_count.set_len(0)
    c.mismatch_count.set_len(0)
    return

  var b = c.high
  while c.base_count[b] < uint32(min_support) and b > a:
    b -= 1

  if a > 0 or b <= c.high:
    c.base_count = c.base_count[a..b]
    c.sequence = c.sequence[a..b]
    c.mismatch_count = c.mismatch_count[a..b]

proc fastq*(c:Contig, s:var string, name:string=""): string =
  ## return the 4 line fastq string for the record.
  ## the base-quality is the number of reads supporting each base.
  s.set_len(0)
  s.add "@" & name & "#"
  s.add "nreads=" & $c.nreads & "#"
  s.add "len=" & $c.len
  s.add "\n"

  s.add c.sequence
  s.add "\n"

  s.add "+\n"

  s.add join(cast[string](map(c.base_count, proc(i:uint32): char = return cast[char](min(90, int(i+33))) )))
  s.add "\n"

  result = s


const line_len = 160

proc `$`*(c: Contig): string =
  var i = 0
  result = new_string_of_cap(256)
  result.add("length: $#, reads: $#\n" % [$c.len, $c.nreads])
  while i < c.len:
    var s = c.sequence[i..<min(c.len, i+line_len)]
    var b = c.base_count[i..<min(c.len, i+line_len)]
    i += line_len
    result.add(s)
    result.add("\n")
    var found = true
    var pow = 0
    while found:
      found = false
      for base in b:
        var bc = intToStr(int(base))
        if len(bc) > pow:
          found = true
          result.add(bc[pow])
        else:
          result.add(" ")
      result.add("\n")
      pow += 1
    result.add("\n")

proc count_matches*(c: Contig, dna: var Contig, min_overlap:int=40, max_mismatch:int=0, p_overlap:float64=0.7): Match =
  ## add a dna sequence to the contig. return indicates if it was added.
  if c.sequence == nil:
    c.sequence = dna.sequence
    c.nreads = dna.nreads
    c.base_count = dna.base_count
    c.mismatch_count = dna.mismatch_count
    return Match(matches: dna.len, added: true)

  var min_mm = max_mismatch + 1
  var max_ma = 0
  var max_ma_offset = 0

  var sub: string
  var max_sub: string

  # slide the read along the contig and find where it matches best
  for read_offset in -(dna.len - min_overlap)..(c.len - min_overlap):
    if read_offset < 0:
      # read overhangs contig on left, so chop
      sub = dna.sequence[-read_offset..dna.high]
    elif (read_offset + dna.len) > c.len:
      # read overhangs contig on right, so chop
      sub = dna.sequence[0..<(c.len - read_offset)]
    else:
      sub = dna.sequence

    var mm = 0
    var ma = 0
    for i, s in sub:
      if i + max(0, read_offset) == c.sequence.len: break
      #echo s, c.sequence[i + max(0, read_offset)]

      if s != c.sequence[i + max(0, read_offset)]:
        mm += 1
      #elif i + read_offset >= 0 and s == c.sequence[i + read_offset]:
      elif s == c.sequence[i + max(0, read_offset)]:
        ma += 1
      if mm > max_mismatch:
        break

    if mm <= max_mismatch and ma > max_ma:
      max_ma_offset = read_offset
      max_ma = ma
      min_mm = mm
      max_sub = sub

  if min_mm > max_mismatch or float64(max_ma) / float64(min(dna.len, c.len)) < p_overlap:
    return nil

  return Match(matches: max_ma, added: false, contig: c,
      offset: max_ma_offset, mm: min_mm, sub: max_sub)

proc insert*(c: Contig, dna: var Contig, match:Match=nil, min_overlap:int=40, max_mismatch:int=0, p_overlap:float64=0.7): bool =
  var matches = match
  if matches == nil:
    matches = c.count_matches(dna, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap)
    if matches == nil: return false
    if matches.added: return true

  if DEBUG:
    var x = ""
    for i in 0..<abs(matches.offset):
      x &= " "
    echo "BEFORE"
    if matches.offset < 0:
      echo "CTG :", x & c.sequence
      echo "READ:", dna.sequence
    else:
      echo "CTG :", c.sequence
      echo "READ:", x & dna.sequence

  if matches.mm > max_mismatch:
    return false

  if matches.offset >= 0:
    # add sequence ot the right end (or middle).
    var extra = 0
    if (matches.offset + dna.len) > c.sequence.len:
      # the added sequence extends the contig 
      extra = matches.offset + dna.len - c.sequence.len
      c.sequence &= dna.sequence[(dna.len - extra)..<dna.len]
    if c.sequence.len > c.base_count.len:
      # this should always be true
      var L = c.base_count.len
      c.base_count.set_len(c.sequence.len)
      c.mismatch_count.set_len(c.sequence.len)
      for i in L..<c.base_count.len:
        c.base_count[i] = 1

    #for i in max_ma_offset..<min(c.len, dna.len - max_ma_offset):
    for i in matches.offset..<min(c.len-extra, dna.len + matches.offset):
      if c.sequence[i] == dna.sequence[i - matches.offset]:
        #echo i, " adding ", c.len, " ", i - max_ma_offset, " in ", dna.len
        c.base_count[i] += dna.base_count[i - matches.offset]
      elif dna.base_count[i - matches.offset] > c.base_count[i]:
        # if the incoming contig has more evidence for a different base, we adopt it.
        c.sequence[i] = dna.sequence[i - matches.offset]
        c.base_count[i] = dna.base_count[i - matches.offset]
  else:
    # add new sequence to the left end
    c.sequence = dna.sequence[0..<(-matches.offset)] & c.sequence
    # if the dna.sequence completely encompasses the contig sequence, we have to added
    # to the right end as well:
    if dna.sequence.len > c.sequence.len:
      var d = dna.sequence.len - c.sequence.len
      c.sequence &= dna.sequence[dna.sequence.len-d..<dna.sequence.len]

    if c.sequence.len > c.base_count.len:
      c.base_count.set_len(c.sequence.len)
      c.mismatch_count.set_len(c.sequence.len)
      # move all the counts to the left
      for i in countdown(c.len-1, -matches.offset):
        c.base_count[i] = c.base_count[i + matches.offset]
        c.mismatch_count[i] = c.mismatch_count[i + matches.offset]
      # zero out the newest.
      for i in countdown(-matches.offset-1, 0):
        c.base_count[i] = 0
        c.mismatch_count[i] = 0
      # increment the match.
      for i in 0..<dna.len:
        if c.sequence[i] != dna.sequence[i]:
          c.mismatch_count[i] += 1
        if c.sequence[i] == dna.sequence[i]:
          # increment by incoming count
          c.base_count[i] += dna.base_count[i]
        elif dna.base_count[i] > c.base_count[i]:
          # in case of mismatch, the most popular base_count determines the sequence.
          c.sequence[i] = dna.sequence[i]
          c.base_count[i] = dna.base_count[i]

  # ERROR correction
  # now track the mismatches so we can do a crude error correction by voting.
  var found = false
  var rseq = c.sequence
  if matches.offset < 0:
    rseq = rseq[-matches.offset..rseq.high]

  for i, s in matches.sub:
      if i == rseq.len: break
      if s != rseq[i]:
        echo max_mismatch
        echo matches.sub
        echo rseq
        #echo c.sequence[max(0, matches.offset)..c.sequence.high]
        echo matches.offset
        echo matches.mm
        quit()
        c.mismatch_count[i + max(0, matches.offset)] += 1
        found = true
  if found:
    for i, s in matches.sub:
      var cidx = i + max(0, matches.offset)
      if cidx == c.sequence.len: break
      if c.mismatch_count[cidx] > c.base_count[cidx]:
        # switch to more likely base
        var tmp = c.mismatch_count[cidx]
        c.mismatch_count[cidx] = c.base_count[cidx]
        c.base_count[cidx] = tmp
        c.sequence[cidx] = s

  # end ERROR correction
  c.nreads += dna.nreads

  if DEBUG:
    echo "AFTER"
    echo c
  return true

proc make_contig(dna: string): Contig =
  ## make a contig from a sequence string.
  var bc = new_seq[uint32](dna.len)
  for i in 0..bc.high:
    bc[i] = 1
  var o = Contig(sequence: dna, base_count:bc, nreads: 1)
  o.mismatch_count = new_seq[uint32](dna.len)
  return o

proc insert*(c: Contig, dna: string, min_overlap:int=48, max_mismatch:int=0, p_overlap:float64=0.75): bool {.inline.} =
  ## add a sequence to a contig.
  ## return value indicates that it was added to an existing contig
  var o = make_contig(dna)
  return c.insert(o, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap)

proc count_matches(ctgs: Contigs, o: var Contig, min_overlap:int=40, max_mismatch:int=0, p_overlap:float64=0.7): seq[Match] =
  var matches = new_seq_of_cap[Match](ctgs.len)
  for c in ctgs:
    var ma = c.count_matches(o, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap)
    if ma == nil: continue
    if ma.added: return @[ma]
    matches.add(ma)
    # we are appending to a seq to find the best match, but if we find one that's very good, we bail early.
    if ma.mm == 0 and ma.matches > (o.len - 5) and ma.matches > 50:
      return @[ma]

  matches.sort(match_sort)
  return matches

proc insert*(ctgs:var Contigs, o:var Contig, min_overlap:int=40, max_mismatch:int=0, p_overlap:float64=0.7): bool =
  ## insert a contig into the best match out of a set of contigs. If a suitable contig is not found, create a new one.
  ## return value indicates that it was added to an existing contig
  var matches = ctgs.count_matches(o, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap)
  if len(matches) == 0:
    ctgs.add(o)
    return false
  if len(matches) == 1 and matches[0].added:
    return false

  if not matches[matches.high].contig.insert(o, match=matches[matches.high]):
    stderr.write_line("[scrutinize] error adding new contig")
    quit(2)
  return true

proc insert*(ctgs:var Contigs, dna:string, min_overlap:int=40, max_mismatch:int=0, p_overlap:float64=0.7): bool =
  ## insert a dna sequence into the best contig or create a new one as needed.
  ## return value indicates that it was added to an existing contig
  var o = make_contig(dna)
  return ctgs.insert(o, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap)

proc combine*(contigs: var Contigs, p_overlap:float64=0.6): seq[Contig] =
  ## merge contigs. note that this modifies the contigs in-place.
  contigs.sort(proc (a, b: Contig): int = return b.len - a.len)
  var used = initSet[int]()
  # TODO: make this use count_matches
  for i in 0..contigs.high:
    for j in 0..<i:
      if contigs[i].insert(contigs[j], p_overlap=p_overlap):
        used.incl(j)
  if used.len == 0:
    return contigs

  var k = 0
  for i in 0..contigs.high:
    if i in used: continue
    contigs[k] = contigs[i]
    k += 1
  contigs = contigs[0..<k]
  contigs.sort(proc (a, b: Contig): int = return b.len - a.len)
  return contigs

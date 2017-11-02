import strutils
import sequtils
import algorithm
import sets

type 
  Contig* = ref object of RootObj
    ## a contig is a collection of sequences that have been
    ## merged. the base_count indicates the number of reads supporting
    ## the contig at each base.
    sequence*: string
    base_count*: seq[uint32]
    mismatch_count*: seq[uint32]
    nreads*: int
    # TODO: remove this. just for debugging
    reads: seq[string]

const complement = ['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']

proc revcomp(c: var Contig) =
  c.base_count.reverse
  c.sequence.reverse
  c.mismatch_count.reverse
  for i, ch in c.sequence:
    c.sequence[i] = complement[ch.int]

var DEBUG* = false

proc len*(c:Contig): int {.inline.} =
  ## the length of the sequence of the contig
  return c.sequence.len

proc high*(c:Contig): int {.inline.} =
  return c.sequence.len - 1

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

proc add*(c: Contig, dna: var Contig, min_overlap:int=40, max_mismatch:int=3, p_overlap:float64=0.7, reverse:bool=false): bool =
  ## add a dna sequence to the contig. return indicates if it was added.
  if c.sequence == nil:
    c.sequence = dna.sequence
    c.nreads = dna.nreads
    c.base_count = dna.base_count
    c.mismatch_count = dna.mismatch_count
    when defined(reads):
      c.reads = dna.reads
    return true

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
    if reverse:
      revcomp(dna)
      return c.add(dna, min_overlap=min_overlap, max_mismatch=max_mismatch, p_overlap=p_overlap, reverse=false)

    return false

  if DEBUG:
    var x = ""
    for i in 0..<abs(max_ma_offset):
      x &= " "
    echo "BEFORE"
    if max_ma_offset < 0:
      echo "CTG :", x & c.sequence
      echo "READ:", dna.sequence
    else:
      echo "CTG :", c.sequence
      echo "READ:", x & dna.sequence


  if max_ma_offset >= 0:
    # add sequence ot the right end (or middle).
    var extra = 0
    if (max_ma_offset + dna.len) > c.sequence.len:
      # the added sequence extends the contig 
      extra = max_ma_offset + dna.len - c.sequence.len
      c.sequence &= dna.sequence[(dna.len - extra)..<dna.len]
    if c.sequence.len > c.base_count.len:
      # this should always be true
      var L = c.base_count.len
      c.base_count.set_len(c.sequence.len)
      c.mismatch_count.set_len(c.sequence.len)
      for i in L..<c.base_count.len:
        c.base_count[i] = 1

    #for i in max_ma_offset..<min(c.len, dna.len - max_ma_offset):
    for i in max_ma_offset..<min(c.len-extra, dna.len + max_ma_offset):
      if c.sequence[i] == dna.sequence[i - max_ma_offset]:
        #echo i, " adding ", c.len, " ", i - max_ma_offset, " in ", dna.len
        c.base_count[i] += dna.base_count[i - max_ma_offset]
      elif dna.base_count[i - max_ma_offset] > c.base_count[i]:
        # if the incoming contig has more evidence for a different base, we adopt it.
        c.sequence[i] = dna.sequence[i - max_ma_offset]
        c.base_count[i] = dna.base_count[i - max_ma_offset]
  else:
    # add new sequence to the left end
    c.sequence = dna.sequence[0..<(-max_ma_offset)] & c.sequence
    # if the dna.sequence completely encompasses the contig sequence, we have to added
    # to the right end as well:
    if dna.sequence.len > c.sequence.len:
      var d = dna.sequence.len - c.sequence.len
      c.sequence &= dna.sequence[dna.sequence.len-d..<dna.sequence.len]

    if c.sequence.len > c.base_count.len:
      c.base_count.set_len(c.sequence.len)
      c.mismatch_count.set_len(c.sequence.len)
      # move all the counts to the left
      for i in countdown(c.len-1, -max_ma_offset):
        c.base_count[i] = c.base_count[i + max_ma_offset]
        c.mismatch_count[i] = c.mismatch_count[i + max_ma_offset]
      # zero out the newest.
      for i in countdown(-max_ma_offset-1, 0):
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
  for i, s in max_sub:
      if i + max(0, max_ma_offset) == c.sequence.len: break
      if s != c.sequence[i + max(0, max_ma_offset)]:
        c.mismatch_count[i + max(0, max_ma_offset)] += 1
        found = true
  if found:
    for i, s in max_sub:
      var cidx = i + max(0, max_ma_offset)
      if cidx == c.sequence.len: break
      if c.mismatch_count[cidx] > c.base_count[cidx]:
        # switch to more likely base
        c.base_count[cidx] += c.mismatch_count[cidx]
        c.mismatch_count[cidx] = 0
        c.sequence[cidx] = s

  # end ERROR correction

  c.nreads += dna.nreads
  when defined(reads):
    if dna.reads != nil and c.reads == nil:
      c.reads = dna.reads
    else:
      for r in dna.reads:
        c.reads.add(r)

  if DEBUG:
    echo "AFTER"
    echo c
  return true

proc add*(c: Contig, dna: string, min_overlap:int=40, max_mismatch:int=3, p_overlap:float64=0.8, read:string=nil): bool {.inline.} =
  ## add a sequence to a contig.
  var bc = new_seq[uint32](dna.len)
  for i in 0..bc.high:
    bc[i] = 1
  var o: Contig
  if read == nil:
    o = Contig(sequence: dna, base_count:bc, nreads: 1)
  else:
    o = Contig(sequence: dna, base_count:bc, nreads: 1, reads: @[read])
  o.mismatch_count = new_seq[uint32](dna.len)
  revcomp(o)
  return c.add(o, min_overlap, max_mismatch, p_overlap, reverse=false)

proc combine*(contigs: var seq[Contig], p_overlap:float64=0.6, reverse:bool=false): seq[Contig] =
  ## merge contigs. note that this modifies the contigs in-place.
  contigs.sort(proc (a, b: Contig): int = return b.len - a.len)
  var used = initSet[int]()
  for i in 0..contigs.high:
    for j in 0..<i:
      if contigs[i].add(contigs[j], p_overlap=p_overlap, reverse=reverse):
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

when isMainModule:
  DEBUG = true
  var c: Contig

  if true:
    var dna =               "AAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC"
    var dna_left = "CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTAT"
    var dna_mid =                     "GCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCG"
    # TODO: check counts

    var a = dna[0..<75]
    var b = dna[(dna.len-75)..<dna.len]

    c = Contig()
    echo "\nA:\n##"
    assert c.add(a, p_overlap=0.6)
    echo "\nB:\n##"
    assert c.add(b, p_overlap=0.6)
    assert c.len == dna.len

    echo "\nC: add to left end\n##"
    assert c.add(dna_left, p_overlap=0.6)
    var sl = c.len()
    assert sl > dna.len
    echo "\nD: add in middle\n##"
    assert c.add(dna_mid, p_overlap=0.6)
    assert sl == c.len(), "should not have longer "

    assert c.sequence == "CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC"
    assert join(map(c.base_count, proc(x:uint32): string = $x)) == "1111111112222222222333333333333333444444444444444444444444333333333333333333333333332222222222222222221111111"

    c = Contig()
    assert c.add(         "AAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC")
    assert c.add("CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTAT", p_overlap=0.6)

    c = Contig()
    assert c.add( "AAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTA")
    assert c.nreads == 1
    assert c.add("AAAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTAC")
    assert c.nreads == 2
    assert c.sequence == "AAAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTAC"
    var bc = join(map(c.base_count, proc(x:uint32): string = $x))
    assert bc ==         "122222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222221"

    c = Contig()
    assert c.add(        "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
    assert c.add(          "ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACGG")
    assert c.sequence == "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACGG"
    bc = join(map(c.base_count, proc(x:uint32): string = $x))
    assert bc ==         "112221222222222222222222222222222222222222222222222222222222222222222222222111"

    c = Contig()
    assert c.add(        "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
    assert c.add(          "ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC")
    assert c.sequence == "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC"
    bc = join(map(c.base_count, proc(x:uint32): string = $x))
    assert bc ==         "1122212222222222222222222222222222222222222222222222222222222222222222222221"

    c = Contig()
    assert c.add(        "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGG")
    assert c.add(                                                 "TTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG")
    assert c.sequence == "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGGG"
    assert c.add(                                                 "TTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG")
    # after we accumulate enough mismatches at a given spot, we flip it to a base that matches.
    assert c.sequence == "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG"
    echo c.sequence

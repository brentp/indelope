
type 
  Contig* = ref object of RootObj
    ## a contig is a collection of sequences that have been
    ## merged. the base_count indicates the number of reads supporting
    ## the contig at each base.
    sequence*: string
    ## number of reads covering each base in the contig.
    support*: seq[uint32]
    nreads*: int
    start: int

  correction_site = tuple[qoff:int, toff:int, qbest:bool]
  ## a Match is a putative match. corrections are sites that were mismatched, but that had low
  ## support in either target or query indicating that they could be overwritten.
  Match = tuple[matches: int, offset: int, mismatches: int, corrections:seq[correction_site]]

const unaligned* = low(int)

proc len*(c:Contig): int {.inline.} =
  return c.sequence.len

proc `[]`*(c:Contig, i:int): char {.inline.} =
  return c.sequence[i]

proc allowable_mismatch(qsup: uint32, tsup: uint32): bool {.inline.} =
  ## default implementation to determine if we can overwrite a mismatch
  return ((qsup < 3'u32 and tsup > 3'u32 * qsup) or 
         (tsup < 3'u32 and qsup > 3'u32 * tsup))

proc slide_align(q: Contig, t:Contig, min_overlap:int=50, max_mismatch:int=0): Match =
  ## align (q)uery to (t)arget requiring a number of bases of overlap and fewer
  ## than the specified mismatches. If unaligned, the constant unaligned is returned.
  ## a negative value indicates that the query extends the target to the left.
  ## a positive value indicates the start position of the query into the reference
  ##
  ## default rule is to allow a mismatch in q if qsupport < 3 and tsupport > 5 * qsupport
  ##
  var omin = -(len(q) - min_overlap)
  var omax = len(t) - min_overlap
  var obest = unaligned
  var best_ma = min_overlap - 1
  var best_mm = max_mismatch + 1
  var best_correction: seq[correction_site]
  var correction = new_seq[correction_site](6)
  for o in 0..omax:
    correction.set_len(0)
    var qo = 0
    var to = o
    var mm = 0
    var ma = 0
    while qo < len(q) and to < len(t):
      if q[qo] != t[to]:
        if not allowable_mismatch(q.support[qo], t.support[to]):
          mm += 1
          if mm == max_mismatch:
            break
        else:
          correction.add((qo, to, q.support[qo] > t.support[to]))

      else:
        ma += 1

      qo += 1
      to += 1

    if mm <= max_mismatch and (ma > best_ma or ma == best_ma and mm < best_mm):
      best_ma = ma
      best_mm = mm
      obest = o
      best_correction = correction

  # skip zero as we did that above.
  for o in 1..abs(omin):
    correction.set_len(0)
    var qo = o
    var to = 0
    var mm = 0
    var ma = 0

    while qo < len(q) and to < len(t):
      if q[qo] != t[to]:
        if not allowable_mismatch(q.support[qo], t.support[to]):
          mm += 1
          if mm == max_mismatch:
            break
        else:
          correction.add((qo, to, q.support[qo] > t.support[to]))
      else:
        ma += 1

      qo += 1
      to += 1

    if mm <= max_mismatch and (ma > best_ma or ma == best_ma and mm < best_mm):
      best_ma = ma
      best_mm = mm
      obest = -o
      best_correction = correction

  return (matches:best_ma, offset:obest, mismatches:best_mm, corrections:best_correction)

proc make_contig(dna: string, start:int, support:uint32=1): Contig =
  ## make a contig from a sequence string.
  ## support is only used for testing, should be left as 1
  var bc = new_seq[uint32](dna.len)
  for i in 0..bc.high:
    bc[i] = support
  var o = Contig(sequence: dna, support:bc, nreads: int(support), start:start)
  return o

proc slide_align(q: string, t:string, qstart:int=0, tstart:int=0, min_overlap:int=50, max_mismatch:int=0): Match =
  return slide_align(make_contig(q, qstart), make_contig(t, qstart), min_overlap, max_mismatch)

when isMainModule:
  import unittest

  suite "contig":
    test "slide_align positive":
      var sa = slide_align(     "ACTGGGTACGGT",
                             "TTAACTGGGTACGGT", min_overlap=5)
      check 3 == sa.offset
      check 12 == sa.matches

    test "slide_align extend":
      check 3 == slide_align(   "ACTGGGTACGGTGGG",
                             "TTAACTGGGTACGGT", min_overlap=5).offset
    test "slide_align inside":
      check 3 == slide_align(   "ACTGGGTACG",
                             "TTAACTGGGTACGGT", min_overlap=5).offset

    test "slide_align same":
      check 0 == slide_align("TTAACTGGGTACGGT",
                             "TTAACTGGGTACGGT", min_overlap=5).offset
    test "slide_align left":
      check -1 == slide_align("ATTAACTGGGTACGGT",
                               "TTAACTGGGTACGGT", min_overlap=5).offset

      check -1 == slide_align("ATTAACTGGGTACGGT",
                               "TTAACTGGGTACGGTTTT", min_overlap=5).offset
    test "slide query contains target":
      check -1 == slide_align("ATTAACTGGGTACGGTTTGGGG",
                               "TTAACTGGGTACGGTTTG", min_overlap=5).offset

    test "min_overlap":
      check unaligned == slide_align("ATTAACTGGGTACGGTTTGGGG",
                                      "TTAACTGGGTACGGTTTG", min_overlap=50).offset

    test "corrections":
      var t = make_contig("ATTAACTGGGTACGGTTTGGGG", 0, 2)
      var q = make_contig( "TTAACTGGGXACGGTTTGG", 0, 6)

      var ma = q.slide_align(t, min_overlap=5)
      check ma.corrections == nil

      q = make_contig( "TTAACTGGGXACGGTTTGG", 0, 7)
      ma = q.slide_align(t, min_overlap=5)
      check ma.corrections.len == 1
      check q.sequence[ma.corrections[0].qoff] == 'X'
      check t.sequence[ma.corrections[0].toff] == 'T'
      check ma.corrections[0].qbest

      t = make_contig(    "ATTAACTGGGAACGGTTTGGGG", 0, 7)
      q = make_contig("GGAGATTAACTGGGXACGGTTTGG", 0, 2)
      ma = q.slide_align(t, min_overlap=5)
      check ma.corrections.len == 1
      check q.sequence[ma.corrections[0].qoff] == 'X'
      check t.sequence[ma.corrections[0].toff] == 'A'
      check (not ma.corrections[0].qbest)

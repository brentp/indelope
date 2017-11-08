
type 
  Contig* = ref object of RootObj
    ## a contig is a collection of sequences that have been
    ## merged. the base_count indicates the number of reads supporting
    ## the contig at each base.
    sequence*: string
    base_count*: seq[uint32]
    mismatch_count*: seq[uint32]
    nreads*: int
    start: int

  Match = tuple[matches: int, offset: int, mismatches: int]

const unaligned* = low(int)

proc slide_align(q: string, t:string, min_overlap:int=50, max_mismatch:int=0): Match =
  ## align (q)uery to (t)arget requiring a number of bases of overlap and fewer
  ## than the specified mismatches. If unaligned, the constant unaligned is returned.
  ## a negative value indicates that the query extends the target to the left.
  ## a positive value indicates the start position of the query into the reference
  var omin = -(len(q) - min_overlap)
  var omax = len(t) - min_overlap
  var obest = unaligned
  var best_ma = min_overlap - 1
  var best_mm = max_mismatch + 1
  for o in 0..omax:
    var qo = 0
    var to = o
    var mm = 0
    var ma = 0
    while qo < len(q) and to < len(t):
      if q[qo] != t[to]:
        mm += 1
        if mm == max_mismatch:
          break
      else:
        ma += 1

      qo += 1
      to += 1

    if mm <= max_mismatch and (ma > best_ma or ma == best_ma and mm < best_mm):
      best_ma = ma
      best_mm = mm
      obest = o

  # skip zero as we did that above.
  for o in 1..abs(omin):
    var qo = o
    var to = 0
    var mm = 0
    var ma = 0

    while qo < len(q) and to < len(t):
      if q[qo] != t[to]:
        mm += 1
        if mm == max_mismatch:
          break
      else:
        ma += 1

      qo += 1
      to += 1

    if mm <= max_mismatch and (ma > best_ma or ma == best_ma and mm < best_mm):
      best_ma = ma
      best_mm = mm
      obest = -o

  return (best_ma, obest, best_mm)

#proc insert*(c:Contig, o:Contig, m:Match) =



when isMainModule:
  import unittest

  suite "contig":
    test "slide_align positive":
      check 3 == slide_align(   "ACTGGGTACGGT",
                             "TTAACTGGGTACGGT", min_overlap=5).offset

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


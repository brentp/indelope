import ksw2_c

type 
  Ez* = ref object of RootObj
    c: ptr ksw_extz_t

  CArray{.unchecked.}[T] = array[0..0, T]
  cigar_pair* = tuple[op:uint32, length: uint32]

iterator full_cigar*(e:Ez): cigar_pair =
  var cig = cast[ptr CArray[uint32]](e.c.cigar)
  for i in 0..<e.c.n_cigar:
    yield (cig[i] and 0xf, cig[i] shr 4)

iterator cigar*(e:Ez): cigar_pair =
  ## report the cigar up to the max score for the query.
  var cig = cast[ptr CArray[uint32]](e.c.cigar)
  var result:cigar_pair
  var max_off = uint32(e.c.max_q)
  var off = 0'u32
  for i in 0..<e.c.n_cigar:
    if off >= max_off: break
    result = (cig[i] and 0xf, cig[i] shr 4)
    if result.op != 2:
      off += result.length
    yield result

proc str*(p:cigar_pair): char {.inline.} =
  "MID"[p.op.int]

proc qstop(e: Ez): int {.inline.} =
  ## 1-based stop of the query alignment
  return e.c.max_q + 1

proc tstop(e: Ez): int {.inline.} =
  ## 1-based stop of the target alignment
  return e.c.max_t + 1

proc max_event_length(e: Ez): uint32 {.inline.} =
  result = 0
  for c in e.cigar:
    if c.op != 0:
      result = max(result, c.length)

proc score*(e:Ez): int {.inline.} =
  ## report the max score at the end of the query
  return e.c.mqe_t

proc draw*(e:Ez, q:string, t:string): string =
  ## string represention of the alignment target\nquery
  var qo = ""
  var to = ""
  var qoff = 0
  var toff = 0
  for c in e.cigar:
    if c.op == 0: # M
      qo &= q[qoff..<(qoff + c.length.int)]
      to &= t[toff..<(toff + c.length.int)]
      toff += c.length.int
      qoff += c.length.int
    elif c.op == 1: # I
      qo &= q[qoff..<(qoff + c.length.int)]
      var s = new_string(c.length.int)
      for i, c in s:
        s[i] = ' '
      to &= s
      qoff += c.length.int
    else: # D
      to &= t[toff..<(toff + c.length.int)]
      var s = new_string(c.length.int)
      for i, c in s:
        s[i] = ' '
      qo &= s
      toff += c.length.int
  return to & "\n" & qo



const lookup = [4'u8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

proc encode(dna:string, enc: var seq[uint8]) =
  enc.set_len(dna.len)
  for i, l in dna:
    enc[i] = lookup[l.int]


proc matrix*(match:int8=1, mismatch:int8=(-2)): seq[int8] =
  result = @[match,    mismatch, mismatch, mismatch, 0,
             mismatch, match,    mismatch, mismatch, 0,
             mismatch, mismatch, match,    mismatch, 0,
             mismatch, mismatch, mismatch,    match, 0,
             0,        0,        0,           0,     0]

proc align_to*(query: var seq[uint8], target: var seq[uint8], matrix: var seq[int8], ez:Ez, gap_open:int8=(3), gap_ext:int8=(1), flag:cint=KSW_EZ_RIGHT) =
  #ez.reset()
  var bw = -1 # TODO
  var z = -1
  assert gap_open > 0 and gap_ext > 0
  ksw_extz2_sse(nil.pointer, query.len.cint, cast[ptr uint8](query[0].addr),
                target.len.cint, cast[ptr uint8](target[0].addr),
                5.int8, cast[ptr int8](matrix[0].addr),
                gap_open, gap_ext, bw.cint, z.cint, flag, ez.c)

when isMainModule:

  import unittest
  import sequtils

  var tgt = "CGAAACTGGGCTACTCCATGACCAGGGGCAAAATAGGCTTTTAGCCGCTGCGTTCTGGGAGCTCCTCCCCCTTCTGGGAGCTCCTCCCCCTCCCCAGAAGGCCAAGGGATGTGGGGGCTGGGGGACTGGGAGGCCTGGCAGTCTT"
  var qry = "CGAAACTGGGCTACTCCATGACCAGGGGCAAAATAGGCTTTTAGCCGCTGCGTTCTGGGAGCTCCTCCCCCTCCCCAGAAGGCCAAGGGATGTTGGGG"
  var tenc = new_seq[uint8](tgt.len)
  var tqry = new_seq[uint8](qry.len)
  tgt.encode(tenc)
  qry.encode(tqry)
  var mat = matrix()
  var cez = ksw_extz_t()
  var ez = Ez(c: cez.addr)
  tqry.align_to(tenc, mat, ez)

  var cig = toSeq(ez.cigar)

  suite "ksw2 suite":
    test "encode":
      check tenc[0] == 1
      check tqry[0] == 1

      check tqry.len == qry.len

    test "matrix":
      check mat == @[1'i8, -2, -2, -2, 0, -2, 1, -2, -2, 0, -2, -2, 1, -2, 0, -2, -2, -2, 1, 0, 0, 0, 0, 0, 0] 


    test "cigar":
      check cig[0].str == 'M'
      check cig[0].length == 72

      check cig[1].str == 'D'
      check cig[1].length == 19

      check cig[2].str == 'M'
      check cig[2].length == 26

      check cig.len == 3
     
    test "ends":
      check ez.qstop == 98
      check ez.tstop == 117

    test "score":
      check ez.score == 116

    test "event":
      check ez.max_event_length == 19

  for cig in ez.full_cigar:
    echo cig, cig.str
  echo "base cigar"
  for cig in ez.cigar:
    echo cig, cig.str

  echo "score:", ez.score
  echo "query end :", ez.qstop
  echo "target end:", ez.tstop

  echo ez.draw(qry, tgt)


{.compile: "csrc/ksw2_extz2_sse.c".}
const
  HAVE_KALLOC* = true

const
  KSW_NEG_INF* = - 0x40000000
  KSW_EZ_SCORE_ONLY* = 0x00000001
  KSW_EZ_RIGHT* = 0x00000002
  KSW_EZ_GENERIC_SC* = 0x00000004
  KSW_EZ_APPROX_MAX* = 0x00000008
  KSW_EZ_APPROX_DROP* = 0x00000010
  KSW_EZ_EXTZ_ONLY* = 0x00000040
  KSW_EZ_REV_CIGAR* = 0x00000080
  KSW_EZ_SPLICE_FOR* = 0x00000100
  KSW_EZ_SPLICE_REV* = 0x00000200

type
  ksw_extz_t* {.bycopy.} = object
    max* {.bitsize: 31.}: uint32
    zdropped* {.bitsize: 1.}: uint32
    max_q*: cint
    max_t*: cint               ##  max extension coordinate
    mqe*: cint
    mqe_t*: cint               ##  max score when reaching the end of query
    mte*: cint
    mte_q*: cint               ##  max score when reaching the end of target
    score*: cint               ##  max score reaching both ends; may be KSW_NEG_INF
    m_cigar*: cint
    n_cigar*: cint
    cigar*: ptr uint32


proc ksw_reset_extz*(ez: ptr ksw_extz_t) {.cdecl, importc: "ksw_reset_extz".}
## *
##  NW-like extension
## 
##  @param km        memory pool, when used with kalloc
##  @param qlen      query length
##  @param query     query sequence with 0 <= query[i] < m
##  @param tlen      target length
##  @param target    target sequence with 0 <= target[i] < m
##  @param m         number of residue types
##  @param mat       m*m scoring mattrix in one-dimension array
##  @param gapo      gap open penalty; a gap of length l cost "-(gapo+l*gape)"
##  @param gape      gap extension penalty
##  @param w         band width (<0 to disable)
##  @param zdrop     off-diagonal drop-off to stop extension (positive; <0 to disable)
##  @param flag      flag (see KSW_EZ_* macros)
##  @param ez        (out) scores and cigar
## 
## void ksw_extz(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int flag, ksw_extz_t *ez);

proc ksw_extz2_sse*(km: pointer; qlen: cint; query: ptr uint8; tlen: cint;
                   target: ptr uint8; m: int8; mat: ptr int8; q: int8; e: int8; w: cint;
                   zdrop: cint; flag: cint; ez: ptr ksw_extz_t) {.cdecl, importc.}

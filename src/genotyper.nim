import math
import strutils

const log2 = ln(2'f64)

type
  GT* {.pure.} = enum
    # the genotype of a (presumably diploid sample)
    HOM_REF
    HET
    HOM_ALT
    UNKNOWN

  Genotype* = tuple[GT:GT, GL: array[3, float64]]

const gl_precision = 4
const gt_encodings = ["0/0", "0/1", "1/1", "./."]

proc `$`*(g:GT): string =
  return gt_encodings[g.int]

proc `$`*(g:Genotype): string =
  return $g.GT & ":" & formatFloat(g.GL[0], format=ffDecimal, precision=gl_precision) & "," & formatFloat(g.GL[1], format=ffDecimal, precision=gl_precision) & "," & formatFloat(g.GL[2], format=ffDecimal, precision=gl_precision)

proc qual*(g:Genotype): float64 =
  if g.GT == GT.HOM_REF:
    return g.GL[0] - max(g.GL[1], g.GL[2])
  if g.GT == GT.HET:
    return g.GL[1] - max(g.GL[0], g.GL[2])
  if g.GT == GT.HOM_ALT:
    return g.GL[2] - max(g.GL[0], g.GL[1])
  return 0

proc genotype*(r:int, a:int, error:float64): Genotype =
  var total = float64(r + a)
  var gls: array[3, float64]
  ## eqn2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/
  if total == 0:
    return (GT:GT.UNKNOWN, GL: gls)

  result = (GT: GT(0), GL: gls)
  for G in 0..2:
    result.GL[G] = -total * log2 + r.float64 * ln(G.float64 * error + (2 - G).float64*(1-error)) + a.float64 * ln(G.float64 * (1 - error) + (2 - G).float64 * error)
    if result.GL[G] > result.GL[int(result.GT)]:
      result.GT = GT(G)

when isMainModule:
  var r = genotype(20 - 10, 10, 1e-4)
  assert r.GT == GT.HET
  assert r.GL[1] > r.GL[0]

  r = genotype(20, 0, 1e-4)
  assert r.GT == GT.HOM_REF

  r = genotype(1, 19, 1e-2)
  assert r.GT == GT.HOM_ALT

  r = genotype(1, 19, 1e-8)
  assert r.GT == GT.HET

  r = genotype(0, 0, 1e-8)
  assert r.GT == GT.UNKNOWN

  r = genotype(1, 19, 1e-8)
  assert $r.GT == "0/1"

  echo $r
  echo r.qual

# Package

version       = "0.0.1"
author        = "Brent Pedersen"
description   = "finding missed medium indels by local assembly"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.2", "hts >= 0.1.0", "docopt"
bin = @["indelope"]
srcDir = "src"

skipDirs = @["tests"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"

task docs, "make docs":
  exec "nim doc2 indelope; mkdir -p docs; mv indelope.html docs/index.html"

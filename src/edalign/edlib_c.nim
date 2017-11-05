{.compile: "edlib.cpp".}
##  Status codes
const
  EDLIB_STATUS_OK* = 0
  EDLIB_STATUS_ERROR* = 1

## *
##  Alignment methods - how should Edlib treat gaps before and after query?
## 

type ## *
    ##  Global method. This is the standard method.
    ##  Useful when you want to find out how similar is first sequence to second sequence.
    ## 
  EdlibAlignMode* = enum
    EDLIB_MODE_NW, ## *
                  ##  Prefix method. Similar to global method, but with a small twist - gap at query end is not penalized.
                  ##  What that means is that deleting elements from the end of second sequence is "free"!
                  ##  For example, if we had "AACT" and "AACTGGC", edit distance would be 0, because removing "GGC" from the end
                  ##  of second sequence is "free" and does not count into total edit distance. This method is appropriate
                  ##  when you want to find out how well first sequence fits at the beginning of second sequence.
                  ## 
    EDLIB_MODE_SHW, ## *
                   ##  Infix method. Similar as prefix method, but with one more twist - gaps at query end and start are
                   ##  not penalized. What that means is that deleting elements from the start and end of second sequence is "free"!
                   ##  For example, if we had ACT and CGACTGAC, edit distance would be 0, because removing CG from the start
                   ##  and GAC from the end of second sequence is "free" and does not count into total edit distance.
                   ##  This method is appropriate when you want to find out how well first sequence fits at any part of
                   ##  second sequence.
                   ##  For example, if your second sequence was a long text and your first sequence was a sentence from that text,
                   ##  but slightly scrambled, you could use this method to discover how scrambled it is and where it fits in
                   ##  that text. In bioinformatics, this method is appropriate for aligning read to a sequence.
                   ## 
    EDLIB_MODE_HW


## *
##  Alignment tasks - what do you want Edlib to do?
## 

type
  EdlibAlignTask* = enum
    EDLIB_TASK_DISTANCE,      ## !< Find edit distance and end locations.
    EDLIB_TASK_LOC,           ## !< Find edit distance, end locations and start locations.
    EDLIB_TASK_PATH           ## !< Find edit distance, end locations and start locations and alignment path.


## *
##  Describes cigar format.
##  @see http://samtools.github.io/hts-specs/SAMv1.pdf
##  @see http://drive5.com/usearch/manual/cigar.html
## 

type
  EdlibCigarFormat* = enum
    EDLIB_CIGAR_STANDARD,     ## !< Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
    EDLIB_CIGAR_EXTENDED      ## !< Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.


##  Edit operations.

const
  EDLIB_EDOP_MATCH* = 0
  EDLIB_EDOP_INSERT* = 1
  EDLIB_EDOP_DELETE* = 2
  EDLIB_EDOP_MISMATCH* = 3

## *
##  @brief Defines two given characters as equal.
## 

type
  EdlibEqualityPair* {.bycopy.} = object
    first*: char
    second*: char


## *
##  @brief Configuration object for edlibAlign() function.
## 

type
  EdlibAlignConfig* {.bycopy.} = object
    k*: cint ## *
           ##  Set k to non-negative value to tell edlib that edit distance is not larger than k.
           ##  Smaller k can significantly improve speed of computation.
           ##  If edit distance is larger than k, edlib will set edit distance to -1.
           ##  Set k to negative value and edlib will internally auto-adjust k until score is found.
           ## 
    ## *
    ##  Alignment method.
    ##  EDLIB_MODE_NW: global (Needleman-Wunsch)
    ##  EDLIB_MODE_SHW: prefix. Gap after query is not penalized.
    ##  EDLIB_MODE_HW: infix. Gaps before and after query are not penalized.
    ## 
    mode*: EdlibAlignMode ## *
                        ##  Alignment task - tells Edlib what to calculate. Less to calculate, faster it is.
                        ##  EDLIB_TASK_DISTANCE - find edit distance and end locations of optimal alignment paths in target.
                        ##  EDLIB_TASK_LOC - find edit distance and start and end locations of optimal alignment paths in target.
                        ##  EDLIB_TASK_PATH - find edit distance, alignment path (and start and end locations of it in target).
                        ## 
    task*: EdlibAlignTask ## *
                        ##  List of pairs of characters, where each pair defines two characters as equal.
                        ##  This way you can extend edlib's definition of equality (which is that each character is equal only
                        ##  to itself).
                        ##  This can be useful if you have some wildcard characters that should match multiple other characters,
                        ##  or e.g. if you want edlib to be case insensitive.
                        ##  Can be set to NULL if there are none.
                        ## 
    additionalEqualities*: ptr EdlibEqualityPair ## *
                                              ##  Number of additional equalities, which is non-negative number.
                                              ##  0 if there are none.
                                              ## 
    additionalEqualitiesLength*: cint


## *
##  Helper method for easy construction of configuration object.
##  @return Configuration object filled with given parameters.
## 

proc edlibNewAlignConfig*(k: cint; mode: EdlibAlignMode; task: EdlibAlignTask;
                         additionalEqualities: ptr EdlibEqualityPair;
                         additionalEqualitiesLength: cint): EdlibAlignConfig {.cdecl, importc: "edlibNewAlignConfig".}

proc edlibDefaultAlignConfig*(): EdlibAlignConfig {.cdecl, importc: "edlibDefaultAlignConfig".}

## *
##  Container for results of alignment done by edlibAlign() function.
## 

type
  EdlibAlignResult* {.bycopy.} = object
    status*: cint ## *
                ##  EDLIB_STATUS_OK or EDLIB_STATUS_ERROR. If error, all other fields will have undefined values.
                ## 
    ## *
    ##  -1 if k is non-negative and edit distance is larger than k.
    ## 
    editDistance*: cint ## *
                      ##  Array of zero-based positions in target where optimal alignment paths end.
                      ##  If gap after query is penalized, gap counts as part of query (NW), otherwise not.
                      ##  Set to NULL if edit distance is larger than k.
                      ##  If you do not free whole result object using edlibFreeAlignResult(), do not forget to use free().
                      ## 
    endLocations*: ptr cint ## *
                         ##  Array of zero-based positions in target where optimal alignment paths start,
                         ##  they correspond to endLocations.
                         ##  If gap before query is penalized, gap counts as part of query (NW), otherwise not.
                         ##  Set to NULL if not calculated or if edit distance is larger than k.
                         ##  If you do not free whole result object using edlibFreeAlignResult(), do not forget to use free().
                         ## 
    startLocations*: ptr cint   ## *
                           ##  Number of end (and start) locations.
                           ## 
    numLocations*: cint ## *
                      ##  Alignment is found for first pair of start and end locations.
                      ##  Set to NULL if not calculated.
                      ##  Alignment is sequence of numbers: 0, 1, 2, 3.
                      ##  0 stands for match.
                      ##  1 stands for insertion to target.
                      ##  2 stands for insertion to query.
                      ##  3 stands for mismatch.
                      ##  Alignment aligns query to target from begining of query till end of query.
                      ##  If gaps are not penalized, they are not in alignment.
                      ##  If you do not free whole result object using edlibFreeAlignResult(), do not forget to use free().
                      ## 
    alignment*: ptr cuchar      ## *
                        ##  Length of alignment.
                        ## 
    alignmentLength*: cint ## *
                         ##  Number of different characters in query and target together.
                         ## 
    alphabetLength*: cint


## *
##  Frees memory in EdlibAlignResult that was allocated by edlib.
##  If you do not use it, make sure to free needed members manually using free().
## 

proc edlibFreeAlignResult*(r: EdlibAlignResult) {.cdecl, importc: "edlibFreeAlignResult".}
## *
##  Aligns two sequences (query and target) using edit distance (levenshtein distance).
##  Through config parameter, this function supports different alignment methods (global, prefix, infix),
##  as well as different modes of search (tasks).
##  It always returns edit distance and end locations of optimal alignment in target.
##  It optionally returns start locations of optimal alignment in target and alignment path,
##  if you choose appropriate tasks.
##  @param [in] query  First sequence.
##  @param [in] queryLength  Number of characters in first sequence.
##  @param [in] target  Second sequence.
##  @param [in] targetLength  Number of characters in second sequence.
##  @param [in] config  Additional alignment parameters, like alignment method and wanted results.
##  @return  Result of alignment, which can contain edit distance, start and end locations and alignment path.
##           Make sure to clean up the object using edlibFreeAlignResult() or by manually freeing needed members.
## 

proc edlibAlign*(query: cstring; queryLength: cint; target: cstring;
                targetLength: cint; config: EdlibAlignConfig): EdlibAlignResult {.cdecl, importc: "edlibAlign".}
## *
##  Builds cigar string from given alignment sequence.
##  @param [in] alignment  Alignment sequence.
##      0 stands for match.
##      1 stands for insertion to target.
##      2 stands for insertion to query.
##      3 stands for mismatch.
##  @param [in] alignmentLength
##  @param [in] cigarFormat  Cigar will be returned in specified format.
##  @return Cigar string.
##      I stands for insertion.
##      D stands for deletion.
##      X stands for mismatch. (used only in extended format)
##      = stands for match. (used only in extended format)
##      M stands for (mis)match. (used only in standard format)
##      String is null terminated.
##      Needed memory is allocated and given pointer is set to it.
##      Do not forget to free it later using free()!
## 

proc edlibAlignmentToCigar*(alignment: ptr cuchar; alignmentLength: cint;
                           cigarFormat: EdlibCigarFormat): cstring {.cdecl, importc: "edlibAlignmentToCigar".}

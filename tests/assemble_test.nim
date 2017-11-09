import unittest
include assemble

var dna =               "AAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC"
var dna_left = "CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTAT"
var dna_mid =                     "GCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCG"
var a = dna[0..<75]
var b = dna[(dna.len-75)..<dna.len]

suite "assemble-suite":
  test "test adding":
    # TODO: check counts

    var c = Contig()
    check c.insert(a, p_overlap=0.6)
    check c.insert(b, p_overlap=0.6)
    check c.len == dna.len

  test "add to left end":
    var c = Contig()
    check c.insert(a, p_overlap=0.6)
    check c.insert(b, p_overlap=0.6)
    check c.insert(dna_left, p_overlap=0.6)
    var sl = c.len()
    check sl > dna.len

  test "add to middle":
    var c = Contig()
    check c.insert(a, p_overlap=0.6)
    check c.insert(b, p_overlap=0.6)
    check c.insert(dna_left, p_overlap=0.6)
    var sl = c.len()
    check c.insert(dna_mid, p_overlap=0.6)
    check sl == c.len() #, "should not have longer "

    check c.sequence == "CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC"

    check join(map(c.base_count, proc(x:uint32): string = $x)) == "1111111112222222222333333333333333444444444444444444444444333333333333333333333333332222222222222222221111111"

  test "mismatches":

    var c = Contig()
    check c.insert(         "AAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC")
    check c.insert("CGCGCGCGTAAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTAT", p_overlap=0.6)

  test "count reads":
    var c = Contig()
    check c.insert( "AAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTA")
    check c.nreads == 1
    check c.insert("AAAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTAC")
    check c.nreads == 2
    check c.sequence == "AAAAAAAAAAAAAATTCCTTTGGGAAGGCCTTCTACATAAAAATCTTCAACATGAGACTGGAAAAAAGGGTATGGGATCATCACCGGACCTTTGGCTTTTAC"
    var bc = join(map(c.base_count, proc(x:uint32): string = $x))
    check bc ==         "122222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222221"

  test "count mismatch":

    var c = Contig()
    check c.insert(        "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
    check c.insert(          "ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACGG")
    check c.sequence == "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACGG"
    var bc = join(map(c.base_count, proc(x:uint32): string = $x))
    check bc ==         "112221222222222222222222222222222222222222222222222222222222222222222222222111"

  test "count overhang":
    var c = Contig()
    check c.insert(        "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
    check c.insert(          "ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC")
    check c.sequence == "TAACCATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC"
    var bc = join(map(c.base_count, proc(x:uint32): string = $x))
    check bc ==         "1122212222222222222222222222222222222222222222222222222222222222222222222221"

  test "flip mismatch":

    var c = Contig()
    check c.insert(        "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGG")
    check c.insert(                                                 "TTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG")
    check c.sequence == "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGGG"
    check c.insert(                                                 "TTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG")
    # after we accumulate enough mismatches at a given spot, we flip it to a base that matches.
    check c.sequence == "CTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTG"

  test "contig trimming":

    var c = Contig()
    check c.insert(a)
    check c.insert(b, p_overlap=0.6)

    var s = $c
    check s.split("\n")[0] == "length: 100, reads: 2"
    c.trim()

  test "fastq":
    var c = Contig()
    check c.insert(a)
    check c.insert(b, p_overlap=0.6)
    var s = ""
    var res = c.fastq(s, name="read")
    check res.split("\n")[0].strip() == "@read#nreads=2#len=100"
    check res.split("\n")[1].strip() == "AAAAACTCTAGCTATATATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGAACTAATCATATTTTATATCTTCTTCGAAACCAC"
    check res.split("\n")[2].strip() == "+"
    check res.split("\n")[3].len == c.len

suite "match-suite":

  test "match sorting":

    var aa = Match(matches: 98, mm: 2)
    var bb = Match(matches: 100, mm: 4)
    var mas = @[aa, bb]
    mas.sort(match_sort)
    check mas[0].matches == 98
    aa.matches = 100
    mas.sort(match_sort)
    check mas[0].mm == 4
    check mas[0].matches == 100
    check mas[1].matches == 100


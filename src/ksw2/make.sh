echo "
#ifdef C2NIM
#mangle uint8_t uint8
#mangle uint16_t uint16
#mangle uint64_t uint64
#mangle uint32_t uint32
#mangle int8_t int8
#mangle int32_t int32
#mangle int64_t int64
#mangle ssize_t int64
#endif
" > ksw2_c.h

echo "#define HAVE_KALLOC
" >> ksw2_c.h
cat csrc/ksw2.h >> ksw2_c.h

c2nim --cdecl --stdcall ksw2_c.h
sed -i 's/cigar_\:/cigar\:/g' ksw2_c.nim
sed -i 's/cdecl/cdecl, importc/g' ksw2_c.nim
mv ksw2_c.nim o.nim
echo '{.compile: "csrc/ksw2_extz2_sse.c".}' > ksw2_c.nim
cat o.nim >> ksw2_c.nim
rm o.nim

#pragma once

// .svb -- the 2-bit read format snape reads natively (src/svb_pack.cpp, src/svb_dump.cpp).
//
// A 32-byte "SVB1" header, then one record per read: a u32 tag whose low 31 bits are the read length
// and whose bit 31 is a RAW escape, followed by the payload. Packed payload is 2 bits per base
// (A=00 C=01 G=10 T=11, base i in byte i/4, high 2 bits first); a read containing any non-ACGT base
// (an N) is escaped and stored as plain bytes instead, so packing is lossless. Quality is not stored
// -- snape never reads it.
//
// There is no index and no footer, and svb_next stops at a clean EOF, so concatenating one file whole
// with the others' 32-byte headers stripped yields a valid single .svb. Measured on HiFi reads the
// format is ~8x smaller than the FASTQ it came from.
//
// The reader lives in graph.cpp (open_svb / svb_next); these two are the command-line ends of it.
// Both take their arguments positionally, as `snape pack` / `snape dump`, and return 0 on success.

// FASTQ on stdin -> <out.svb>. Streaming, O(1) memory.
int cmd_pack(int argc, char** argv);

// <in.svb> -> FASTQ on stdout. Header and quality are synthesised; only the sequence is meaningful.
int cmd_dump(int argc, char** argv);

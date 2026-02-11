#include <iostream>
#include <cassert>
#include <string>
#include <vector>

#include "../src/sv.h"

using namespace std;

Bubble makeBubble(const string& normal, const string& tumor) {
    Bubble b; b.start_node = 0; b.end_node = 0;
    b.normal_seq = normal; b.tumor_seq = tumor;
    return b;
}

void testDetectSV() {
    cout << "Testing detectSV\n";

    // indel: single char gap
    SVRegion r1 = detectSV("ABCDE", "ABDE");
    cout << "  ABCDE vs ABDE → " << r1.start << " " << r1.length << " " << r1.call << "\n";
    assert(r1.start == 2 && r1.length == 1 && r1.call == "indel");

    // indel: multi char gap
    SVRegion r2 = detectSV("ABCXYZDE", "ABCDE");
    cout << "  ABCXYZDE vs ABCDE → " << r2.start << " " << r2.length << " " << r2.call << "\n";
    assert(r2.start == 3 && r2.length == 3 && r2.call == "indel");

    // snv: single base difference, same length
    SVRegion r3 = detectSV("ABXDE", "ABYDE");
    cout << "  ABXDE vs ABYDE → " << r3.start << " " << r3.length << " " << r3.call << "\n";
    assert(r3.start == 2 && r3.length == 1 && r3.call == "snv");

    // multi: gaps in both
    SVRegion r4 = detectSV("ABCXDE", "ABYDE");
    cout << "  ABCXDE vs ABYDE → " << r4.call << "\n";
    assert(r4.call == "multi");

    cout << "PASS: detectSV\n\n";
}

void testCallSV() {
    cout << "Testing callSV\n";
    int64_t ref_offset = 1000;

    // deletion: XYZ deleted from normal
    auto c1 = callSV(makeBubble("ABCXYZDE", "ABCDE"), ref_offset);
    cout << "  ABCXYZDE vs ABCDE → " << c1[0].svtype << " " << c1[0].start << " " << c1[0].length << " ref=" << c1[0].normal_seq << " alt=" << c1[0].tumor_seq << " ref_pos=" << c1[0].ref_pos << "\n";
    assert(c1.size() == 1 && c1[0].svtype == "DEL" && c1[0].start == 3 && c1[0].length == 3);
    assert(c1[0].normal_seq == "XYZ" && c1[0].tumor_seq == "");
    assert(c1[0].ref_pos == ref_offset + 3);

    // insertion: XX inserted in tumor
    auto c2 = callSV(makeBubble("AGCTTAG", "AGCTXXTTAG"), ref_offset);
    cout << "  AGCTTAG vs AGCTXXTTAG → " << c2[0].svtype << " " << c2[0].start << " " << c2[0].length << " ref=" << c2[0].normal_seq << " alt=" << c2[0].tumor_seq << " ref_pos=" << c2[0].ref_pos << "\n";
    assert(c2.size() == 1 && c2[0].svtype == "INS" && c2[0].start == 4 && c2[0].length == 2);
    assert(c2[0].normal_seq == "" && c2[0].tumor_seq == "XX");
    assert(c2[0].ref_pos == ref_offset + 4);

    // single base del: T deleted from normal
    auto c3 = callSV(makeBubble("AGTCG", "AGCG"), ref_offset);
    cout << "  AGTCG vs AGCG → " << c3[0].svtype << " " << c3[0].start << " " << c3[0].length << " ref=" << c3[0].normal_seq << " alt=" << c3[0].tumor_seq << " ref_pos=" << c3[0].ref_pos << "\n";
    assert(c3.size() == 1 && c3[0].svtype == "DEL" && c3[0].length == 1);
    assert(c3[0].normal_seq == "T" && c3[0].tumor_seq == "");
    assert(c3[0].ref_pos == ref_offset + 2);

    // multi: 5 variants (3 SNP + 1 DEL + 1 INS) via parasail alignment
    // Variants separated by 6-base matching blocks to force unambiguous alignment.
    //
    // normal: ACGTAC G GATCGA CCC TGACGT A CCTTAG       GCATGC T AGTCAG   (42 bp)
    // tumor:  ACGTAC T GATCGA     TGACGT G CCTTAG AATT  GCATGC C AGTCAG   (43 bp)
    //           6=   1X   6=   3D   6=   1X   6=   4I     6=   1X  6=
    //
    // Expected calls (positions in normal):
    //   SNP start=6  len=1   (G→T)
    //   DEL start=13 len=3   (CCC deleted)
    //   SNP start=22 len=1   (A→G)
    //   INS start=29 len=4   (AATT inserted)
    //   SNP start=35 len=1   (T→C)
    string multi_normal = "ACGTACGGATCGACCCTGACGTACCTTAGGCATGCTAGTCAG";
    string multi_tumor  = "ACGTACTGATCGATGACGTGCCTTAGAATTGCATGCCAGTCAG";
    auto cm = callSV(makeBubble(multi_normal, multi_tumor), ref_offset);
    cout << "  multi (42bp vs 43bp) → " << cm.size() << " calls:\n";
    for (auto& c : cm) cout << "    " << c.svtype << " start=" << c.start << " len=" << c.length << " ref=" << c.normal_seq << " alt=" << c.tumor_seq << " ref_pos=" << c.ref_pos << "\n";

    assert(cm.size() == 5);

    // SNP G→T at pos 6
    assert(cm[0].svtype == "SNP" && cm[0].start == 6  && cm[0].length == 1);
    assert(cm[0].normal_seq == "G" && cm[0].tumor_seq == "T");
    assert(cm[0].ref_pos == ref_offset + 6);

    // DEL CCC at pos 13
    assert(cm[1].svtype == "DEL" && cm[1].start == 13 && cm[1].length == 3);
    assert(cm[1].normal_seq == "CCC" && cm[1].tumor_seq == "");
    assert(cm[1].ref_pos == ref_offset + 13);

    // SNP A→G at pos 22
    assert(cm[2].svtype == "SNP" && cm[2].start == 22 && cm[2].length == 1);
    assert(cm[2].normal_seq == "A" && cm[2].tumor_seq == "G");
    assert(cm[2].ref_pos == ref_offset + 22);

    // INS AATT at pos 29
    assert(cm[3].svtype == "INS" && cm[3].start == 29 && cm[3].length == 4);
    assert(cm[3].normal_seq == "" && cm[3].tumor_seq == "AATT");
    assert(cm[3].ref_pos == ref_offset + 29);

    // SNP T→C at pos 35
    assert(cm[4].svtype == "SNP" && cm[4].start == 35 && cm[4].length == 1);
    assert(cm[4].normal_seq == "T" && cm[4].tumor_seq == "C");
    assert(cm[4].ref_pos == ref_offset + 35);

    cout << "PASS: callSV\n\n";
}

int main() {
    cout << "=== SVCall Tests ===\n\n";
    testDetectSV();
    testCallSV();
    cout << "All tests passed!\n";
    return 0;
}

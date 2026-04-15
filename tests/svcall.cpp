#include <iostream>
#include <cassert>
#include <string>
#include <vector>

#include "../src/sv.h"
#include "../src/utils.h"

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

    // multi: DEL + INS separated by a long unique anchor (>=50bp)
    // normal: [anchor 60bp] CCC [anchor 60bp]          [anchor 60bp]
    // tumor:  [anchor 60bp]     [anchor 60bp] AATTGGCC [anchor 60bp]
    string anc = "ACGTAACCGGTTGATCACGTAACCGGTTGATCACGTAACCGGTTGATCACGTAACCGGTTGA"; // 62bp unique
    string multi_normal = anc + "CCC" + anc + anc;
    string multi_tumor  = anc +         anc + "AATTGGCC" + anc;
    auto cm = callSV(makeBubble(multi_normal, multi_tumor), ref_offset);
    cout << "  multi (" << multi_normal.size() << "bp vs " << multi_tumor.size() << "bp) → " << cm.size() << " calls:\n";
    for (auto& c : cm) cout << "    " << c.svtype << " start=" << c.start << " len=" << c.length << " ref=" << c.normal_seq << " alt=" << c.tumor_seq << " ref_pos=" << c.ref_pos << "\n";

    assert(cm.size() >= 2);
    // At least a DEL and an INS must be found
    bool has_del = false, has_ins = false;
    for (auto& c : cm) {
        if (c.svtype == "DEL") has_del = true;
        if (c.svtype == "INS") has_ins = true;
    }
    assert(has_del && has_ins);

    cout << "PASS: callSV\n\n";
}

void testReverseComplement() {
    cout << "Testing reverse_complement\n";

    assert(reverse_complement("") == "");
    assert(reverse_complement("A") == "T");
    assert(reverse_complement("T") == "A");
    assert(reverse_complement("C") == "G");
    assert(reverse_complement("G") == "C");
    assert(reverse_complement("AAAA") == "TTTT");
    assert(reverse_complement("AACG") == "CGTT");
    assert(reverse_complement("ACGT") == "ACGT");  // palindrome
    assert(reverse_complement("AACGTAA") == "TTACGTT");

    cout << "PASS: reverse_complement\n\n";
}

void testINV() {
    cout << "Testing INV detection in callSV\n";
    int64_t ref_offset = 1000;

    // Pure inversion: tumor_seq == reverse_complement(normal_seq)
    // Use non-palindromic sequence >= MIN_INV_SIZE (10)
    string normal = "AACCGGTTAACC";  // 12bp
    string tumor = reverse_complement(normal);  // GGTTAACCGGTT
    cout << "  normal=" << normal << " tumor=" << tumor << "\n";

    auto c1 = callSV(makeBubble(normal, tumor), ref_offset);
    cout << "  Pure INV: " << c1.size() << " calls, type=" << c1[0].svtype << "\n";
    assert(c1.size() == 1);
    assert(c1[0].svtype == "INV");
    assert(c1[0].length == (int64_t)normal.size());
    assert(c1[0].ref_pos == ref_offset);

    // Noisy inversion: ~92% identity (1 mismatch out of 12)
    string noisy_tumor = tumor;
    noisy_tumor[0] = 'A';  // was G → now A
    auto c2 = callSV(makeBubble(normal, noisy_tumor), ref_offset);
    cout << "  Noisy INV (~92%): " << c2.size() << " calls, type=" << c2[0].svtype << "\n";
    assert(c2.size() == 1);
    assert(c2[0].svtype == "INV");

    // Below INV threshold: 3 mismatches out of 12 = 75% identity → not INV
    string low_tumor = tumor;
    low_tumor[0] = 'A'; low_tumor[1] = 'A'; low_tumor[2] = 'A';
    auto c3 = callSV(makeBubble(normal, low_tumor), ref_offset);
    cout << "  Low identity (~75%): " << c3.size() << " calls, type=" << c3[0].svtype << "\n";
    assert(c3[0].svtype != "INV");

    // Too short: 6bp < MIN_INV_SIZE (10) → not INV
    string short_n = "AACCGG";
    string short_t = reverse_complement(short_n);
    auto c4 = callSV(makeBubble(short_n, short_t), ref_offset);
    cout << "  Short seq (6bp): " << c4.size() << " calls, type=" << c4[0].svtype << "\n";
    assert(c4[0].svtype != "INV");

    // Longer inversion (50bp)
    string long_normal = "ACGTACGTACGTACGTAACCGGTTAACCGGTTACGTACGTACGTACGTAA";  // 50bp
    string long_tumor = reverse_complement(long_normal);
    auto c5 = callSV(makeBubble(long_normal, long_tumor), ref_offset);
    cout << "  Long INV (50bp): " << c5.size() << " calls, type=" << c5[0].svtype << "\n";
    assert(c5.size() == 1);
    assert(c5[0].svtype == "INV");
    assert(c5[0].length == 50);

    cout << "PASS: INV detection\n\n";
}

int main() {
    cout << "=== SVCall Tests ===\n\n";
    testDetectSV();
    testCallSV();
    testReverseComplement();
    testINV();
    cout << "All tests passed!\n";
    return 0;
}

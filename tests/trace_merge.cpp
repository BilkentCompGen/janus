#include <iostream>
#include <string>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/utils.h"
using namespace std;
int main(){
    // Reproduce the exact merge sequence from the bubble
    // Start node sequence (121bp)
    string start_seq = "CACCTCGGAACCTGAAGCAGAGCCGTGCTCAGTCACCAGGGAGCCAGCTCTGAAGCACCCCCTACCCTGAGGTCTGGCTGGGAGTGGTCGTATCCAGGGATGATGAAGGGAATATGGGCCTG";
    // T1 sequence (133bp)
    string t1_seq = "GCCAGCTCTGAAGCACCCCCTACCCTGAGGTCTGGCTGGGAGTGGTCGTATCCAGGGATGATGAAGGGAATATGGGCCTGGCTGTCTGGCTGGGAGTGGTCGTATCCAGGGATGATGAAGGGAATATGGGCCTG";
    // T2 sequence (142bp)
    string t2_seq = "TGGGCCTGGCTGTCTGGCTGGGAGTGGTCGTATCCAGGGATGATGAAGGGAATATGGGCCTGGCTGTCCTAAGGCTGAAGGCGGCTCCTGGTTCAGTCACTGACTTCCCCTCCCCAAGCTTCAGTTTCCTCATCTGGTTTGTGG";
    // End node sequence (145bp)
    string end_seq = "GAAGGGAATATGGGCCTGGCTGTCCTAAGGCTGAAGGCGGCTCCTGGTTCAGTCACTGacttcccctccccaagcttcagtttcctcatctggtttgtgggaatggtccttctgcctgcctgctccctgctcgctctctcagtta";

    size_t start_len = start_seq.size();
    cout << "start_seq len=" << start_len << endl;

    // Step 1: acc = start, merge T1
    Node acc;
    acc.sequence = start_seq;
    Node t1_node; t1_node.sequence = t1_seq;
    acc.merge(t1_node);
    cout << "\nAfter merge(T1):" << endl;
    cout << "  acc len=" << acc.sequence.size() << endl;
    cout << "  acc seq=" << acc.sequence << endl;

    // Step 2: next_acc = acc, merge T2
    Node next_acc;
    next_acc.sequence = acc.sequence;
    Node t2_node; t2_node.sequence = t2_seq;
    next_acc.merge(t2_node);
    cout << "\nAfter merge(T2):" << endl;
    cout << "  next_acc len=" << next_acc.sequence.size() << endl;
    cout << "  next_acc seq=" << next_acc.sequence << endl;

    // Step 3: merge end node
    size_t before = next_acc.sequence.size();
    Node end_node; end_node.sequence = end_seq;
    next_acc.merge(end_node);
    size_t end_overlap = end_node.sequence.size() - (next_acc.sequence.size() - before);
    cout << "\nAfter merge(end):" << endl;
    cout << "  next_acc len=" << next_acc.sequence.size() << endl;
    cout << "  before=" << before << " end_overlap=" << end_overlap << endl;

    string tumor_seq = next_acc.sequence.substr(start_len, before - start_len - end_overlap);
    cout << "\nExtracted tumor_seq:" << endl;
    cout << "  len=" << tumor_seq.size() << endl;
    cout << "  seq=" << tumor_seq << endl;

    // Also show: what does merge find as overlap at each step?
    cout << "\n=== OVERLAP ANALYSIS ===" << endl;

    // Start vs T1: what suffix of start matches prefix of T1?
    for (size_t len = min(start_seq.size(), t1_seq.size()); len > 0; len--) {
        if (start_seq.compare(start_seq.size()-len, len, t1_seq, 0, len) == 0) {
            cout << "start→T1 overlap: " << len << "bp" << endl;
            cout << "  overlap: " << start_seq.substr(start_seq.size()-len) << endl;
            break;
        }
    }

    // After merge(T1), what suffix of acc matches prefix of T2?
    string acc_after_t1 = start_seq;
    { Node a; a.sequence=acc_after_t1; Node b; b.sequence=t1_seq; a.merge(b); acc_after_t1=a.sequence; }
    for (size_t len = min(acc_after_t1.size(), t2_seq.size()); len > 0; len--) {
        if (acc_after_t1.compare(acc_after_t1.size()-len, len, t2_seq, 0, len) == 0) {
            cout << "acc(after T1)→T2 overlap: " << len << "bp" << endl;
            cout << "  overlap: " << acc_after_t1.substr(acc_after_t1.size()-len, min(len,(size_t)80)) << (len>80?"...":"") << endl;
            break;
        }
    }

    // After merge(T2), what suffix matches prefix of end?
    string acc_after_t2 = acc_after_t1;
    { Node a; a.sequence=acc_after_t2; Node b; b.sequence=t2_seq; a.merge(b); acc_after_t2=a.sequence; }
    for (size_t len = min(acc_after_t2.size(), end_seq.size()); len > 0; len--) {
        if (acc_after_t2.compare(acc_after_t2.size()-len, len, end_seq, 0, len) == 0) {
            cout << "acc(after T2)→end overlap: " << len << "bp" << endl;
            cout << "  overlap: " << acc_after_t2.substr(acc_after_t2.size()-len, min(len,(size_t)80)) << (len>80?"...":"") << endl;
            break;
        }
    }

    return 0;
}

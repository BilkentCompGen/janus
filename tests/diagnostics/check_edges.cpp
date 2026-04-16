#include <iostream>
#include <string>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/utils.h"
using namespace std;
string colorStr(uint32_t c) {
    bool n=(c&0x1),t=(c&0x2);
    if(n&&t)return"N+T";if(n)return"N";if(t)return"T";return"?";
}
int main(){
    LCP_INIT();
    LCPGraph graph;
    string ref_seq,tumor_seq;
    {FastaStreamer s("simulation/ref.fasta",false);
     for(const FastaPair&p:s){
        string c=p.first;size_t hp=c.find("_hap");if(hp!=string::npos)c=c.substr(0,hp);
        if(c!="19")continue;
        ref_seq=p.second;uint32_t sid=graph.seq_names.size();graph.seq_names.push_back(p.first);
        struct lps l;init_lps(&l,ref_seq.c_str(),ref_seq.size());lps_deepen(&l,5);
        graph.addString(&l,0,ref_seq.c_str(),sid);free_lps(&l);
    }}
    {FastaStreamer s("simulation/duplication_simulation.fa",false);
     for(const FastaPair&p:s){
        string c=p.first;size_t hp=c.find("_hap");if(hp!=string::npos)c=c.substr(0,hp);
        if(c!="19")continue;
        tumor_seq=p.second;uint32_t sid=graph.seq_names.size();graph.seq_names.push_back(p.first);
        struct lps l;init_lps(&l,tumor_seq.c_str(),tumor_seq.size());lps_deepen(&l,5);
        graph.addString(&l,1,tumor_seq.c_str(),sid);free_lps(&l);
    }}

    uint64_t t1 = 7856366017321890231ULL;   // first T-only core (tumor 616685)
    uint64_t t2 = 15783397869307142715ULL;   // second T-only core (tumor 616686)
    uint64_t end_node = 17856515116681298403ULL; // end node (N+T)
    uint64_t start_node = 4996526221544266035ULL; // start node (N+T)
    uint64_t n_only = 13404391996925811288ULL; // N-only core

    cout << "=== EDGE CHECK ===" << endl;

    // t1 successors
    auto it = graph.next_edges.find(t1);
    if(it!=graph.next_edges.end()){
        cout << "T1 successors:" << endl;
        for(auto nh:it->second){
            auto nit=graph.nodes.find(nh);
            cout << "  hash=" << nh << " color=" << colorStr(nit->second.colors)
                 << " seq_len=" << nit->second.sequence.size() << endl;
            if(nh==t2) cout << "  ^^^ THIS IS T2!" << endl;
            if(nh==end_node) cout << "  ^^^ THIS IS END NODE!" << endl;
        }
    } else cout << "T1 has NO successors" << endl;

    // t2 successors
    it = graph.next_edges.find(t2);
    if(it!=graph.next_edges.end()){
        cout << "T2 successors:" << endl;
        for(auto nh:it->second){
            auto nit=graph.nodes.find(nh);
            cout << "  hash=" << nh << " color=" << colorStr(nit->second.colors)
                 << " seq_len=" << nit->second.sequence.size() << endl;
            if(nh==end_node) cout << "  ^^^ THIS IS END NODE!" << endl;
        }
    } else cout << "T2 has NO successors" << endl;

    // t1 predecessors
    it = graph.prev_edges.find(t1);
    if(it!=graph.prev_edges.end()){
        cout << "T1 predecessors:" << endl;
        for(auto nh:it->second){
            auto nit=graph.nodes.find(nh);
            cout << "  hash=" << nh << " color=" << colorStr(nit->second.colors) << endl;
            if(nh==start_node) cout << "  ^^^ THIS IS START NODE!" << endl;
        }
    }

    // t2 predecessors
    it = graph.prev_edges.find(t2);
    if(it!=graph.prev_edges.end()){
        cout << "T2 predecessors:" << endl;
        for(auto nh:it->second){
            auto nit=graph.nodes.find(nh);
            cout << "  hash=" << nh << " color=" << colorStr(nit->second.colors) << endl;
            if(nh==t1) cout << "  ^^^ THIS IS T1!" << endl;
            if(nh==start_node) cout << "  ^^^ THIS IS START NODE!" << endl;
        }
    }

    return 0;
}

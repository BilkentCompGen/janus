#include <iostream>
#include <string>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/utils.h"
using namespace std;
int main(){
    LCP_INIT();
    LCPGraph graph;
    string ref_seq;
    {FastaStreamer s("simulation/ref.fasta",false);
     for(const FastaPair&p:s){
        string c=p.first;size_t hp=c.find("_hap");if(hp!=string::npos)c=c.substr(0,hp);
        if(c!="19")continue;ref_seq=p.second;
        uint32_t sid=graph.seq_names.size();graph.seq_names.push_back(p.first);
        struct lps l;init_lps(&l,ref_seq.c_str(),ref_seq.size());lps_deepen(&l,5);
        graph.addString(&l,0,ref_seq.c_str(),sid);free_lps(&l);
    }}
    {FastaStreamer s("simulation/duplication_simulation.fa",false);
     for(const FastaPair&p:s){
        string c=p.first;size_t hp=c.find("_hap");if(hp!=string::npos)c=c.substr(0,hp);
        if(c!="19")continue;
        string tumor_seq=p.second;
        uint32_t sid=graph.seq_names.size();graph.seq_names.push_back(p.first);
        struct lps l;init_lps(&l,tumor_seq.c_str(),tumor_seq.size());lps_deepen(&l,5);
        graph.addString(&l,1,tumor_seq.c_str(),sid);free_lps(&l);
    }}
    BubbleResult markers=graph.detectBubbles();
    uint64_t n_only=13404391996925811288ULL;
    cout << "N-only core (hash 1340..1288):" << endl;
    cout << "  in bubble_ends? " << markers.bubble_ends.count(n_only) << endl;
    cout << "  in bubble_starts? " << markers.bubble_starts.count(n_only) << endl;
    auto it=graph.nodes.find(n_only);
    cout << "  color=" << it->second.colors << " (1=N,2=T,3=N+T)" << endl;
    // Check its predecessors and successors
    auto pit=graph.prev_edges.find(n_only);
    if(pit!=graph.prev_edges.end()){
        cout << "  predecessors:" << endl;
        for(auto ph:pit->second){
            auto pn=graph.nodes.find(ph);
            cout << "    hash=" << ph << " color=" << pn->second.colors << endl;
        }
    }
    return 0;
}

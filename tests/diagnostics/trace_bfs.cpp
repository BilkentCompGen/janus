#include <iostream>
#include <string>
#include <queue>
#include <unordered_set>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/utils.h"
using namespace std;
string colorStr(uint32_t c){bool n=(c&0x1),t=(c&0x2);if(n&&t)return"N+T";if(n)return"N";if(t)return"T";return"?";}
int main(){
    LCP_INIT();
    LCPGraph graph;
    RefIndex ref_index;
    string ref_seq;
    {FastaStreamer s("simulation/ref.fasta",false);
     for(const FastaPair&p:s){
        string c=p.first;size_t hp=c.find("_hap");if(hp!=string::npos)c=c.substr(0,hp);
        if(c!="19")continue;ref_seq=p.second;
        struct lps l;init_lps(&l,ref_seq.c_str(),ref_seq.size());lps_deepen(&l,5);
        ref_index.addString(&l,ref_seq.c_str());
        uint32_t sid=graph.seq_names.size();graph.seq_names.push_back(p.first);
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

    BubbleResult markers = graph.detectBubbles();
    uint64_t start_hash = 4996526221544266035ULL;
    uint64_t t1_hash = 7856366017321890231ULL;
    uint64_t t2_hash = 15783397869307142715ULL;
    uint64_t end_hash = 17856515116681298403ULL;

    cout << "start_hash in bubble_starts? " << markers.bubble_starts.count(start_hash) << endl;
    cout << "end_hash in bubble_ends? " << markers.bubble_ends.count(end_hash) << endl;
    cout << "t1_hash in bubble_ends? " << markers.bubble_ends.count(t1_hash) << endl;
    cout << "t2_hash in bubble_ends? " << markers.bubble_ends.count(t2_hash) << endl;
    cout << "t1_hash in bubble_starts? " << markers.bubble_starts.count(t1_hash) << endl;
    cout << "t2_hash in bubble_starts? " << markers.bubble_starts.count(t2_hash) << endl;

    // Manually simulate findColoredPath for tumor (color_bit=1)
    uint32_t target = 0x2, avoid = 0x1;
    cout << "\n=== Manual BFS trace (tumor path) ===" << endl;

    auto it_start = graph.next_edges.find(start_hash);
    const Node& sn = graph.nodes.at(start_hash);
    cout << "Start node successors:" << endl;
    for (node_hash n : it_start->second) {
        auto it = graph.nodes.find(n);
        bool is_end = markers.bubble_ends.count(n) > 0;
        bool has_target = (it->second.colors & target) != 0;
        bool has_avoid = (it->second.colors & avoid) != 0;
        cout << "  hash=" << n << " color=" << colorStr(it->second.colors)
             << " is_bubble_end=" << is_end
             << " has_target=" << has_target
             << " has_avoid=" << has_avoid;
        if (is_end) cout << " → RETURN IMMEDIATELY (empty seq)";
        else if (has_target && !has_avoid) cout << " → ENQUEUE";
        else cout << " → SKIP";
        cout << endl;
    }

    return 0;
}

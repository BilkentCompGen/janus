#include <iostream>
#include <string>
#include <iomanip>
#include "lib/lcptools/lps.h"
#include "lib/lcptools/encoding.h"
#include "src/utils.h"

using namespace std;

void print_bits(ublock hash) {
    for (int i = 63; i >= 0; i--) {
        cout << ((hash >> i) & 1);
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <string> <depth>" << endl;
        return 1;
    }
    
    string input_string = argv[1];
    int depth = stoi(argv[2]);
    
    LCP_INIT();
    
    struct lps lcp_string;
    init_lps(&lcp_string, input_string.c_str(), input_string.size());
    lps_deepen(&lcp_string, depth);
    
    cout << "LABEL\tBIT_REP\tLOCATIONS\tCHARS" << endl;
    cout << string(80, '-') << endl;
    
    for (int i = 0; i < lcp_string.size; i++) {
        ublock hash = hash_core(&lcp_string.cores[i]);
        uint64_t start = lcp_string.cores[i].start;
        uint64_t end = lcp_string.cores[i].end;
        
        // Extract characters
        string chars = input_string.substr(start, end - start);
        
        cout << i << "\t";
        print_bits(hash);
        cout << "\t" << start << "-" << end;
        cout << "\t" << chars << endl;
    }
    
    free_lps(&lcp_string);
    return 0;
}

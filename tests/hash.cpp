// src/test_lcp.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <cstring>

// LCP headers
#include "../lib/lcptools/encoding.h"
#include "../lib/lcptools/lps.h"

#include "../src/utils.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <depth_k>" << std::endl;
        return 1;
    }
    
    std::string fasta_file = argv[1];
    int depth = std::atoi(argv[2]);
    
    LCP_INIT();

    std::unordered_map<uint64_t, std::unordered_set<std::string>> hash_to_seqs;
    std::unordered_map<std::string, int> seq_count;
    size_t total_cores = 0;
    size_t total_length = 0;

    FastaStreamer streamer(fasta_file, false);
    for (auto& entry : streamer) {
        const std::string& name = entry.first;
        const std::string& sequence = entry.second;
        if (sequence.empty()) continue;
        
        std::cout << "Processing record: " << name << " (" << sequence.size() << " bp)" << std::endl;
        total_length += sequence.size();
        
        struct lps lcp_string;
        init_lps(&lcp_string, sequence.c_str(), sequence.size());
        lps_deepen(&lcp_string, depth);
        
        total_cores += lcp_string.size;
        
        for (int i = 0; i < lcp_string.size; ++i) {
            const core& c = lcp_string.cores[i];
            uint64_t length = c.end - c.start;
            
            std::string region(sequence.c_str() + c.start, length);
            
            uint64_t hash = hash_string(region, 0);
            
            seq_count[region]++;
            hash_to_seqs[hash].insert(region);
        }
        
        free_lps(&lcp_string);
    }
    
    std::cout << "Total sequence length: " << total_length << std::endl;
    std::cout << "Depth k: " << depth << std::endl;
    std::cout << "Total number of cores: " << total_cores << std::endl;
    
    int diff_seq_same_hash = 0;
    for (const auto& [hash, seqs] : hash_to_seqs) {
        if (seqs.size() > 1) {
            diff_seq_same_hash++;
        }
    }
    
    int repeated_sequence = 0;
    for (const auto& [seq, count] : seq_count) {
        if (count > 1) {
            repeated_sequence++;
        }
    }
    
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "1. Different sequences with same hash: " << diff_seq_same_hash << std::endl;
    std::cout << "2. Same sequence repeated >1 time: " << repeated_sequence << std::endl;
    std::cout << "\nUnique sequences: " << seq_count.size() << std::endl;
    std::cout << "Unique hashes: " << hash_to_seqs.size() << std::endl;
    
    return 0;
}
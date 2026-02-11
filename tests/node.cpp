#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <cassert>
#include <bitset>
#include <cstring>

inline uint64_t rotl64(uint64_t x, int8_t r) {
    return (x << r) | (x >> (64 - r));
}

inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

inline uint64_t load_u64_le(const uint8_t* p) {
    // Load little-endian 64-bit value
    uint64_t result;
    memcpy(&result, p, sizeof(result));
    return result;
}

// Returns h1 (low 64 bits) of MurmurHash3_x64_128
uint64_t MurmurHash3_64(const void* key, size_t len, uint64_t seed) {
    const uint8_t* data = (const uint8_t*)key;
    size_t nblocks = len / 16;
    
    uint64_t h1 = seed;
    uint64_t h2 = seed;
    
    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;
    
    for (size_t i = 0; i < nblocks; i++) {
        uint64_t k1 = load_u64_le(data + i * 16 + 0);
        uint64_t k2 = load_u64_le(data + i * 16 + 8);
        
        k1 *= c1;
        k1 = rotl64(k1, 31);
        k1 *= c2;
        h1 ^= k1;
        
        h1 = rotl64(h1, 27);
        h1 += h2;
        h1 = h1 * 5 + 0x52dce729;
        
        k2 *= c2;
        k2 = rotl64(k2, 33);
        k2 *= c1;
        h2 ^= k2;
        
        h2 = rotl64(h2, 31);
        h2 += h1;
        h2 = h2 * 5 + 0x38495ab5;
    }
    
    const uint8_t* tail = data + nblocks * 16;
    uint64_t k1 = 0;
    uint64_t k2 = 0;
    
    switch (len & 15) {
        case 15: k2 ^= ((uint64_t)tail[14]) << 48; [[fallthrough]];
        case 14: k2 ^= ((uint64_t)tail[13]) << 40; [[fallthrough]];
        case 13: k2 ^= ((uint64_t)tail[12]) << 32; [[fallthrough]];
        case 12: k2 ^= ((uint64_t)tail[11]) << 24; [[fallthrough]];
        case 11: k2 ^= ((uint64_t)tail[10]) << 16; [[fallthrough]];
        case 10: k2 ^= ((uint64_t)tail[9]) << 8; [[fallthrough]];
        case  9: k2 ^= ((uint64_t)tail[8]) << 0; k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2; [[fallthrough]];
        case  8: k1 ^= ((uint64_t)tail[7]) << 56; [[fallthrough]];
        case  7: k1 ^= ((uint64_t)tail[6]) << 48; [[fallthrough]];
        case  6: k1 ^= ((uint64_t)tail[5]) << 40; [[fallthrough]];
        case  5: k1 ^= ((uint64_t)tail[4]) << 32; [[fallthrough]];
        case  4: k1 ^= ((uint64_t)tail[3]) << 24; [[fallthrough]];
        case  3: k1 ^= ((uint64_t)tail[2]) << 16; [[fallthrough]];
        case  2: k1 ^= ((uint64_t)tail[1]) << 8; [[fallthrough]];
        case  1: k1 ^= ((uint64_t)tail[0]) << 0; k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1; [[fallthrough]];
        case  0: break;
    }
    
    h1 ^= (uint64_t)len;
    h2 ^= (uint64_t)len;
    
    h1 += h2;
    h2 += h1;
    
    h1 = fmix64(h1);
    h2 = fmix64(h2);
    
    h1 += h2;
    
    return h1;
}

// Wrapper for std::string
uint64_t murmurHash64(const std::string& data, uint64_t seed = 0) {
    return MurmurHash3_64(data.data(), data.size(), seed);
}

using node_hash = uint64_t;

struct Node {
    node_hash hash;
    uint32_t colors;
    std::string sequence;
    std::vector<uint64_t> ref_starts;
    
    Node() : hash(0), colors(0) {}
    
    Node(const std::string& seq, uint64_t seed = 0) : colors(0) {
        sequence = seq;
        hash = murmurHash64(sequence, seed);
    }
    
    Node(node_hash h, const std::string& seq) : hash(h), colors(0), sequence(seq) {}
    
    bool merge(const Node& next, uint64_t seed) {
        const std::string& pat = next.sequence;
        size_t n = pat.size();

        if (n == 0 || colors != next.colors) return false;
        
        std::vector<int> pi(n);
        for (size_t i = 1, j = 0; i < n; i++) {
            while (j > 0 && pat[i] != pat[j]) j = pi[j - 1];
            if (pat[i] == pat[j]) j++;
            pi[i] = j;
        }
        
        size_t m = sequence.size();
        size_t check_len = std::min(m, n);
        
        size_t j = 0;
        for (size_t i = m - check_len; i < m; i++) {
            while (j > 0 && sequence[i] != pat[j]) j = pi[j - 1];
            if (sequence[i] == pat[j]) j++;
        }
        
        sequence += pat.substr(j);
        hash = murmurHash64(sequence, seed);
        return true;
    }
    
    bool merge(const Node& next) {
        return merge(next, 0);
    }
    
    bool hasColor(int color) const {
        return (colors >> color) & 1;
    }
    
    void setColor(int color) {
        colors |= (1U << color);
    }
    
    void updateHash(uint64_t seed = 0) {
        hash = murmurHash64(sequence, seed);
    }
    
    void printInfo(const std::string& label = "") const {
        if (!label.empty()) {
            std::cout << label << ":\n";
        }
        std::cout << "  Hash: 0x" << std::hex << hash << std::dec << "\n";
        std::cout << "  Colors: " << std::bitset<32>(colors) << "\n";
        std::cout << "  Sequence: " << sequence << " (length: " << sequence.length() << ")\n";
        std::cout << "  Ref starts: ";
        for (size_t i = 0; i < ref_starts.size(); i++) {
            std::cout << ref_starts[i];
            if (i < ref_starts.size() - 1) std::cout << ", ";
        }
        std::cout << "\n" << std::endl;
    }
};

// Test function to verify MurmurHash consistency
void testMurmurHashConsistency() {
    std::cout << "Testing MurmurHash3_64 consistency\n";
    std::cout << "==================================\n\n";
    
    // Test 1: Same input should give same hash
    std::string seq1 = "AGCTAGCT";
    uint64_t hash1 = murmurHash64(seq1);
    uint64_t hash2 = murmurHash64(seq1);
    assert(hash1 == hash2);
    std::cout << "✓ Same sequence → same hash: 0x" << std::hex << hash1 << std::dec << "\n";
    
    // Test 2: Different seeds give different hashes
    uint64_t hash3 = murmurHash64(seq1, 1);
    assert(hash1 != hash3);
    std::cout << "✓ Different seeds → different hashes\n";
    
    // Test 3: Concatenation hash vs merge hash
    Node node1("AGCT", 42);
    Node node2("CTAG", 42);
    
    // Merge them
    node1.merge(node2, 42);
    
    // Create a node with the full sequence
    Node full("AGCTAG", 42);
    
    std::cout << "\nMerge test:\n";
    std::cout << "Node1 after merge hash: 0x" << std::hex << node1.hash << std::dec << "\n";
    std::cout << "Full sequence hash:     0x" << std::hex << full.hash << std::dec << "\n";
    
    // They should match!
    assert(node1.hash == full.hash);
    std::cout << "✓ Merge produces correct hash!\n\n";
}

// Test the KMP overlap finding with hash updates
void testMergeWithOverlap() {
    std::cout << "Testing merge with overlap detection\n";
    std::cout << "====================================\n\n";
    
    // Case 1: Partial overlap
    Node node1("AGCTAGCT", 0);
    Node node2("AGCTTAGG", 0);
    
    std::cout << "Before merge:\n";
    node1.printInfo("Node 1 (AGCTAGCT)");
    node2.printInfo("Node 2 (AGCTTAGG)");
    
    node1.merge(node2);
    
    std::cout << "After merge:\n";
    node1.printInfo("Merged Node");
    
    // Verify sequence
    assert(node1.sequence == "AGCTAGCTTAGG");
    
    // Verify hash matches direct computation
    Node expected("AGCTAGCTTAGG", 0);
    assert(node1.hash == expected.hash);
    std::cout << "✓ Partial overlap merge correct!\n\n";
    
    // Case 2: Full overlap
    Node node3("ABCDEF", 0);
    Node node4("ABCDEF", 0);
    
    node3.merge(node4);
    assert(node3.sequence == "ABCDEF");
    Node expected2("ABCDEF", 0);
    assert(node3.hash == expected2.hash);
    std::cout << "✓ Full overlap merge correct!\n\n";
    
    // Case 3: No overlap
    Node node5("AAAAA", 0);
    Node node6("BBBBB", 0);
    
    node5.merge(node6);
    assert(node5.sequence == "AAAAABBBBB");
    Node expected3("AAAAABBBBB", 0);
    assert(node5.hash == expected3.hash);
    std::cout << "✓ No overlap merge correct!\n\n";
}

// Test with colors and ref_starts (uncomment merge code to use)
void testAdvancedFeatures() {
    std::cout << "Testing with colors and ref_starts\n";
    std::cout << "===================================\n\n";
    
    Node node1("ACGT", 0);
    node1.setColor(0);
    node1.setColor(2);
    node1.ref_starts.push_back(100);
    node1.ref_starts.push_back(200);
    
    Node node2("GTAC", 0);
    node2.setColor(1);
    node2.setColor(3);
    node2.ref_starts.push_back(300);
    
    std::cout << "Before merge:\n";
    node1.printInfo("Node 1");
    node2.printInfo("Node 2");
    
    // Merge (without merging colors/ref_starts in current implementation)
    node1.merge(node2);
    
    std::cout << "After merge (colors and ref_starts not merged):\n";
    node1.printInfo("Merged Node");
    
    // Colors and ref_starts should remain unchanged
    assert(node1.hasColor(0));
    assert(node1.hasColor(2));
    assert(!node1.hasColor(1));  // Not merged
    assert(!node1.hasColor(3));  // Not merged
    assert(node1.ref_starts.size() == 2);  // Original ref_starts
    
    std::cout << "✓ Colors and ref_starts not merged (current behavior)\n";
    std::cout << "  To merge them, uncomment the merge code in Node::merge()\n\n";
}
// Test functions
void testColorMatchingRequirement() {
    std::cout << "Testing color matching requirement\n";
    std::cout << "===================================\n\n";
    
    // Test 1: Same colors - should merge
    Node node1("AGCT", 0);
    node1.setColor(1);
    node1.setColor(3);
    node1.ref_starts.push_back(100);
    
    Node node2("CTAG", 0);
    node2.setColor(1);
    node2.setColor(3);
    node2.ref_starts.push_back(200);
    
    std::cout << "Test 1 - Same colors (should merge):\n";
    node1.printInfo("Node 1 before");
    node2.printInfo("Node 2");
    
    bool result = node1.merge(node2);
    assert(result == true);
    assert(node1.sequence == "AGCTAG");  // Overlap "CT"
    assert(node1.colors == 0b1010);      // Colors unchanged
    assert(node1.ref_starts.size() == 1); // ref_starts unchanged
    assert(node1.ref_starts[0] == 100);   // Original ref_start
    
    std::cout << "After successful merge:\n";
    node1.printInfo("Node 1 after");
    std::cout << "✓ Same colors - merge successful\n\n";
    
    // Test 2: Different colors - should NOT merge
    Node node3("AAAA", 0);
    node3.setColor(0);
    node3.setColor(2);
    
    Node node4("BBBB", 0);
    node4.setColor(1);
    node4.setColor(3);
    
    std::cout << "Test 2 - Different colors (should NOT merge):\n";
    node3.printInfo("Node 3 before");
    node4.printInfo("Node 4");
    
    std::string original_seq = node3.sequence;
    uint32_t original_colors = node3.colors;
    
    result = node3.merge(node4);
    assert(result == false);  // Merge should fail
    assert(node3.sequence == original_seq);  // Sequence unchanged
    assert(node3.colors == original_colors); // Colors unchanged
    
    std::cout << "After failed merge:\n";
    node3.printInfo("Node 3 after");
    std::cout << "✓ Different colors - merge rejected\n\n";
    
    // Test 3: Subset colors - should NOT merge
    Node node5("CCCC", 0);
    node5.setColor(0);
    node5.setColor(1);
    node5.setColor(2);
    
    Node node6("DDDD", 0);
    node6.setColor(0);
    node6.setColor(1);
    // Missing color 2
    
    std::cout << "Test 3 - Subset colors (should NOT merge):\n";
    node5.printInfo("Node 5");
    node6.printInfo("Node 6");
    
    result = node5.merge(node6);
    assert(result == false);  // Colors must be exactly equal
    
    std::cout << "✓ Subset colors - merge rejected\n\n";
    
    // Test 4: Empty colors - should merge if both empty
    Node node7("EEEE", 0);
    Node node8("FFFF", 0);
    
    std::cout << "Test 4 - Both empty colors (should merge):\n";
    node7.printInfo("Node 7");
    node8.printInfo("Node 8");
    
    result = node7.merge(node8);
    assert(result == true);
    
    std::cout << "After merge:\n";
    node7.printInfo("Node 7 after");
    std::cout << "✓ Empty colors - merge successful\n\n";
}

void testRefStartsUnchanged() {
    std::cout << "Testing ref_starts remain unchanged\n";
    std::cout << "====================================\n\n";
    
    Node node1("GATTACA", 0);
    node1.setColor(1);
    node1.ref_starts.push_back(1000);
    node1.ref_starts.push_back(2000);
    node1.ref_starts.push_back(3000);
    
    Node node2("ACAGATT", 0);
    node2.setColor(1);
    node2.ref_starts.push_back(4000);
    node2.ref_starts.push_back(5000);
    
    std::cout << "Before merge:\n";
    node1.printInfo("Node 1");
    node2.printInfo("Node 2");
    
    bool result = node1.merge(node2);
    assert(result == true);
    
    // Verify ref_starts unchanged
    assert(node1.ref_starts.size() == 3);
    assert(node1.ref_starts[0] == 1000);
    assert(node1.ref_starts[1] == 2000);
    assert(node1.ref_starts[2] == 3000);
    
    std::cout << "After merge:\n";
    node1.printInfo("Node 1");
    std::cout << "✓ ref_starts remain unchanged\n\n";
}

int main() {
    std::cout << "=== Node Merge with MurmurHash3 Testing ===\n\n";
    
    testMurmurHashConsistency();
    testMergeWithOverlap();
    testAdvancedFeatures();
    testColorMatchingRequirement();
    testRefStartsUnchanged();
    
    std::cout << "All tests passed! ✓\n";
    return 0;
}
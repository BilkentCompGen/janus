#include <algorithm>
#include <chrono>
#include <cctype>
#include <ctime>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "zstr.hpp"

#include "core.h"
#include "utils.h"

using namespace std;
using namespace chrono;

//////////////////////////////////////////////////////////// Logger ////////////////////////////////////////////////////////////
Logger::Logger(const string& name) : name(name) {}

void Logger::log_impl(const string& msg, const char* file, int line) const {
    auto now = chrono::system_clock::now();
    auto time_t = chrono::system_clock::to_time_t(now);
    auto ms = chrono::duration_cast<chrono::milliseconds>(now.time_since_epoch()) % 1000;
    
    stringstream ss;
    ss << put_time(localtime(&time_t), "%H:%M:%S");
    ss << "." << setfill('0') << setw(3) << ms.count();
    
    cout << ss.str() << " | " << name << " | " << file << ":" << line << " | " << msg << endl;
}

/////////////////////////////////////////////////////////// Gzip I/O ///////////////////////////////////////////////////////////
unique_ptr<istream> open_input_stream(const string& file_path, bool is_gzipped) {
    if (is_gzipped) return make_unique<gzifstream>(file_path);
    else return make_unique<ifstream>(file_path);
}

////////////////////////////////////////////////////////// FASTA Reader ////////////////////////////////////////////////////////
static inline void ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !isspace(ch);
    }));
}

static inline void rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !isspace(ch);
    }).base(), s.end());
}

static inline void trim(string &s) {
    ltrim(s);
    rtrim(s);
}

FastaStreamer::FastaStreamer(const string& file_path, bool is_gzipped) {
    input_stream = open_input_stream(file_path, is_gzipped);
    if (!*input_stream) {
        throw runtime_error("Could not open file: " + file_path);
    }
}

FastaIterator FastaStreamer::begin() {
    return FastaIterator(input_stream.get());
}

FastaIterator FastaStreamer::end() {
    return FastaIterator(nullptr);
}

FastaIterator::FastaIterator(istream* in) : input(in), at_end(in == nullptr) {
    if (!at_end) read_next_sequence();
}

FastaIterator& FastaIterator::operator++() {
    read_next_sequence();
    return *this;
}

const FastaPair& FastaIterator::operator*() const {
    return current_pair;
}

bool FastaIterator::operator!=(const FastaIterator& other) const {
    return at_end != other.at_end;
}

void FastaIterator::read_next_sequence() {
    string line;
    current_pair.first.clear();
    current_pair.second.clear();

    if (input->peek() == EOF) {
        at_end = true;
        return;
    }

    while (getline(*input, line)) {
        trim(line);
        if (line.empty()) continue;
        if (line[0] == '>') {
            current_pair.first = line.substr(1);
            break;
        }
    }

    if (current_pair.first.empty()) {
        at_end = true;
        return;
    }

    while (input->peek() != '>' && input->peek() != EOF) {
        getline(*input, line);
        trim(line);
        if(!line.empty()) current_pair.second += line;
    }
}

//////////////////////////////////////////////////////// Argument Parser /////////////////////////////////////////////////////
ArgParser::ArgParser(const string& description) : program_description(description) {}

void ArgParser::add(const string& name, const string& flag, const string& short_flag, const string& description, const string& type, bool required) {
    definitions[name] = {flag, short_flag, description, type, required};
    flag_to_name[flag] = name;
    flag_to_name[short_flag] = name;
}

void ArgParser::print_help() {
    cerr << "Usage: " << program_name << " [options]" << endl;
    if (!program_description.empty()) cerr << endl << program_description << endl;
    cerr << endl << "Options:" << endl;
    for (const pair<const string, Argument>& pair : definitions) {
        const Argument& arg = pair.second;
        cerr << "  " << arg.flag << ", " << arg.short_flag << "\t" << arg.description << " (type: " << arg.type << ")" << (arg.required ? " [required]" : "") << endl;
    }
    cerr << "  --help, -h\tShow this help message and exit" << endl;
    exit(0);
}

map<string, any> ArgParser::parse(int argc, char** argv) {
    program_name = argv[0];
    if (argc == 1) print_help();
    
    map<string, any> results;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") print_help();
        
        if (!flag_to_name.count(arg)) continue;
        
        string name = flag_to_name[arg];
        const Argument& def = definitions[name];
        
        if (def.type == "bool") {
            results[name] = true;
        } else {
            if (i + 1 >= argc) throw runtime_error("Missing value for argument: " + arg);
            
            string value_str = argv[++i];
            
            if (def.type == "string") {
                results[name] = value_str;
            } else if (def.type == "int") {
                try {
                    results[name] = stoi(value_str);
                } catch (const invalid_argument&) {
                    throw runtime_error("Invalid integer value for argument: " + arg);
                }
            } else if (def.type == "float") {
                try {
                    results[name] = stof(value_str);
                } catch (const invalid_argument&) {
                    throw runtime_error("Invalid float value for argument: " + arg);
                }
            }
        }
    }
    
    for (const auto& [name, def] : definitions) {
        if (def.required && !results.count(name)) {
            throw runtime_error("Missing required argument: " + def.flag);
        }
    }
    
    return results;
}

//////////////////////////////////////////////////////// LCP Core Hashing ///////////////////////////////////////////////////////////

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

uint64_t hash_string(const std::string& data, uint64_t seed) {
    return MurmurHash3_64(data.data(), data.size(), seed);
}

//////////////////////////////////////////////////////// CSV Export ///////////////////////////////////////////////////////////

void writeCSV(const std::map<std::string, std::vector<std::string>>& data, const std::string& filename, const std::vector<std::string>& headers, char separator) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Could not open file for writing: " + filename);
    
    size_t max_cols = 0;
    for (const auto& [key, vec] : data) max_cols = std::max(max_cols, vec.size());
    
    if (!headers.empty()) {
        for (size_t i = 0; i < headers.size(); i++) {
            if (i > 0) out << separator;
            out << headers[i];
        }
        out << "\n";
    } else throw std::runtime_error("Headers are required for CSV export");
    
    for (const auto& [key, vec] : data) {
        out << key;
        for (size_t i = 0; i < max_cols; i++) {
            out << separator;
            if (i < vec.size()) out << vec[i];
            else out << "";
        }
        out << "\n";
    }
}
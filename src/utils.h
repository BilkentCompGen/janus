#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <any>
#include <memory>

class Logger {
public:
    Logger(const std::string& name);
    void log_impl(const std::string& msg, const char* file, int line) const;
private:
    std::string name;
};

#define debug(msg) log_impl(msg, __FILE__, __LINE__)

std::unique_ptr<std::istream> open_input_stream(const std::string& file_path, bool is_gzipped = false);

using FastaPair = std::pair<std::string, std::string>;

class FastaIterator;

class FastaStreamer {
public:
    FastaStreamer(const std::string& file_path, bool is_gzipped = false);
    FastaIterator begin();
    FastaIterator end();
private:
    std::unique_ptr<std::istream> input_stream;
};

class FastaIterator {
public:
    FastaIterator(std::istream* in);
    FastaIterator& operator++();
    const FastaPair& operator*() const;
    bool operator!=(const FastaIterator& other) const;

    using iterator_category = std::input_iterator_tag;
    using value_type = FastaPair;
    using difference_type = std::ptrdiff_t;
    using pointer = FastaPair*;
    using reference = FastaPair&;

private:
    std::istream* input;
    FastaPair current_pair;
    bool at_end;
    void read_next_sequence();
};

std::vector<FastaPair> read_fasta_file(const std::string& file_path, bool is_gzipped = false);

struct Argument {
    std::string flag;
    std::string short_flag;
    std::string description;
    std::string type;
    bool required;
};

class ArgParser {
public:
    ArgParser(const std::string& description = "");
    void add(const std::string& name, const std::string& flag, const std::string& short_flag, const std::string& description, const std::string& type, bool required = false);
    std::map<std::string, std::any> parse(int argc, char** argv);
private:
    std::string program_name;
    std::string program_description;
    std::map<std::string, Argument> definitions;
    std::map<std::string, std::string> flag_to_name;
    void print_help();
};

uint32_t hash_core(const struct core* c);

uint32_t hash_k_cores(const struct core* cores, int k);

uint32_t update_k_cores_hash(uint32_t old_hash, const struct core* old_first, const struct core* new_last, int k);

#endif // UTILS_H

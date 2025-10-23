#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <any>
#include <memory>

// Logger
// Usage:
// Logger logger = get_logger("my_logger");
class Logger {
public:
    Logger(const std::string& name);
    void log_impl(const std::string& msg, const char* file, int line) const;
private:
    std::string name;
};

// Use as logger.log("message")
#define debug(msg) log_impl(msg, __FILE__, __LINE__)

// Gzip I/O
// Returns a unique_ptr to an istream that can read from a plain or gzipped file.
std::unique_ptr<std::istream> open_input_stream(const std::string& file_path, bool is_gzipped = false);

// FASTA Reader
using FastaPair = std::pair<std::string, std::string>;

// FastaStreamer provides an iterator-based way to read a FASTA file.
// It does not load the entire file into memory.
// Usage:
// FastaStreamer streamer("path/to/file.fasta");
// for (const FastaPair& pair : streamer) {
//     std::cout << "Header: " << pair.first << ", Sequence: " << pair.second << std::endl;
// }
class FastaIterator; // Forward declaration

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

    // Standard iterator traits
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

// Reads an entire FASTA file into a vector.
std::vector<FastaPair> read_fasta_file(const std::string& file_path, bool is_gzipped = false);


// Argument Parser
// A simple argument parser.
// Usage:
// ArgParser parser;
// parser.add("input", "--input", "-i", "Input file", "string", true);
// parser.add("threads", "--threads", "-t", "Number of threads", "int", false);
// std::map<std::string, std::any> args = parser.parse(argc, argv);
// std::string input_file = std::any_cast<std::string>(args["input"]);
// int threads = args.count("threads") ? std::any_cast<int>(args["threads"]) : 1;
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

#endif // UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <any>

#include "src/utils.h"

using namespace std;

Logger logger("parser_and_reader_logger");

int main(int argc, char** argv) {
    ArgParser parser("A simple program to demonstrate the ArgParser and FastaStreamer classes.");

    parser.add("input", "--input", "-i", "Input FASTA file", "string", true);
    parser.add("gzip", "--gzip", "-g", "Gzipped input FASTA file", "bool", false);

    map<string, any> args = parser.parse(argc, argv);
    
    string input_file = any_cast<string>(args["input"]);
    bool gzip = args.count("gzip") ? any_cast<bool>(args["gzip"]) : false;
    
    logger.debug("Input file: " + input_file);
    logger.debug("Gzip: " + to_string(gzip));
    
    logger.debug("Testing FASTA streamer...");
    FastaStreamer streamer(input_file, gzip);
    for (const FastaPair& pair : streamer) {
        cout << pair.first << ": " << pair.second.size() << endl;
    }
    
    logger.debug("Testing FASTA file reader...");
    vector<FastaPair> fasta_data = read_fasta_file(input_file, gzip);
    logger.debug("Read " + to_string(fasta_data.size()) + " sequences.");

    return 0;
}

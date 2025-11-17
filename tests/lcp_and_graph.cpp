#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <any>
#include <unordered_set>

#include "src/utils.h"
#include "src/graph.h"
#include "lps.h"

using namespace std;

Logger main_logger("main_logger");
Logger parser_logger("parser_logger");
Logger lcp_logger("lcp_logger");
Logger graph_logger("graph_logger");

int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("A simple program to demonstrate the ArgParser and FastaStreamer classes.");

    parser.add("input", "--input", "-i", "Input FASTA file", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("output", "--output", "-o", "Output file name", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string input_file = any_cast<string>(args["input"]);
    int depth = any_cast<int>(args["depth"]);
    string output_file = any_cast<string>(args["output"]);

    main_logger.debug("Input file: " + input_file);
    main_logger.debug("Depth: " + to_string(depth));
    main_logger.debug("Output file: " + output_file);
    
    graph_logger.debug("Building graph...");
    LCPGraph graph;
    FastaStreamer streamer(input_file, false);
    for (const FastaPair& pair : streamer) {
        main_logger.debug(pair.first + ": " + to_string(pair.second.size()));
        lcp_logger.debug("Initializing LCP...");
        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        lcp_logger.debug("LCP size: " + to_string(lcp_string.size));
        graph_logger.debug("Adding string to graph...");
        graph.addString(&lcp_string, 0);
        graph_logger.debug("String added to graph.");
        free_lps(&lcp_string);
        break;
    }
    
    graph_logger.debug("Graph built.");
    graph.exportGFA(output_file);
    graph_logger.debug("Graph exported to " + output_file + ".");

    main_logger.debug("Done!");

    return 0;
}

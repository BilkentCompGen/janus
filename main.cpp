#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <any>

#include "src/utils.h"
#include "lps.h"

using namespace std;

Logger logger("main_logger");

int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("A simple program to demonstrate the ArgParser and FastaStreamer classes.");

    parser.add("str", "--str", "-s", "Input string", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", false);

    map<string, any> args = parser.parse(argc, argv);

    string str = any_cast<string>(args["str"]);
    int depth = args.count("depth") ? any_cast<int>(args["depth"]) : 3;
    logger.debug("Input string: " + str);
    logger.debug("Input string size: " + to_string(str.size()));
    
    logger.debug("Testing LCP...");

    struct lps lcp_string;

    init_lps(&lcp_string, str.c_str(), str.size());

    int isSuccess = lps_deepen(&lcp_string, depth);

    logger.debug("LCP size: " + to_string(lcp_string.size));
    logger.debug("LCP isSuccess: " + to_string(isSuccess));

    logger.debug("Printing LCP...");
    for (int i = 0; i < lcp_string.size; i++) {
        print_core(&(lcp_string.cores[i]));
        printf(" ");
    }
    printf("\n");

    for (int i = 0; i < lcp_string.size; i++) {
        printf("(%d %d)", lcp_string.cores[i].start, lcp_string.cores[i].end);
        printf(" ");
    }
    printf("\n");

    logger.debug("Freeing LCP...");
    free_lps(&lcp_string);

    logger.debug("Done!");

    return 0;
}

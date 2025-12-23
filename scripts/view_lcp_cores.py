#!/usr/bin/env python3
"""
Visualize LCP cores for a given string.

Usage:
    python3 view_lcp_cores.py <string> [depth]
    
Example:
    python3 view_lcp_cores.py "ACGTACGT" 3
"""

import sys
import subprocess
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 view_lcp_cores.py <string> [depth]", file=sys.stderr)
        print("Example: python3 view_lcp_cores.py 'ACGTACGT' 3", file=sys.stderr)
        sys.exit(1)
    
    input_string = sys.argv[1]
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    cpp_binary = os.path.join(project_root, "view_lcp_cores")
    
    # Compile if needed
    cpp_source = os.path.join(script_dir, "view_lcp_cores.cpp")
    utils_cpp = os.path.join(project_root, "src", "utils.cpp")
    if not os.path.exists(cpp_binary) or os.path.getmtime(cpp_source) > os.path.getmtime(cpp_binary):
        print("Compiling C++ helper...", file=sys.stderr)
        compile_cmd = [
            "g++", "-std=c++20", "-O3", "-I", project_root,
            "-o", cpp_binary,
            cpp_source, utils_cpp,
            "-L", os.path.join(project_root, "lib", "lcptools"),
            "-llcptools", "-lz"
        ]
        result = subprocess.run(compile_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Compilation failed: {result.stderr}", file=sys.stderr)
            sys.exit(1)
    
    # Run the C++ program
    result = subprocess.run([cpp_binary, input_string, str(depth)], 
                          capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    
    # Parse output
    lines = result.stdout.strip().split('\n')
    if len(lines) < 2:
        print("No cores found", file=sys.stderr)
        return
    
    header = lines[0]
    data_lines = lines[2:]  # Skip header and separator
    
    cores = []
    for line in data_lines:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 4:
            label = parts[0]
            bit_rep = parts[1]
            locations = parts[2]
            chars = parts[3]
            cores.append((label, bit_rep, locations, chars))
    
    # Print formatted output
    print("=" * 80)
    print("LCP CORE VISUALIZATION")
    print("=" * 80)
    print(f"Input string: {input_string}")
    print(f"Depth: {depth}")
    print(f"Number of cores: {len(cores)}")
    print()
    
    # Print cores in horizontal format (columns)
    if cores:
        num_cores = len(cores)
        
        # Calculate column widths
        max_label_len = max(len(c[0]) for c in cores)
        max_bit_len = max(len(c[1]) for c in cores)
        max_loc_len = max(len(c[2]) for c in cores)
        max_char_len = max(len(c[3]) for c in cores)
        
        # Print LABEL row
        print("LABEL" + " " * (max_label_len - 4), end="")
        for label, _, _, _ in cores:
            print(f" {label:>{max_label_len}}", end="")
        print()
        
        # Print BIT REP row
        print("BIT" + " " * (max_label_len - 2), end="")
        for _, bit_rep, _, _ in cores:
            print(f" {bit_rep:<{max_bit_len}}", end="")
        print()
        
        # Print LOCATIONS row
        print("LOC" + " " * (max_label_len - 2), end="")
        for _, _, locations, _ in cores:
            print(f" {locations:<{max_loc_len}}", end="")
        print()
        
        # Print CHARS row
        print("CHAR" + " " * (max_label_len - 3), end="")
        for _, _, _, chars in cores:
            print(f" {chars:<{max_char_len}}", end="")
        print()
        print()
    
    # Print sequence visualization at the end
    print("=" * 80)
    print("SEQUENCE VISUALIZATION (4-line thickness):")
    print("=" * 80)
    
    # Create 4-line visualization
    seq = input_string
    line1 = "┌" + "─" * len(seq) + "┐"
    line2 = "│" + seq + "│"
    line3 = "│" + " " * len(seq) + "│"
    line4 = "└" + "─" * len(seq) + "┘"
    
    print(line1)
    print(line2)
    print(line3)
    print(line4)
    
    # Add position markers
    pos_line = " "
    for i in range(len(seq)):
        if i % 10 == 0:
            pos_line += str(i)
            pos_line += " " * (10 - len(str(i)))
        else:
            pos_line += " "
    print(pos_line[:len(seq)+2])

if __name__ == "__main__":
    main()

#!/bin/bash

# Generate PNG images for all GFAs under output/gfa using Bandage.
# Skips images that already exist.

BASE_DIR="outputs"
GFA_DIR="${BASE_DIR}/gfa"
IMG_DIR="${BASE_DIR}/images"

mkdir -p "$IMG_DIR"

if [ ! -d "$GFA_DIR" ]; then
    echo "GFA directory not found: $GFA_DIR"
    exit 1
fi

for gfa in "$GFA_DIR"/*.gfa; do
    [ -f "$gfa" ] || continue
    basename="$(basename "${gfa%.gfa}")"
    output_png="${IMG_DIR}/${basename}.png"

    if [ -f "$output_png" ]; then
        echo "Skipping $gfa (image already exists: $output_png)"
        continue
    fi

    # Find matching node colours CSV for this GFA:
    #   1) output/csv/<basename>.csv (preferred, written by main.cpp)
    node_colours_arg=()
    csv_match="${BASE_DIR}/csv/${basename}.csv"
    if [ -f "$csv_match" ]; then
        node_colours_arg=(--color "$csv_match")
    fi

    echo "Generating image for $gfa -> $output_png"
    Bandage image "$gfa" "$output_png" \
        --nodewidth 100 \
        --depwidth 0 \
        --width 2048 \
        --height 2048

    if [ $? -eq 0 ]; then
        echo "  ✓ Success: $output_png"
    else
        echo "  ✗ Failed: $gfa"
    fi
done

echo ""
echo "All images generated in $IMG_DIR/"
ls -lh "$IMG_DIR"/*.png 2>/dev/null

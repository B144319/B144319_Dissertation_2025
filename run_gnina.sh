#!/bin/bash

CONFIG_DIR=""
OUTPUT_DIR=""
mkdir -p "$OUTPUT_DIR"

for CONFIG in "$CONFIG_DIR"/*.txt; do
    BASENAME=$(basename "$CONFIG" .txt)
    echo "Running $CONFIG"
    gnina --config "$CONFIG" --out "$OUTPUT_DIR/${BASENAME}_out.pdb"
done
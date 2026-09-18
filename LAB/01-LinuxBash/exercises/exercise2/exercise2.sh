#!/bin/bash
# Exercise 2: count the fields of the ATOM lines, then extract atom numbers and names.
# Usage: ./exercise2.sh out-2RVD.pdb    (writes out-2RVD.txt)
#
# awk reads a file line by line and splits each line into fields $1, $2, ...
# NF is the number of fields on the current line. On an ATOM line, $2 is the
# atom number and $3 the atom name.

input=$1
output="$(basename "$input" .pdb).txt"

echo "Number of fields on the ATOM lines (count, value):"
awk '/^ATOM/ {print NF}' "$input" | sort | uniq -c

awk '/^ATOM/ {print $2, $3}' "$input" > "$output"

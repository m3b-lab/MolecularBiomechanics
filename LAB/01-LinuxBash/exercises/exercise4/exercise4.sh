#!/bin/bash
# Exercise 4: keep only the backbone atoms (N, CA, C) of an XYZ trajectory.
# Usage: ./exercise4.sh out-2RVD.xyz    (writes backbone-2RVD.xyz)
#
# awk walks through the file frame by frame: the first line of a frame gives
# the number of atoms n, the second is the title, then come n atom lines.
# It stores the backbone lines of each frame and prints them with the new count.

input=$1
output="backbone-$(basename "$input" | sed 's/^out-//')"

awk '
  state == 0 { n = $1; state = 1; next }
  state == 1 { title = $0; state = 2; k = 0; nb = 0; next }
  state == 2 {
      k++
      if ($1 == "N" || $1 == "CA" || $1 == "C") bb[++nb] = $0
      if (k == n) {
          print nb; print title
          for (i = 1; i <= nb; i++) print bb[i]
          state = 0
      }
  }
' "$input" > "$output"

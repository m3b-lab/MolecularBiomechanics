#!/bin/bash
# Exercise 6: RMSD of the CA atoms of every frame from the average structure.
# Usage: ./exercise6.sh backbone-2RVD.xyz average-2RVD.xyz    (writes rmsd.txt)
#
# RMSD = sqrt( sum over CA atoms of the squared distance / number of CA atoms )
# awk reads the average structure first (FNR == NR is true only for the first
# file) and stores its CA coordinates, then reads the trajectory frame by frame.

trajectory=$1
average=$2

awk '
  FNR == NR { if (FNR > 2 && $1 == "CA") { na++; ax[na] = $2; ay[na] = $3; az[na] = $4 }; next }
  state == 0 { n = $1; state = 1; next }
  state == 1 { frame++; state = 2; k = 0; i = 0; sum = 0; next }
  state == 2 {
      k++
      if ($1 == "CA") { i++; sum += ($2 - ax[i])^2 + ($3 - ay[i])^2 + ($4 - az[i])^2 }
      if (k == n) { printf "frame %2d  RMSD %.3f A\n", frame, sqrt(sum / i); total += sqrt(sum / i); state = 0 }
  }
  END { printf "mean      RMSD %.3f A\n", total / frame }
' "$average" "$trajectory" > rmsd.txt

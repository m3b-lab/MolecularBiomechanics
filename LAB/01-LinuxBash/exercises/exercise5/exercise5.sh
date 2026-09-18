#!/bin/bash
# Exercise 5: average position of every backbone atom over all frames.
# Usage: ./exercise5.sh backbone-2RVD.xyz    (writes average-2RVD.xyz)
#
# The 20 models of 2RVD are already superimposed on each other, so averaging
# the raw coordinates makes sense. For an MD trajectory you would first have
# to fit every frame onto a reference (LAB 4).

input=$1
output="average-$(basename "$input" | sed 's/^backbone-//')"

awk '
  state == 0 { n = $1; state = 1; next }
  state == 1 { frames++; state = 2; k = 0; next }
  state == 2 {
      k++
      name[k] = $1; x[k] += $2; y[k] += $3; z[k] += $4
      if (k == n) state = 0
  }
  END {
      print n
      print "average of", frames, "frames"
      for (i = 1; i <= n; i++) printf "%s %.3f %.3f %.3f\n", name[i], x[i] / frames, y[i] / frames, z[i] / frames
  }
' "$input" > "$output"

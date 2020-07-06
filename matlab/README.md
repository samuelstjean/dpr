# Diffusion profile realignment (dpr)

This is the matlab / [octave](https://www.gnu.org/software/octave/) version of the code as now available in [ExploreDTI](http://www.exploredti.com/).

See the example on how to run it as a standalone version and adapt it to your need.
In this case, we start from individual .mat files per subject which contain already extracted metrics.
The original algorithm works with a text file (one subject per line, each column is a different point), which are created from these mat files using ExploreDTI.

The example goes through loading the data, realiging, resampling, truncating and producing figures after statistical testing
and every step can be adapted by hand as needed.

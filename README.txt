This directory contains programs that implement the Reed-Muench formula for determining the TCID50 (tissue-culture infectious dose 50%) from viral titering in columns of a 96-well plate.

The Reed-Muench formula is implemented as described in the following references:
http://whqlibdoc.who.int/monograph/WHO_MONO_23_(3ed)_appendices.pdf
http://www.fao.org/docrep/005/ac802e/ac802e0w.htm

To use this program, simply run the self contained Python executable script reedmuenchcalculator.py. This script is tested with Python 2.6.7, but should work for most other versions of python 2.*. The input is read from a text file, and an example is shown in the file example_input.txt. An output text file is generated.

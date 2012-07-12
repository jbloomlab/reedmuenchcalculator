"""Calculates titer using Reed-Muench formula.

This script takes a text input file. It calculates the titer for
all of the samples, and prints them out to a new text file with the
name based on the input file. 

Written by Jesse Bloom, 2012."""


import os
import sys
import re
import math
import traceback



def Titer(infectedwells, volume, dilution):
    """Calculates the titer using the Reed-Muench formula.

    infectedwells -> This is a list of lists. The number of lists
        gives the number of replicates, so there are len(infectedwells)
        replicates. Each of the entry lists describes the wells with
        observed infection in rows of a 96-well plate. So for example,
            [[A, B, C, D], [A, B, C], [A, B, C, D]]
        corresponds to three replicates having cytopathic effect
        in the first four rows (first replicate), the first three rows
        (second replicate) and the first four rows (third replicate).
        There need to be at least two replicates.
    volume -> This is the infection volume in the first row (row A).
    dilution -> This is the dilution factor between successive rows.
        For example, 10 is a typical dilution factor for this assay.

    This method returns a number which gives the titer as TCID50
        per unit volume in whatever units are used to specify
        the input variable volume.

    The Reed-Muench formula is implemented as described in
        http://whqlibdoc.who.int/monograph/WHO_MONO_23_(3ed)_appendices.pdf
        http://www.fao.org/docrep/005/ac802e/ac802e0w.htm
    """
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] # row labels in order
    reverserows = [row for row in rows]
    reverserows.reverse()
    nreplicates = len(infectedwells)
    if nreplicates < 2:
        raise ValueError("This implementation of the Reed-Muench formula requires at least two replicates. Only %d are provided." % nreplicates)
    counts = dict([(r, 0) for r in rows]) # counts of infected wells at each dilution
    for replicatewells in infectedwells:
        for well in replicatewells:
            if well not in rows:
                raise ValueError("One of the rows is specified as %s, which is not a valid row." % well)
            counts[well] += 1
    infected = {} # cumulative totals of infected wells going up plate
    uninfected = {} # cumulative totals of uninfected wells going down plate
    n = 0
    for row in rows:
        uninfected[row] = n + nreplicates - counts[row]
        n = uninfected[row]
    n = 0
    for row in reverserows:
        infected[row] = n + counts[row]
        n = infected[row]
    percentinfected = {} # cumulative percent infected
    for row in rows:
        percentinfected[row] = 100.0 * infected[row] / (infected[row] + uninfected[row])
    for irow in range(len(rows)):
        if percentinfected[rows[irow]] < 50:
            if irow == 0:
                raise ValueError("Even the first dilution has < 50% infected.")
            else:
                rowabove50 = rows[irow - 1]
                break
    else:
        raise ValueError("No dilutions have < 50% infected.")
    percentrowabove50 = percentinfected[rowabove50]
    if rowabove50 != rows[-1]:
        percentrowbelow50 = percentinfected[rows[rows.index(rowabove50) + 1]]
    else:
        percentrowbelow50 = 0
    index = (percentrowabove50 - 50.0) / (percentrowabove50 - percentrowbelow50)
    startdilution = rows.index(rowabove50)
    titer = dilution**(startdilution + index) / volume
    return titer


def ParseInput(infile):
    """Reads an input text file.

    This file should be in the format specified by the example input file.
    The returned variable is the following tuple: (sampledata, volume, dilution)
    samplenames -> This is a list of the sample names in the order in
        which they appear in the input file.
    sampledata -> This is a dictionary. It is keyed by sample names.
        Each entry is a list of lists. The number of lists
        gives the number of replicates, so there are len(infectedwells)
        replicates. Each of the entry lists describes the wells with
        observed infection in rows of a 96-well plate. So for example,
        [[A, B, C, D], [A, B, C], [A, B, C, D]]
        corresponds to three replicates having cytopathic effect
        in the first four rows (first replicate), the first three rows
        (second replicate) and the first four rows (third replicate).
        There need to be at least two replicates.
    volume -> This is the infection volume in the first row (row A).
    dilution -> This is the dilution factor between successive rows.
        For example, 10 is a typical dilution factor for this assay.
    This method raises and IOError if the input file is not in the
        expected format or specifies invalid values.
    """
    lines = [line for line in open(infile).readlines() if line[0] != '#' and not line.isspace()]
    line1match = re.compile('^\s*VOLUME\s+(?P<volume>\d+\.{0,1}\d*)\s*\n$')
    m = line1match.search(lines[0])
    if not m:
        raise IOError("Failed to parse VOLUME from the first line.")
    volume = float(m.group('volume'))
    line2match = re.compile('^\s*DILUTION\s+(?P<dilution>\d+\.{0,1}\d*)\s*\n$')
    m = line2match.search(lines[1])
    if not m:
        raise IOError("Failed to parse DILUTION from the second line.")
    dilution = float(m.group('dilution'))
    if dilution <= 1:
        raise IOError("The dilution factor must be > 1, but read a value of %f" % dilution)
    line3match = re.compile('^\s*NREPLICATES\s+(?P<nreplicates>\d+)\s*\n$')
    m = line3match.search(lines[2])
    if not m:
        raise IOError("Failed to parse an integer value for NREPLICATES from the third line.")
    nreplicates = int(m.group('nreplicates'))
    if nreplicates < 2:
        raise IOError("There must be at least two replicates, but read a value of %d." % nreplicates)
    lines = lines[3 : ] # the remaining lines
    # there should be nreplicates + 1 line for each sample
    linespersample = nreplicates + 1
    if len(lines) % linespersample != 0:
        raise IOError("The sample data is not specified correctly. There should be a total of %d lines for each sample (the sample name plus a line for each of the %d replicates), but the number additional lines is not divisible by %d." % (linespersample, nreplicates, linespersample))
    nsamples = len(lines) / linespersample
    sampledata = {}
    namematch = re.compile('^\s*SAMPLE\s+(?P<name>.+)\n$')
    validrows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] 
    samplenames = []
    for isample in range(nsamples):
        nameline = lines[isample * linespersample]
        samplelines = lines[isample * linespersample + 1 : (isample + 1) * linespersample]
        assert len(samplelines) == nreplicates
        m = namematch.search(nameline)
        if not m:
            raise IOError("Failed to match sample name from line: %s" % nameline)
        sample = m.group('name').strip()
        if sample in sampledata:
            raise IOError("Duplicate sample name of %s" % sample)
        sampledata[sample] = []
        samplenames.append(sample)
        for line in samplelines:
            if line.strip() == 'na':
                sampledata[sample].append([]) # no rows with effect
            else:
                rows = [x.strip() for x in line.split(',')]
                for x in rows:
                    if x not in validrows:
                        raise IOError("Invalid row specification of %s in the following line: %s\nValid row labels are A to H." % (x, line))
                    if rows.count(x) != 1:
                        raise IOError("Row identifier of %s appears more than once in the following line: %s" % (x, line))
                sampledata[sample].append(rows)
    return (samplenames, sampledata, volume, dilution)


def AskOverwrite(filename):
    """If a file already exists, determines whether it should be overwritten.

    filename -> The name of the file which we are examining.
    If filename does not already exist, then this function returns
        True as we can write it without overwriting anything.
    If filename does already exist, then we ask the user to specify
        whether it should be overwritten. Returns True if the user's
        answer is yes, and False if the answer is no.
    """
    if os.path.isfile(filename):
        ans = None
        while True:
            ans = raw_input("File %s already exists. Should we overwrite it [Y, N]? " % filename).strip()
            if ans in ['Y', 'y']:
                print "The existing file will be overwritten."
                return True
            elif ans in ['N', 'n']:
                print "The existing file will not be overwritten."
                return False
            else:
                print "Invalid entry. Try again."
    else:
        return True


def main():
    """Main body of the script."""
    print "\nWelcome to the Bloom lab Reed-Muench calculator.\n"
    infile = None
    while not infile:
        infile = raw_input("Enter the name of the input file in text format: ").strip()
        if os.path.isfile(infile):
            break
        elif infile in ['Q', 'q']:
            print "Quitting."
            sys.exit()
        else:
            infile = None
            print "Failed to find the specified input file of %s. Try again to enter a valid file name, or enter Q to quit." % infile
    print "Reading input from the file %s." % infile
    (samplenames, sampledata, volume, dilution) = ParseInput(infile)
    print "Read data for %d samples." % len(sampledata)
    titers = {}
    for (sample, data) in sampledata.iteritems():
        titers[sample] = Titer(data, volume, dilution)
    print "\nHere are the computed titers in TCID50 per ul:"
    for sample in samplenames:
        print "%s: %.3f" % (sample, titers[sample])
    (base, ext) = os.path.splitext(infile)
    outfile = '%s-titers.txt' % base
    print "\nNow we will write these titers to the output file %s." % outfile
    if AskOverwrite(outfile):
        f = open(outfile, 'w')
        f.write("Here are the computed titers in TCID50 per ul.\n")
        for sample in samplenames:
            f.write("%s:\t%.3f\n" % (sample, titers[sample]))
        f.close()

# run the script in a try - except loop so the terminal doesn't close before the error can be diagnosed in windows.
try: 
    main() # run the script.
except Exception, e:
    print "ERROR."
    print e
finally:
    raw_input('Press any key to close the program.')

"""Module written to implement aspects of the Reed-Muench formula.

Jesse Bloom, 2012."""


import math


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


"""
Description of program
"""
import numpy as np
from astropy.table import Table as table
from astropy.io import fits
import matplotlib.pyplot as plt

# ==============================================================================
# Define variables
# ==============================================================================
# luminosity class determines size
lumorder = [-100, 0, 100]
lumfmt = {100: 10, 0: 5, -100: 2}
lumstr = {100: 'Giants', 0: 'Main Sequence', -100: 'Subdwarfs'}
# spectral subclass determines marker
subtypeorder = [0, 2, 4, 6, 8]
sptsubfmt = {0: 'o', 1: 'o', 2: 's', 3: 's', 4: 'd', 5: 'd',
             6: '^', 7: '^', 8: 'v', 9: 'v'}
sptsubstr = {0: r'$X0-X2$', 2: r'$X2-X4$', 4: r'$X4-X6$', 6: r'$X6-X8$',
             8: r'$X8+$'}
sptsubalp = {0: 0.40, 1: 0.40, 2: 0.55, 3: 0.55, 4: 0.70, 5: 0.70,
             6: 0.85, 7: 0.85, 8: 1.00, 9: 1.00}
# spectral class determines colour
spts_inorder = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T', 'Y']
sptclassfmt = {'O': 'red', 'B': 'crimson', 'A': 'darkorange',
               'F': 'darkgreen', 'G': 'lawngreen',
               'K': 'palegreen', 'M': 'blue', 'L': 'saddlebrown',
               'T': 'silver', 'Y': 'gray'}


# ==============================================================================
# Define functions
# ==============================================================================
def simbad_spt_to_num_spt(sptcol):
    """
    Converts SIMBAD spectral type to numerical spectral type (M0.0 = 0.0) and
    luminosity class and provides any binary flags

    :param sptcol: list or array containing string of SIMBAD spectral types

    :return: numerical spectral type
                - M0.0 --> +00.0
                - L4.5 --> +14.5
                - G0.0 --> -20.0
                - K7.0 --> -03.0

             luminosity class, array of float:
                - giant star (I, II, III, IV...) = 100
                - main sequence star (V) = 0
                - sub dwarf star (sd) = -100

             binary flag, array of boolean (whether flagged as binary)
                - binary is identified by a '+' in the string
                - the numerical spectral type and luminosity class will
                  only be calculated for the first spectral type found

    luminosity classes and spectral type information from:
        http://simbad.u-strasbg.fr/simbad/sim-display?data=sptypes

    """
    # check sptcol for length and string
    if not hasattr(sptcol, '__len__') or type(sptcol) == str:
        if type(sptcol) == str:
            sptcol = [sptcol]
        else:
            sptcol = [str(sptcol)]
    # check that all elemets are strings
    for si in range(len(sptcol)):
        if type(sptcol[si]) != str:
            sptcol[si] = str(si)
    # get spcs
    spckey, spcval, spckey1, spcval1 = get_spcs()
    # get lum classes
    lcckey, lccval = get_lum_classes()
    # get cont
    cont = get_cont()
    # don't edit original spt col
    strspts = np.array(sptcol)
    numspts = np.repeat(np.nan, len(sptcol))
    lumclass = np.repeat(np.nan, len(sptcol))
    binary = np.zeros(len(sptcol), dtype=bool)
    # loop round the spt rows
    for row in range(len(sptcol)):
        rawspt = strspts[row].replace(' ', '').upper()
        # if NaN then skip
        if len(rawspt) == 0:
            continue
        # remove the cont
        for crow in cont:
            rawspt = rawspt.replace(crow.upper(), '')
        # deal with binaries (i.e. a +) assume spectral type is primary
        # but flag as binary
        if '+' in rawspt:
            rawspt = rawspt.split('+')[0]
            binary[row] = True
        # deal with '...' (substrings that need to be split)
        for char in ['...']:
            if char in rawspt:
                rawspt = rawspt.split(char)[0]
        # strip and save luminosity class
        for lrow in range(len(lcckey)):
            if lcckey[lrow].upper() in rawspt:
                rawspt = rawspt.replace(lcckey[lrow].upper(), '')
                lumclass[row] = lccval[lrow]
        # remove any characters that aren't numbers or in spckey or spckey1
        for char in range(len(rawspt)):
            cond1 = rawspt[char].isdigit()
            cond2 = rawspt[char] in spckey
            cond3 = rawspt[char] in ['.', '/', '-']
            if (not cond1) and (not cond2) and (not cond3):
                rawspt = rawspt.replace(rawspt[char], ' ')
        # finally we should be able to extract a spectral type
        found = False
        # first need to take out annoying X-Y and X/Y --> assume these are the
        # second spectral type (i.e. Y)
        for srow in range(len(spckey1)):
            sp0, sp1, sp2 = spckey1[srow].upper()
            if sp0 in rawspt and sp1 in rawspt and sp2 in rawspt:
                rawspt = rawspt.split(sp1)[-1]
                rawspt = rawspt.replace(sp2.upper(), '')
                # if we are only given a letter assume it is X0.0
                if rawspt.replace(' ', '') == '':
                    numspts[row] = float(spcval1[srow])
                else:
                    numspts[row] = float(rawspt) + float(spcval1[srow])
                found = True
        # now remove all "/" and "-" left (should not be in here after last
        # step)
        if '/' in rawspt:
            rawspt = rawspt.split('/')[0]
        if '-' in rawspt:
            rawspt = rawspt.split('-')[0]
        # if we haven't yet found spectral try finding it for signal objects
        if found is False:
            for srow in range(len(spckey)):
                if spckey[srow].upper() in rawspt:
                    rawspt = rawspt.split(spckey[srow].upper())[-1]
                    # may still be other spt chars in
                    if np.sum(np.in1d(list(rawspt), spckey)) > 0:
                        continue
                    # if we are only given a letter assume it is X0.0
                    if rawspt.replace(' ', '') == '':
                        numspts[row] = float(spcval[srow])
                    else:
                        numspts[row] = float(rawspt) + float(spcval[srow])
    # finally return the numerical spectral types and the luminosity class
    return numspts, lumclass, binary


def get_spcs():
    """
    Defines the SIMBAD Spectral classes and then sorts them largest to
    smallest

    from http://simbad.u-strasbg.fr/simbad/sim-display?data=sptypes
    :return:
    """
    # define string spectral type convertion i.e. M0.0 --> 0.0
    spckey = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T', 'Y']
    spcval = [-60, -50, -40, -30, -20, -10, 0, 10, 20, 30]
    # define uncertain spts as the lower spectral type
    spckey1 = ['o-b', 'o/b', 'b-a', 'b/a', 'a-f', 'a/f', 'f-g', 'f/g',
               'g-k', 'g/k', 'k-m', 'k/m', 'm-l', 'm/l', 'l-t', 'l/t', 't-y',
               't/y']
    spcval1 = [-50, -50, -40, -40, -30, -30, -20, -20, -10, -10, 0, 0, 10, 10,
               20, 20, 30, 30]
    spckey, spcval = np.array(spckey), np.array(spcval)
    spckey1, spcval1 = np.array(spckey1), np.array(spcval1)
    mask = sort_str_array_by_length(spckey)
    mask1 = sort_str_array_by_length(spckey1)
    return spckey[mask], spcval[mask], spckey1[mask1], spcval1[mask1]


def get_lum_classes():
    """
    Defines the SIMBAD luminosity classes and then sorts them largest to
    smallest

    from http://simbad.u-strasbg.fr/simbad/sim-display?data=sptypes
    :return:
    """
    # define luminosity classes
    lcckey = ['Ia0', 'II-III', 'VI', 'Iab-Ib', 'IIb-IIIa', 'Ia-0', 'IIIa', 'I',
              'Iab/b', 'III/I', 'Ia0-Ia', 'III', 'I-II', 'Iab/II', 'III-IIIa',
              'Ia-0/Ia', 'IIIb', 'I/II', 'Ib/III', 'III-IIIb', 'Ia', 'III-IV',
              'I-III', 'Ib-IIIa', 'III/III-IV', 'Ia-ab', 'IVa', 'I/III', 'Ib/V',
              'III/IV', 'Ia/Iab', 'IV', 'III-V', 'II-IIb', 'IV-V/V', 'Iab',
              'IVb', 'III/V', 'II/II-III', 'V/III', 'Iab/Ib', 'IVa/V', '0-Ia',
              'II/III', 'c', 'Ib', 'IV-V', 'I/V', 'II-III/III', 'Ib-II',
              'IV/V', 'Ia/ab', 'II-IV', 'd', 'Ib/II', 'Va', 'Ia-Iab', 'II/IV',
              'sd', 'IIa', 'V', 'Ia-Ib', 'II/V', 'esd ', 'II', 'Vb', 'Ia/Ib',
              'IIab-IIb', 'usd', 'IIb ', 'V-VI', 'Iab-b', 'IIb-III']
    lccval = [100, 100, -100, 100, 100, 100, 100, 100,
              100, 100, 100, 100, 100, 100, 100,
              100, 100, 100, 100, 100, 100, 100,
              100, 100, 100, 100, 100, 100, 100,
              100, 100, 100, 100, 100, 100, 100,
              100, 100, 100, 0, 100, 100, 100, 100,
              100, 100, 100, 100, 100, 100,
              100, 100, 100, 0, 100, 0, 100, 100,
              -100, 100, 0, 100, 100, -100, 100, 0, 100,
              100, -100, 100, 0, 100, 100]
    # additional not in SIMBAD ref list
    lcckey += ['Vk', 'V(k)']
    lccval += [0, 0]

    lcckey, lccval = np.array(lcckey), np.array(lccval)
    mask = sort_str_array_by_length(lcckey)
    return lcckey[mask], lccval[mask]


def get_cont():
    """
    Defines the SIMBAD spectral types currently deemed contamination and
    then sorts them largest to smallest
    :return:

    from http://simbad.u-strasbg.fr/simbad/sim-display?data=sptypes
    """
    # remove these as contamination
    cont = ['OB+', 'OB', 'OB-', 'Of', 'Of*', 'Of+', 'O(f)', 'O(f+)', 'O((f*))',
            'O((f+))', 'O((f))', 'WN', 'WNE', 'WR', 'WC', 'WC+WN', 'R', 'N',
            'CH', 'C', 'DA', 'WD', 'DB', 'DO', 'DC', 'PG1159', 'DQ', 'DZ',
            'DX', 'SN', 'SN.I', 'SN.Ia', 'SN.Ia/c', 'SN.Ib', 'SN.Ib', 'SN.II',
            'SN.II/Ib', 'SN.II/Ic', 'SN.II/IIb', 'SN.IIn', 'SN.IIb', 'SN.IIL',
            'SN.IIP']
    # additional not in SIMBAD ref list
    cont += ['Fe', 'Am+', 'h', '_CN0.5']

    cont = np.array(cont)
    mask = sort_str_array_by_length(cont)
    return cont[mask]


def sort_str_array_by_length(array, first='longest'):
    """
    Takes a list or array of strings (array) and returns a mask of the sorting
    order (sorted by the length of each string) if first=shortest then shortest
    first else longest first
    :param array: list or array of string objects (or any object with length)
    :param first: if shortest returned mask has shortest first, else longest
    :return: sorting mask (longest to shortest or shortest to longest string)
    """
    lengths = []
    for row in array:
        lengths = np.append(lengths, len(row))
    if first == 'shortest':
        return np.argsort(lengths)
    else:
        return np.argsort(lengths)[::-1]


def get_spt_class_subclass(x):
    """
    Turn a numerical spectral type into a class and subclass

    i.e. 0.0 --> M0.0    -22.5 --> F7.5       22.5 --> T2.5

    :param x: float, numerical spectral type
    :return:
    """
    # get spectral type conversion information
    spckey, spcval, spckey1, spcval1 = get_spcs()
    # convert to dict
    spc = dict(zip(spcval, spckey))
    # find the nearest spectral type class to x
    nclass = np.floor(x / 10.0) * 10
    # use the spc dictionary to select this spectral class string
    sclass = spc[nclass]
    # spectral subclass is just the remainder
    ssubclass = x - nclass
    # return spectral class and subclass
    return sclass, ssubclass

# ==============================================================================
# Start of code
# ==============================================================================
if __name__ == '__main__':
    # Test 1
    test1 = 'M4.0V'
    spt, lumt, b = simbad_spt_to_num_spt(test1)
    args = [test1, spt[0], lumt[0], b[0]]
    print 'String = {0}\tSpt = {1}\tLtype = {2}\tbinary={3}'.format(*args)
    # Test 2
    test2 = 'K7.5II'
    spt, lumt, b = simbad_spt_to_num_spt(test2)
    args = [test2, spt[0], lumt[0], b[0]]
    print 'String = {0}\tSpt = {1}\tLtype = {2}\tbinary={3}'.format(*args)
    # Test 3
    test3 = 'L4.5 esd + T8.0'
    spt, lumt, b = simbad_spt_to_num_spt(test3)
    args = [test3, spt[0], lumt[0], b[0]]
    print 'String = {0}\tSpt = {1}\tLtype = {2}\tbinary={3}'.format(*args)
    # Test 4
    test4 = 4
    spt, lumt, b = simbad_spt_to_num_spt(test4)
    args = [test4, spt[0], lumt[0], b[0]]
    print 'String = {0}\tSpt = {1}\tLtype = {2}\tbinary={3}'.format(*args)
# ------------------------------------------------------------------------------

# ==============================================================================
# End of code
# ==============================================================================

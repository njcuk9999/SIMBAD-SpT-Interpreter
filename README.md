# SIMBAD-SpT-Interpreter

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

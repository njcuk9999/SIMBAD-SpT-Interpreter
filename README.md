# SIMBAD-SpT-Interpreter

```python
simbad_spt_to_num_spt(sptcol)
```

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

### Example code

```python
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
```

db.CH2CO = {}
db.CH2CO.atomicConstituents = {C=2,H=2,O=1,}
db.CH2CO.charge = 0
db.CH2CO.M = {
   value = 42.036680e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CO.gamma = {
   value = 1.1908e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      2.13583630E+00,
      1.81188721E-02,
     -1.73947474E-05,
      9.34397568E-09,
     -2.01457615E-12,
     -7.04291804E+03,
      1.22156480E+01,
   },
   segment1 = {
      0,
      0,
      4.51129732E+00,
      9.00359745E-03,
     -4.16939635E-06,
      9.23345882E-10,
     -7.94838201E-14,
     -7.55105311E+03,
      6.32247205E-01,
   }
}
db.CH2CO.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      3.549598e+04,
     -4.063063e+02,
      3.718922e+00,
      1.583502e-02,
     -1.726196e-05,
      1.157377e-08,
     -3.305843e-12,
     -5.209993e+03,
      3.839604e+00,
   },
   segment1 = {
      2.013565e+06,
     -8.200887e+03,
      1.759694e+01,
     -1.464545e-03,
      2.695887e-07,
     -2.665675e-11,
      1.094205e-15,
      4.177777e+04,
     -8.725804e+01
   }
}

db.CH2CO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.193913326467e+01,
      B = 2.737948202373e+00,
      C = -1.862634718958e-01,
      D = 4.621015645441e-03,
   }
}
db.CH2CO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.799368137710e+01,
      B = 3.225542866669e+00,
      C = -1.223022955920e-01,
      D = -2.848863216417e-03,
   }
}

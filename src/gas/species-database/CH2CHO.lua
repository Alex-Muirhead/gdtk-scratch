db.CH2CHO = {}
db.CH2CHO.atomicConstituents = {C=2,H=3,O=1,}
db.CH2CHO.charge = 0
db.CH2CHO.M = {
   value = 43.044620e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2CHO.gamma = {
   value = 1.1776e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2CHO.sigma = {
   value = 3.970,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.epsilon = {
   value = 436.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2CHO.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {300.0, 1000.0, 5000.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      0.03409062E+02,
      0.10738574E-01,
      0.01891492E-04,
     -0.07158583E-07,
      0.02867385E-10,
      0.15214766E+04,
      0.09558290E+02,
   },
   segment1 = {
      0,
      0,
      0.05975670E+02,
      0.08130591E-01,
     -0.02743624E-04,
      0.04070304E-08,
     -0.02176017E-12,
      0.04903218E+04,
     -0.05045251E+02,
   }
}
-- Same as GRIMech
db.CH2CHO.ceaThermoCoeffs = db.CH2CHO.grimechThermoCoeffs

db.CH2CHO.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -2.720702096187e+01,
      B = 4.986894409273e+00,
      C = -5.038737864236e-01,
      D = 1.949669378655e-02,
   }
}
db.CH2CHO.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 300.000,
      T_upper = 5000.000,
      A = -2.577625433343e+01,
      B = 6.373198044217e+00,
      C = -5.427655805016e-01,
      D = 1.604699040543e-02,
   }
}

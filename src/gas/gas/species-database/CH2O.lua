db.CH2O = {}
db.CH2O.atomicConstituents = {C=1,H=2,O=1}
db.CH2O.charge = 0
db.CH2O.M = {
   value = 30.025980e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
db.CH2O.gamma = {
   value = 1.3065e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
db.CH2O.sigma = {
   value = 3.590,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.epsilon = {
   value = 498.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH2O.Lewis = {
   value = 1.329
}
db.CH2O.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2, 
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0.000000000e+00, 
      0.000000000e+00, 
      4.793723150e+00, 
     -9.908333690e-03, 
      3.732200080e-05, 
     -3.792852610e-08, 
      1.317726520e-11, 
     -1.430895670e+04, 
      6.028129000e-01, 
  },
  segment1 = {
     0.000000000e+00, 
     0.000000000e+00, 
     1.760690080e+00, 
     9.200000820e-03, 
    -4.422588130e-06, 
     1.006412120e-09, 
    -8.838556400e-14, 
    -1.399583230e+04, 
     1.365632300e+01, 
   }
}
-- Same as GRIMech
db.CH2O.ceaThermoCoeffs = db.CH2O.grimechThermoCoeffs

db.CH2O.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.944420616305e+01,
      B = 1.583460638757e+00,
      C = -1.418900703767e-02,
      D = -3.702993576843e-03,
   }
}
db.CH2O.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -3.458688239117e-01,
      B = -4.719241554972e+00,
      C = 1.048457032022e+00,
      D = -5.935855918886e-02,
   }
}

db.CH2O.Hf = {
   value = -108580.0,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
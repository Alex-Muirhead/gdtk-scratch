db.CH3 = {}
db.CH3.atomicConstituents = {C=1,H=3,}
db.CH3.charge = 0
db.CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
db.CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
db.CH3.sigma = {
   value = 3.800,
   units = 'Angstrom',
   description = 'Lennard-Jones potential distance',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.epsilon = {
   value = 144.000,
   units = 'K',
   description = 'Lennard-Jones potential well depth.',
   reference = 'GRI-Mech 3.0 transport file.'
}
db.CH3.Lewis = {
   value = 1.049
}
db.CH3.grimechThermoCoeffs = {
   notes = 'data from GRIMECH 3.0',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 3500.0},
   T_blend_ranges = {400.0},
   segment0 = {
      0,
      0,
      3.67359040E+00,
      2.01095175E-03,
      5.73021856E-06,
     -6.87117425E-09,
      2.54385734E-12,
      1.64449988E+04,
      1.60456433E+00,
   },
   segment1 = {
      0,
      0,
      2.28571772E+00,
      7.23990037E-03,
     -2.98714348E-06,
      5.95684644E-10,
     -4.67154394E-14,
      1.67755843E+04,
      8.48007179E+00,
   }
}
db.CH3.ceaThermoCoeffs = {
   notes = 'NASA/TP—2002-211556',
   nsegments = 2,
   T_break_points = {200.0, 1000.0, 6000.0},
   T_blend_ranges = {400.0},
   segment0 = {
     -2.876189e+04,
      5.093269e+02,
      2.002144e-01,
      1.363606e-02,
     -1.433989e-05,
      1.013557e-08,
     -3.027332e-12,
      1.408272e+04,
      2.022773e+01,
   },
   segment1 = {
      2.760803e+06,
     -9.336531e+03,
      1.487730e+01,
     -1.439430e-03,
      2.444478e-07,
      -2.224556e-11,
      8.395066e-16,
      7.481809e+04,
     -7.919682e+01
   }
}

db.CH3.chemkinViscosity = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -2.270870995825e+01,
      B = 3.690876383840e+00,
      C = -4.016022429639e-01,
      D = 1.764922953865e-02,
   }
}
db.CH3.chemkinThermCond = {
   notes = 'Generated by species-generator.py',
   nsegments = 1, 
   segment0 ={
      T_lower = 200.000,
      T_upper = 3500.000,
      A = -1.243228188431e-01,
      B = -3.501395522712e+00,
      C = 7.418831266650e-01,
      D = -3.973086278935e-02,
   }
}

db.CH3.Hf = {
   value = 146658.04,
   units = 'J/mol',
   description = 'Molar Heat of Formation at 298.15K',
   reference = 'CEA2::thermo.inp'
}
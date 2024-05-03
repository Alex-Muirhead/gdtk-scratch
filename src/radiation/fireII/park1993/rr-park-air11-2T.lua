species = {[0]='N2', [1]='O2', [2]='N', [3]='O', [4]='NO', [5]='N2+', [6]='O2+', [7]='N+', [8]='O+', [9]='NO+', [10]='e-', }
config = {
  tempLimits = {lower=300.000000, upper=30000.000000},
  odeStep = {method='rkf', errTol=1.000000e-09},
  tightTempCoupling = true,
  maxSubcycles = 10000,
  maxAttempts = 4
}

reaction = {}
reaction[1] = {
  equation = "N2 + M <=> N + N + M",
  type = "anonymous_collider",
  frc = {model='Park', A=7.000000000000e+15, n=-1.600000, C=1.132000000000e+05, s=0.500000, mode=0 },
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0,},
  reacCoeffs = { 1,},
  prodIdx = { 2,},
  prodCoeffs = { 2,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=4.285700e+00,
    [3]=4.285700e+00,
    [4]=1.000000e+00,
    [5]=1.000000e+00,
    [6]=1.000000e+00,
    [7]=4.285700e+00,
    [8]=4.285700e+00,
    [9]=1.000000e+00,
    [10]=1.714286e+03,
  },
}

reaction[2] = {
  equation = "O2 + M <=> O + O + M",
  type = "anonymous_collider",
  frc = {model='Park', A=2.000000000000e+15, n=-1.500000, C=5.950000000000e+04, s=0.500000, mode=0 },
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 1,},
  reacCoeffs = { 1,},
  prodIdx = { 3,},
  prodCoeffs = { 2,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=5.000000e+00,
    [3]=5.000000e+00,
    [4]=1.000000e+00,
    [5]=1.000000e+00,
    [6]=1.000000e+00,
    [7]=5.000000e+00,
    [8]=5.000000e+00,
    [9]=1.000000e+00,
  },
}

reaction[3] = {
  equation = "NO + M <=> N + O + M",
  type = "anonymous_collider",
  frc = {model='Park', A=5.000000000000e+09, n=0.000000, C=7.550000000000e+04, s=0.500000, mode=0 },
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 4,},
  reacCoeffs = { 1,},
  prodIdx = { 2, 3,},
  prodCoeffs = { 1, 1,},
  efficiencies = {
    [0]=1.000000e+00,
    [1]=1.000000e+00,
    [2]=2.000000e+01,
    [3]=2.000000e+01,
    [4]=1.000000e+00,
    [5]=1.000000e+00,
    [6]=1.000000e+00,
    [7]=2.000000e+01,
    [8]=2.000000e+01,
    [9]=1.000000e+00,
  },
}

reaction[4] = {
  equation = "NO + O <=> O2 + N",
  type = "elementary",
  frc = {model='Arrhenius', A=8.400000000000e+06, n=0.000000, C=1.945000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 4,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 2,},
  prodCoeffs = { 1, 1,},
}

reaction[5] = {
  equation = "N2 + O <=> NO + N",
  type = "elementary",
  frc = {model='Arrhenius', A=6.400000000000e+11, n=-1.000000, C=3.840000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 3,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 2, 4,},
  prodCoeffs = { 1, 1,},
}

reaction[6] = {
  equation = "N + O <=> NO+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=8.800000000000e+02, n=1.000000, C=3.190000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2, 3,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 9, 10,},
  prodCoeffs = { 1, 1,},
}

reaction[7] = {
  equation = "O + O <=> O2+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=7.100000000000e-04, n=2.700000, C=8.060000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3,},
  reacCoeffs = { 2,},
  prodIdx = { 6, 10,},
  prodCoeffs = { 1, 1,},
}

reaction[8] = {
  equation = "N + N <=> N2+ + e-",
  type = "elementary",
  frc = {model='Arrhenius', A=4.400000000000e+01, n=1.500000, C=6.750000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2,},
  reacCoeffs = { 2,},
  prodIdx = { 5, 10,},
  prodCoeffs = { 1, 1,},
}

reaction[9] = {
  equation = "NO+ + O <=> N+ + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+06, n=0.500000, C=7.720000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 9,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 7,},
  prodCoeffs = { 1, 1,},
}

reaction[10] = {
  equation = "N+ + N2 <=> N2+ + N",
  type = "elementary",
  frc = {model='Arrhenius', A=1.000000000000e+06, n=0.500000, C=1.220000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 7,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 2, 5,},
  prodCoeffs = { 1, 1,},
}

reaction[11] = {
  equation = "O2+ + N <=> N+ + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=8.700000000000e+07, n=0.140000, C=2.860000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2, 6,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 7,},
  prodCoeffs = { 1, 1,},
}

reaction[12] = {
  equation = "O+ + NO <=> N+ + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=1.400000000000e-01, n=1.900000, C=2.660000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 4, 8,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 7,},
  prodCoeffs = { 1, 1,},
}

reaction[13] = {
  equation = "O2+ + N2 <=> N2+ + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=9.900000000000e+06, n=0.000000, C=4.070000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 6,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 5,},
  prodCoeffs = { 1, 1,},
}

reaction[14] = {
  equation = "O2+ + O <=> O+ + O2",
  type = "elementary",
  frc = {model='Arrhenius', A=4.000000000000e+06, n=-0.090000, C=1.800000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 6,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 1, 8,},
  prodCoeffs = { 1, 1,},
}

reaction[15] = {
  equation = "NO+ + N <=> O+ + N2",
  type = "elementary",
  frc = {model='Arrhenius', A=3.400000000000e+07, n=-1.080000, C=1.280000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2, 9,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 0, 8,},
  prodCoeffs = { 1, 1,},
}

reaction[16] = {
  equation = "NO+ + O2 <=> O2+ + NO",
  type = "elementary",
  frc = {model='Arrhenius', A=2.400000000000e+07, n=0.410000, C=3.260000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 1, 9,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 4, 6,},
  prodCoeffs = { 1, 1,},
}

reaction[17] = {
  equation = "NO+ + O <=> O2+ + N",
  type = "elementary",
  frc = {model='Arrhenius', A=7.200000000000e+06, n=0.290000, C=4.860000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 9,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 2, 6,},
  prodCoeffs = { 1, 1,},
}

reaction[18] = {
  equation = "O+ + N2 <=> N2+ + O",
  type = "elementary",
  frc = {model='Arrhenius', A=9.100000000000e+05, n=0.360000, C=2.280000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 0, 8,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 3, 5,},
  prodCoeffs = { 1, 1,},
}

reaction[19] = {
  equation = "NO+ + N <=> N2+ + O",
  type = "elementary",
  frc = {model='Arrhenius', A=7.200000000000e+07, n=0.000000, C=3.550000000000e+04, rctIndex=-1},
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2, 9,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 3, 5,},
  prodCoeffs = { 1, 1,},
}

reaction[20] = {
  equation = "O + e- <=> O+ + e- + e-",
  type = "elementary",
  frc = {model='Park', A=3.900000000000e+27, n=-3.780000, C=1.585000000000e+05, s=0.000000, mode=0 },
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 3, 10,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 8, 10,},
  prodCoeffs = { 1, 2,},
}

reaction[21] = {
  equation = "N + e- <=> N+ + e- + e-",
  type = "elementary",
  frc = {model='Park', A=2.500000000000e+28, n=-3.820000, C=1.686000000000e+05, s=0.000000, mode=0 },
  brc = {model='fromEqConst', rctIndex=-1},
  ec = {},
  reacIdx = { 2, 10,},
  reacCoeffs = { 1, 1,},
  prodIdx = { 7, 10,},
  prodCoeffs = { 1, 2,},
}


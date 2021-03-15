species = {'N2', 'O2', 'NO', 'N', 'O', }

db = {}
db.N2 = {}
db.N2.type = "molecule"
db.N2.molecule_type = "linear"
db.N2.M = 0.0280134
db.N2.charge = 0
db.N2.sigma = 3.621
db.N2.epsilon = 97.530
db.N2.Lewis = 1.152
db.N2.viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.62526577,
      C = -1640.7983,
      B = -31.779652,
      T_upper = 1000,
      T_lower = 200,
      D = 1.7454992,
    },
   segment1 = {
      A = 0.87395209,
      C = -173948.09,
      B = 561.52222,
      T_upper = 5000,
      T_lower = 1000,
      D = -0.39335958,
    },
   segment2 = {
      A = 0.88503551,
      C = -731290.61,
      B = 909.02171,
      T_upper = 15000,
      T_lower = 5000,
      D = -0.53503838,
    },
}
db.N2.thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.85439436,
      C = -12347.848,
      B = 105.73224,
      T_upper = 1000,
      T_lower = 200,
      D = 0.47793128,
    },
    segment1 =  {
      A = 0.88407146,
      C = -11429.64,
      B = 133.57293,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.24417019,
    },
    segment2 = {
      A = 2.4176185,
      C = 3105580.2,
      B = 8047.7749,
      T_upper = 15000,
      T_lower = 5000,
      D = -14.517761,
    },
}
db.N2.thermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      22103.71497,
      -381.846182,
      6.08273836,
      -0.00853091441,
      1.384646189e-05,
      -9.62579362e-09,
      2.519705809e-12,
      710.846086,
      -10.76003744,
   },
   segment1 = {
      587712.406,
      -2239.249073,
      6.06694922,
      -0.00061396855,
      1.491806679e-07,
      -1.923105485e-11,
      1.061954386e-15,
      12832.10415,
      -15.86640027,

   },
   segment2 =  {
      831013916,
      -642073.354,
      202.0264635,
      -0.03065092046,
      2.486903333e-06,
      -9.70595411e-11,
      1.437538881e-15,
      4938707.04,
      -1672.09974,
  },
}
db.O2 = {}
db.O2.M = 0.0319988
db.O2.type = "molecule"
db.O2.molecule_type = "linear"
db.O2.charge = 0
db.O2.sigma = 3.458
db.O2.epsilon = 107.400
db.O2.Lewis = 1.086
db.O2.viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.6091618,
      C = -599.74009,
      B = -52.244847,
      T_upper = 1000,
      T_lower = 200,
      D = 2.0410801,
   },
   segment1 = {
      A = 0.72216486,
      C = -57974.816,
      B = 175.50839,
      T_upper = 5000,
      T_lower = 1000,
      D = 1.0901044,
   },
   segment2 = {
      A = 0.73981127,
      C = -378331.68,
      B = 391.94906,
      T_upper = 15000,
      T_lower = 5000,
      D = 0.9093178,
   },
}
db.O2.thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.77229167,
      C = -5893.3377,
      B = 6.846321,
      T_upper = 1000,
      T_lower = 200,
      D = 1.2210365,
    },
   segment1 = {
      A = 0.90917351,
      C = -79650.171,
      B = 291.24182,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.064851631,
    },
   segment2 = {
      A = -1.1218262,
      C = 23295011,
      B = -19286.378,
      T_upper = 15000,
      T_lower = 5000,
      D = 20.342043,
    },
}
db.O2.thermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      -34255.6342,
      484.700097,
      1.119010961,
      0.00429388924,
      -6.83630052e-07,
      -2.0233727e-09,
      1.039040018e-12,
      -3391.45487,
      18.4969947,

   },
   segment1 = {
     -1037939.022,
      2344.830282,
      1.819732036,
      0.001267847582,
      -2.188067988e-07,
      2.053719572e-11,
      -8.19346705e-16,
      -16890.10929,
      17.38716506,
   },
   segment2 = {
      497529430,
      -286610.6874,
      66.9035225,
      -0.00616995902,
      3.016396027e-07,
      -7.4214166e-12,
      7.27817577e-17,
      2293554.027,
      -553.062161,
   }
}
db.NO = {}
db.NO.M = 0.0300061
db.NO.type = "molecule"
db.NO.molecule_type = "linear"
db.NO.charge = 0
db.NO.sigma = 3.621
db.NO.epsilon = 97.530
db.NO.Lewis = 1.0
db.NO.viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.60262029,
      C = -139.54524,
      B = -62.017783,
      T_upper = 1000,
      T_lower = 200,
      D = 2.0268332,
    },
   segment1 = {
      A = 0.7800905,
      C = -94847.722,
      B = 304.86891,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.52873381,
    },
   segment2 = {
      A = 0.80580582,
      C = -578792.1,
      B = 624.27878,
      T_upper = 15000,
      T_lower = 5000,
      D = 0.2651645,
    },
}
db.NO.thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      A = 0.95028758,
      C = -9989.4764,
      B = 76.667058,
      T_upper = 1000,
      T_lower = 200,
      D = -0.0062776717,
    },
   segment1 = {
      A = 0.86215238,
      C = -238564.66,
      B = 445.68223,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.46209876,
    },
   segment2 = {
      A = -1.0377865,
      C = 67451187,
      B = -34486.864,
      T_upper = 15000,
      T_lower = 5000,
      D = 20.928749,
    }
}
db.NO.thermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      -11439.16503,
      153.6467592,
      3.43146873,
      -0.002668592368,
      8.48139912e-06,
      -7.68511105e-09,
      2.386797655e-12,
      9098.21441,
      6.72872549,
   },
   segment1 = {
      223901.8716,
      -1289.651623,
      5.43393603,
      -0.00036560349,
      9.88096645e-08,
      -1.416076856e-11,
      9.38018462e-16,
      17503.17656,
      -8.50166909,
   },
   segment2 = {
     -957530354,
      591243.448,
      -138.4566826,
      0.01694339403,
      -1.007351096e-06,
      2.912584076e-11,
      -3.29510935e-16,
      -4677501.24,
      1242.081216,
  },
}
db.N = {}
db.N.M = 0.0140067
db.N.type = "atom"
db.N.charge = 0
db.N.sigma = 3.298 
db.N.epsilon = 71.400 
db.N.Lewis = 1.0
db.N.viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      A = 0.83724737,
      C = -174507.53,
      B = 439.9715,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.10365689,
    },
   segment1 = {
      A = 0.89986588,
      C = -1820047.8,
      B = 1411.2801,
      T_upper = 15000,
      T_lower = 5000,
      D = -0.55811716,
    },
}
db.N.thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      A = 0.83771661,
      C = -175784.46,
      B = 442.4327,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.89942915,
    },
   segment1 = {
      A = 0.9000171,
      C = -1826240.3,
      B = 1414.1175,
      T_upper = 15000,
      T_lower = 5000,
      D = 0.24048513,
    },
}
db.N.thermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      0,
      0,
      2.5,
      0,
      0,
      0,
      0,
      56104.6378,
      4.193905036,
  },
   segment1 = {
      88765.0138,
      -107.12315,
      2.362188287,
      0.0002916720081,
      -1.7295151e-07,
      4.01265788e-11,
      -2.677227571e-15,
      56973.5133,
      4.865231506,
   },
   segment2 =  {
      547518105,
      -310757.498,
      69.1678274,
      -0.00684798813,
      3.8275724e-07,
      -1.098367709e-11,
      1.277986024e-16,
      2550585.618,
      -584.8769753,
   },
}
db.O = {}
db.O.M = 0.0159994
db.O.type = "atom"
db.O.charge = 0
db.O.sigma = 2.750
db.O.epsilon = 80.000
db.O.Lewis = 0.712
db.O.viscosity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      A = 0.77269241,
      C = -58502.098,
      B = 83.842977,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.85100827,
   },
   segment1 = {
      A = 0.87669586,
      C = -1088456.6,
      B = 1015.842,
      T_upper = 15000,
      T_lower = 5000,
      D = -0.18001077,
   },
}
db.O.thermal_conductivity = {
   model = 'CEA',
   nsegments = 2,
   segment0 = {
      A = 0.77271664,
      C = -58580.966,
      B = 83.9891,
      T_upper = 5000,
      T_lower = 1000,
      D = 1.51799,
    },
   segment1 = {
      A = 0.87676666,
      C = -1090669,
      B = 1017.0744,
      T_upper = 15000,
      T_lower = 5000,
      D = 0.48644232,
    },
}
db.O.thermoCoeffs = {
   nsegments = 3,
   T_break_points = {200.0, 1000.0, 6000.0, 20000.0},
   T_blend_ranges = {400.0, 1000.0},
   segment0 = {
      -7953.6113,
      160.7177787,
      1.966226438,
      0.00101367031,
      -1.110415423e-06,
      6.5175075e-10,
      -1.584779251e-13,
      28403.62437,
      8.40424182,
  },
   segment1 = {
      261902.0262,
      -729.872203,
      3.31717727,
      -0.000428133436,
      1.036104594e-07,
      -9.43830433e-12,
      2.725038297e-16,
      33924.2806,
      -0.667958535,
  },
   segment2 = {
      177900426.4,
      -108232.8257,
      28.10778365,
      -0.002975232262,
      1.854997534e-07,
      -5.79623154e-12,
      7.191720164e-17,
      889094.263,
      -218.1728151,
  },
}
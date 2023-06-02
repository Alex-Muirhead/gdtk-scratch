# Step through a steady isentropic expansion,
# from stagnation condition to sonic condition.
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 isentropic-air-expansion2.py
#
# Python port, PJ, 2019-11-21
#                  2023-06-02 PyGasState flavour
#
import math
from gdtk.gas import GasModel, PyGasState

gmodel = GasModel('ideal-air-gas-model.lua')
gs = PyGasState(gmodel)
gs.p = 500e3 # Pa
gs.T = 300.0 # K
gs.update_thermo_from_pT()
# Compute enthalpy and entropy at stagnation conditions
h0 = gs.enthalpy
s0 = gs.entropy
# Set up for stepping process
dp = 1.0 # Pa, use 1 Pa as pressure step size
gs.p = gs.p - dp
mach_tgt = 1.0
# Begin stepping until Mach = mach_tgt
while True:
    gs.update_thermo_from_ps(s0)
    h1 = gs.enthalpy
    v1 = math.sqrt(2*(h0 - h1))
    m1 = v1/gs.a
    if m1 >= mach_tgt:
        print("Stopping at Mach=%g" % m1)
        break
    gs.p = gs.p - dp

print("Gas properties at sonic point are:")
print("p=%g T=%g" % (gs.p, gs.T))

job_title = "Fire II Inviscid Flow"
config.dimensions = 2
config.axisymmetric = false

setGasModel('ideal-air-gas-model.lua')

flow = FlowState:new{
  p=1000.0, T=500.0,
  velx=0.0, vely=0.0
};

flowDict = { initial = flow }
bcDict = {}

makeFluidBlocks(bcDict, flowDict)

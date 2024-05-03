job_title = "Fire II Inviscid Flow"
config.dimensions = 2
config.axisymmetric = true

nsp, nmodes, gmodel = setGasModel('gm-air11-2T.lua')
config.reacting = true
config.reactions_file = 'rr-park-air11-2T.lua'
config.energy_exchange_file = 'ee-park-air11-2T.lua'

-- 1636 Second Trajectory Point, at 71.04 km Altitude
Ri = 0.9347
T_inf = 210      -- K
T_wall= 810
rho_inf = 8.57e-5  -- kg/m^3
mass_fraction = {N2=0.767, O2=0.233}
u_inf = 11.31e3    -- m/s

-- Compute full inflow state from T,rho, and u
Q = GasState:new{gmodel}
Q.T = T_inf;
Q.rho = rho_inf;
Q.massf = mass_fraction
Q.T_modes = {T_inf}
gmodel:updateThermoFromRHOT(Q);
gmodel:updateSoundSpeed(Q)
p_inf = Q.p
M_inf = u_inf/Q.a
gmodel:updateTransCoeffs(Q);

initial = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, vely=0.0, massf=mass_fraction, T_modes={T_inf}}

flowDict = {
   initial=initial,
   inflow=inflow
}

bcDict = {
   inflow_sf=InFlowBC_ShockFitting:new{flowState=inflow},
   outflow=OutFlowBC_Simple:new{}
}
--
makeFluidBlocks(bcDict, flowDict)

-- Configuration options
N_flows = 12.0
N_solutions = 24.0
flow_time = Ri/u_inf
config.max_time = N_flows*flow_time
config.dt_plot = config.max_time/N_solutions
config.dt_init = 1.0e-16
config.cfl_value = 0.2
config.max_step = 2147483647
config.print_count=100
config.flux_calculator = "hanel"
config.gasdynamic_update_scheme = "moving_grid_2_stage"
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 1.0*flow_time
config.shock_fitting_scale_factor = 0.1
config.flowstate_limits_max_temp = 100000
config.interpolation_delay = 4.0*flow_time
config.allow_reconstruction_for_species = false

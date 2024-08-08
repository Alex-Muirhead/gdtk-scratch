# gas_files.mk
#
# We list the source files that are needed to build the gas models package.
#
SRC_DIR := $(GAS_DIR)/gas
GAS_MODEL_FILES := $(SRC_DIR)/package.d \
	$(SRC_DIR)/composite_gas.d \
	$(SRC_DIR)/gas_model.d \
	$(SRC_DIR)/gas_state.d \
	$(SRC_DIR)/init_gas_model.d \
	$(SRC_DIR)/ideal_gas.d \
	$(SRC_DIR)/ideal_helium.d \
	$(SRC_DIR)/cubic_gas.d \
	$(SRC_DIR)/cea_gas.d \
	$(SRC_DIR)/physical_constants.d \
	$(SRC_DIR)/therm_perf_gas.d \
	$(SRC_DIR)/therm_perf_gas_equil.d \
	$(SRC_DIR)/very_viscous_air.d \
	$(SRC_DIR)/uniform_lut.d \
	$(SRC_DIR)/uniform_lut_plus_ideal.d \
	$(SRC_DIR)/adaptive_lut_CEA.d \
	$(SRC_DIR)/ideal_gas_ab.d \
	$(SRC_DIR)/two_temperature_reacting_argon.d \
	$(SRC_DIR)/two_temperature_argon_plus_ideal.d \
	$(SRC_DIR)/ideal_dissociating_gas.d \
	$(SRC_DIR)/two_temperature_air.d \
	$(SRC_DIR)/two_temperature_nitrogen.d \
	$(SRC_DIR)/two_temperature_dissociating_nitrogen.d \
	$(SRC_DIR)/vib_specific_nitrogen.d \
	$(SRC_DIR)/vib_specific_co.d \
	$(SRC_DIR)/fuel_air_mix.d \
	$(SRC_DIR)/equilibrium_gas.d \
	$(SRC_DIR)/two_temperature_gasgiant.d

THERMO_FILES := $(SRC_DIR)/thermo/package.d \
	$(SRC_DIR)/thermo/cea_thermo_curves.d \
	$(SRC_DIR)/thermo/evt_eos.d \
	$(SRC_DIR)/thermo/perf_gas_mix_eos.d \
	$(SRC_DIR)/thermo/pvt_eos.d \
	$(SRC_DIR)/thermo/therm_perf_gas_mix_eos.d \
	$(SRC_DIR)/thermo/thermo_model.d \
	$(SRC_DIR)/thermo/therm_perf_gas_mix.d \
	$(SRC_DIR)/thermo/two_temperature_gas.d \
	$(SRC_DIR)/thermo/three_temperature_gas.d \
	$(SRC_DIR)/thermo/multi_temperature_gas.d \
	$(SRC_DIR)/thermo/energy_modes.d

DIFFUSION_FILES := $(SRC_DIR)/diffusion/package.d \
	$(SRC_DIR)/diffusion/cea_therm_cond.d \
	$(SRC_DIR)/diffusion/cea_viscosity.d \
	$(SRC_DIR)/diffusion/chemkin_therm_cond.d \
	$(SRC_DIR)/diffusion/chemkin_viscosity.d \
	$(SRC_DIR)/diffusion/gas_mixtures.d \
	$(SRC_DIR)/diffusion/sutherland_therm_cond.d \
	$(SRC_DIR)/diffusion/sutherland_viscosity.d \
	$(SRC_DIR)/diffusion/therm_cond.d \
	$(SRC_DIR)/diffusion/transport_properties_model.d \
	$(SRC_DIR)/diffusion/two_temperature_trans_props.d \
	$(SRC_DIR)/diffusion/multi_temperature_trans_props.d \
	$(SRC_DIR)/diffusion/three_temperature_trans_props.d \
	$(SRC_DIR)/diffusion/viscosity.d \
	$(SRC_DIR)/diffusion/wilke_mixing_therm_cond.d \
	$(SRC_DIR)/diffusion/wilke_mixing_viscosity.d \
	$(SRC_DIR)/diffusion/gasgiant_transport_properties.d \
	$(SRC_DIR)/diffusion/binary_diffusion_coefficients.d \
	$(SRC_DIR)/diffusion/rps_diffusion_coefficients.d

GAS_FILES := $(GAS_MODEL_FILES) $(THERMO_FILES) $(DIFFUSION_FILES)

GAS_LUA_FILES := $(SRC_DIR)/luagas_model.d

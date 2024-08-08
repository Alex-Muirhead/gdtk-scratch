SRC_DIR := $(KINETICS_DIR)/kinetics
KINETICS_FILES := $(SRC_DIR)/package.d \
	$(SRC_DIR)/thermochemical_reactor.d \
	$(SRC_DIR)/init_thermochemical_reactor.d \
	$(SRC_DIR)/chemistry_update.d \
	$(SRC_DIR)/energy_exchange_mechanism.d \
	$(SRC_DIR)/energy_exchange_system.d \
	$(SRC_DIR)/equilibrium_update.d \
	$(SRC_DIR)/ideal_dissociating_gas_kinetics.d \
	$(SRC_DIR)/fuel_air_mix_kinetics.d \
	$(SRC_DIR)/powers_aslam_kinetics.d \
	$(SRC_DIR)/yee_kotov_kinetics.d \
	$(SRC_DIR)/rate_constant.d \
	$(SRC_DIR)/reaction.d \
	$(SRC_DIR)/reaction_mechanism.d \
	$(SRC_DIR)/relaxation_time.d \
	$(SRC_DIR)/exchange_cross_section.d \
	$(SRC_DIR)/exchange_chemistry_coupling.d \
	$(SRC_DIR)/multi_temperature_thermochemical_reactor.d \
	$(SRC_DIR)/two_temperature_air_kinetics.d \
	$(SRC_DIR)/two_temperature_argon_kinetics.d \
	$(SRC_DIR)/two_temperature_argon_with_ideal_gas.d \
	$(SRC_DIR)/two_temperature_nitrogen_kinetics.d \
	$(SRC_DIR)/two_temperature_dissociating_nitrogen_kinetics.d \
	$(SRC_DIR)/vib_specific_nitrogen_kinetics.d \
	$(SRC_DIR)/vib_specific_co_kinetics.d \
	$(SRC_DIR)/two_temperature_gasgiant_kinetics.d

KINETICS_LUA_FILES := $(SRC_DIR)/luathermochemical_reactor.d \
	$(SRC_DIR)/luachemistry_update.d \
	$(SRC_DIR)/luaequilibrium_calculator.d \
	$(SRC_DIR)/luareaction_mechanism.d \
	$(SRC_DIR)/luatwo_temperature_air_kinetics.d \
	$(SRC_DIR)/luavib_specific_nitrogen_kinetics.d

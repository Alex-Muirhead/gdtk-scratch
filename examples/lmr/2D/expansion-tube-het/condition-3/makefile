LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep"
	@echo "   run"
	@echo "   post"

# Prepare a simulation with a 1T gas model.
prep: gas1T grid init-sim

# We have our own specification of carbon dioxide species plus helium
# and we want to use it with Rowan's 2T carbon dioxide kinetics files.
gas2T: co2-he-2T.lua
	prep-gas co2-he-2T.lua co2-he-2T.gas
	cp $(DGD_REPO)/examples/kinetics/co2-chemistry-neutrals-2T/co2-5sp-2T-chemistry.lua .
	prep-chem co2-he-2T.gas co2-5sp-2T-chemistry.lua co2-he-2T-chemistry.chem
	cp $(DGD_REPO)/examples/kinetics/co2-chemistry-neutrals-2T/co2-mixture-energy-exchange.lua .
	prep-kinetics co2-he-2T.gas co2-he-2T-chemistry.chem co2-mixture-energy-exchange.lua co2-mixture-energy-exchange.kin

# For a 1T simulation, we use the reactions file with Park's rates, as prepared by Rowan.
gas1T: co2-he-1T.lua
	prep-gas co2-he-1T.lua co2-he-1T.gas
	cp $(DGD_REPO)/examples/kinetics/co2-chemistry-neutrals-1T/park-co2-neutrals-5sp-1T-chemistry.lua .
	prep-chem co2-he-1T.gas park-co2-neutrals-5sp-1T-chemistry.lua co2-he-1T-chemistry.chem

grid: grid.lua
	lmr prep-grid --job=grid.lua

init-sim: het.lua
	lmr prep-sim --job=het.lua

run:
	# lmr run
	mpirun -np 14 --oversubscribe lmr-mpi-run

post: vtk xtdiagram histories

vtk:
	lmr snapshot2vtk --all

xtdiagram:
	./lmr-xtdata.rb --output=xtlogp.data --log10
	gnuplot < contour.gnuplot

histories:
	gnuplot < plot-p-history.gnuplot
	gnuplot < plot-T-history.gnuplot
	gnuplot < plot-rv2-history.gnuplot
	gnuplot < plot-velx-history.gnuplot
	gnuplot < plot-species-history.gnuplot

clean:
	- rm -rf $(LMR_DIR)
	- rm *.gas *.chem *.kin

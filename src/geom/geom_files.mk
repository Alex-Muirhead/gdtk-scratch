SRC_DIR := $(GEOM_DIR)/geom

GEOM_D_FILES := \
	$(SRC_DIR)/elements/nomenclature.d \
	$(SRC_DIR)/elements/package.d \
	$(SRC_DIR)/elements/projection.d \
	$(SRC_DIR)/elements/properties.d \
	$(SRC_DIR)/elements/vector3.d \
	\
	$(SRC_DIR)/gpath/arc.d \
	$(SRC_DIR)/gpath/bezier.d \
	$(SRC_DIR)/gpath/gpath_utils.d \
	$(SRC_DIR)/gpath/helix.d \
	$(SRC_DIR)/gpath/line.d \
	$(SRC_DIR)/gpath/modifiedpath.d \
	$(SRC_DIR)/gpath/nurbs.d \
	$(SRC_DIR)/gpath/package.d \
	$(SRC_DIR)/gpath/path.d \
	$(SRC_DIR)/gpath/polyline.d \
	$(SRC_DIR)/gpath/polynomial.d \
	$(SRC_DIR)/gpath/svgpath.d \
	$(SRC_DIR)/gpath/xspline.d \
	$(SRC_DIR)/gpath/xsplinelsq.d \
	\
	$(SRC_DIR)/grid/grid.d \
	$(SRC_DIR)/grid/package.d \
	$(SRC_DIR)/grid/paver.d \
	$(SRC_DIR)/grid/paver2d.d \
	$(SRC_DIR)/grid/sgrid.d \
	$(SRC_DIR)/grid/usgrid.d \
	\
	$(SRC_DIR)/misc/kdtree.d \
	$(SRC_DIR)/misc/nurbs_utils.d \
	$(SRC_DIR)/misc/package.d \
	$(SRC_DIR)/misc/sketch.d \
	$(SRC_DIR)/misc/svg.d \
	$(SRC_DIR)/misc/univariatefunctions.d \
	\
	$(SRC_DIR)/surface/aopatch.d \
	$(SRC_DIR)/surface/bezierpatch.d \
	$(SRC_DIR)/surface/beziertrianglepatch.d \
	$(SRC_DIR)/surface/channelpatch.d \
	$(SRC_DIR)/surface/controlpointpatch.d \
	$(SRC_DIR)/surface/coonspatch.d \
	$(SRC_DIR)/surface/cubepatch.d \
	$(SRC_DIR)/surface/gmopatch.d \
	$(SRC_DIR)/surface/meshpatch.d \
	$(SRC_DIR)/surface/nozzleexpansionpatch.d \
	$(SRC_DIR)/surface/nurbssurface.d \
	$(SRC_DIR)/surface/package.d \
	$(SRC_DIR)/surface/parametricsurface.d \
	$(SRC_DIR)/surface/ruledsurface.d \
	$(SRC_DIR)/surface/spherepatch.d \
	$(SRC_DIR)/surface/subrangedsurface.d \
	$(SRC_DIR)/surface/sweptpathpatch.d \
	\
	$(SRC_DIR)/volume/meshvolume.d \
	$(SRC_DIR)/volume/nurbsvolume.d \
	$(SRC_DIR)/volume/package.d \
	$(SRC_DIR)/volume/parametricvolume.d \
	$(SRC_DIR)/volume/slabvolume.d \
	$(SRC_DIR)/volume/subrangedvolume.d \
	$(SRC_DIR)/volume/sweptsurfacevolume.d \
	$(SRC_DIR)/volume/tfivolume.d \
	$(SRC_DIR)/volume/twosurfacevolume.d \
	$(SRC_DIR)/volume/wedgevolume.d \
	\
	$(SRC_DIR)/geometry_exception.d \
	\
	$(SRC_DIR)/package.d \

# Look for the static library libplot
LIBPLOT := $(strip $(wildcard /usr/lib/libplot.a) \
                   $(wildcard $(LIBRARY_PATH)/libplot.a))
LIBPLOT_VERSION_STR :=
ifeq ($(findstring libplot,$(LIBPLOT)), libplot)
    $(warning Found libplot:$(LIBPLOT).)
    LIBPLOT_VERSION_STR := with_libplot
    GEOM_D_FILES := $(GEOM_D_FILES) $(SRC_DIR)/misc/libplot.d
endif


GEOM_LUAWRAP_FILES := \
	$(SRC_DIR)/luawrap/package.d \
	$(SRC_DIR)/luawrap/luaunifunction.d \
	$(SRC_DIR)/luawrap/luageom.d \
	$(SRC_DIR)/luawrap/luanomenclature.d \
	$(SRC_DIR)/luawrap/luagpath.d \
	$(SRC_DIR)/luawrap/luagpath_utils.d \
	$(SRC_DIR)/luawrap/luasurface.d \
	$(SRC_DIR)/luawrap/luavolume.d \
	$(SRC_DIR)/luawrap/luagrid.d \
	$(SRC_DIR)/luawrap/luasgrid.d \
	$(SRC_DIR)/luawrap/luausgrid.d \
	$(SRC_DIR)/luawrap/luasketch.d

GEOM_FILES := $(GEOM_D_FILES) $(GEOM_LUAWRAP_FILES)

GEOM_LUA_FILES := $(SRC_DIR)/foam-mesh.lua


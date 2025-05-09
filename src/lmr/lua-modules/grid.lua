-- grid.lua
--
-- Each RegisteredGrid object, in the Lua domain, will hold a reference to a Dlang
-- StructuredGrid or UnstructuredGrid, plus some metadata needed to assign
-- boundary conditions (and the like) at a later point in time.
--
-- PJ, 2021-10-05 pulled out of prep-grids.lua
-- RJG, 2024-03-03 add handling of solid domains
--

local RegisteredGrid = {
   myType = "RegisteredGrid"
}

function RegisteredGrid:new(o)
   -- Input:
   -- A single table with named items.
   -- grid: a StructuredGrid or UnstructuredGrid object that has been generated
   --    or imported.
   -- tag: a string to identify the grid later in the user's script
   -- fieldType: a string labelling the intended domain as 'fluid' or 'solid'
   -- active: a boolean to indicate if the domain on the grid is active (default: true)
   -- omegaz: value of the angular velocity about the z-axis for a rotating frame
   --    Will be 0.0 for a non-rotating-frame.
   -- fsTag: a string that will be used to select the initial flow condition from
   --    a dictionary when the FluidBlock is later constructed.
   -- bcTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the FluidBlock is later constructed.
   -- gridArrayId: needs to be supplied only if the grid is part of a larger array.
   -- ssTag: a string that will be used to select the initial solid condition from
   --    a dictionary when the SolidBlock is later constructed.
   -- solidPropsTag: a string that will be used to select the properties model
   --    for the solid from a dictionary when the SolidBlock is later constructed
   -- solidBCTags: a table of strings that will be used to attach boundary conditions
   --    from a dictionary when the SolidBlock is later constructed
   local flag = type(self)=='table' and self.myType=='RegisteredGrid'
   if not flag then
      error("Make sure that you are using RegisteredGrid:new{} and not RegisteredGrid.new{}", 2)
   end
   local flag = type(o)=='table'
   if not flag then
      error("RegisteredGrid constructor expects a single table with named items.", 2)
   end
   flag = checkAllowedNames(o, {"grid", "tag", "fieldType", "active", "omegaz",
                                "fsTag", "bcTags", "gridArrayId",
                                "ssTag", "solidModelTag", "solidBCTags"})
   if not flag then
      error("Invalid name for item supplied to Grid constructor.", 2)
   end
   setmetatable(o, self)
   self.__index = self
   -- Make a record of the new object, for later use.
   -- Note that we want id to start at zero for the D code.
   o.id = #gridsList
   gridsList[#gridsList+1] = o
   o.tag = o.tag or string.format("Grid-%d", o.id)
   if gridsDict[o.tag] then
      error('Have previously defined a Grid with tag "' .. o.tag .. '"', 2)
   end
   gridsDict[o.tag] = o.id
   -- Set the field type
   -- (we expect the call from registerFluidGrid of registerSolidGrid to set this correctly)
   o.fieldType = o.fieldType or "fluid"
   -- Most common is that grids represent active parts of the domain, however there may be reasons
   -- to deactivate certain regions
   if (o.active == nil) then
      o.active = true
   end
   -- Set to -1 if NOT part of a grid-array, otherwise use supplied value
   o.gridArrayId = o.gridArrayId or -1
   -- Initial FlowState tag
   o.fsTag = o.fsTag or ""
   -- Initial SolidState tag
   o.ssTag = o.ssTag or ""
   -- Solid properties tag
   o.solidModelTag = o.solidModelTag or ""
   -- Default non-rotating frame
   o.omegaz = o.omegaz or 0.0
   -- Must have a grid.
   assert(o.grid, "need to supply a grid")
   -- Check the grid information.
   if config.dimensions ~= o.grid:get_dimensions() then
      local msg = string.format("Mismatch in geometric dimensions, config %d grid %d.",
				config.dimensions, o.grid:get_dimensions())
      error(msg)
   end
   o.bcTags = o.bcTags or {}
   o.solidBCTags = o.solidBCTags or {}
   o.type = o.grid:get_type()
   if o.type == "structured_grid" then
      -- Extract some information from the StructuredGrid
      -- Note 0-based indexing for vertices and cells.
      o.nic = o.grid:get_niv() - 1
      o.njc = o.grid:get_njv() - 1
      if config.dimensions == 3 then
         o.nkc = o.grid:get_nkv() - 1
      else
         o.nkc = 1
      end
      o.ncells = o.nic * o.njc * o.nkc
      -- The following table p for the corner locations,
      -- is to be used later for testing for grid connections.
      o.p = {}
      if config.dimensions == 3 then
         o.p[0] = o.grid:get_vtx(0, 0, 0)
         o.p[1] = o.grid:get_vtx(o.nic, 0, 0)
         o.p[2] = o.grid:get_vtx(o.nic, o.njc, 0)
         o.p[3] = o.grid:get_vtx(0, o.njc, 0)
         o.p[4] = o.grid:get_vtx(0, 0, o.nkc)
         o.p[5] = o.grid:get_vtx(o.nic, 0, o.nkc)
         o.p[6] = o.grid:get_vtx(o.nic, o.njc, o.nkc)
         o.p[7] = o.grid:get_vtx(0, o.njc, o.nkc)
      else
         o.p[0] = o.grid:get_vtx(0, 0)
         o.p[1] = o.grid:get_vtx(o.nic, 0)
         o.p[2] = o.grid:get_vtx(o.nic, o.njc)
         o.p[3] = o.grid:get_vtx(0, o.njc)
      end
      -- Attach default boundary conditions for those not specified.
      for _,face in ipairs(faceList(config.dimensions)) do
         o.bcTags[face] = o.bcTags[face] or o.grid:get_tag(face)
      end
   end
   if o.type == "unstructured_grid" then
      -- Extract some information from the UnstructuredGrid
      o.ncells = o.grid:get_ncells()
      o.nvertices = o.grid:get_nvertices()
      o.nfaces = o.grid:get_nfaces()
      o.nboundaries = o.grid:get_nboundaries()
      -- Attach boundary conditions from list or from the dictionary of conditions.
      for i = 0, o.nboundaries-1 do
         o.bcTags[i] = o.bcTags[i] or o.grid:get_boundaryset_tag(i)
      end
      -- [TODO] RJG, 2024-03-03
      -- Presently, we don't handle unstructured solid domains.
      -- When we do, we'll need to set that information on boundaries here.
   end
   return o
end -- RegisteredGrid:new

function RegisteredGrid:tojson()
   str = '{\n'
   str = str .. string.format('  "tag": "%s",\n', self.tag)
   str = str .. string.format('  "fsTag": "%s",\n', self.fsTag)
   str = str .. string.format('  "type": "%s",\n', self.type)
   str = str .. string.format('  "fieldType": "%s",\n', self.fieldType)
   str = str .. string.format('  "active": %s,\n', self.active)
   str = str .. string.format('  "omegaz": %g,\n', self.omegaz)
   if self.type == "structured_grid" then
      str = str .. string.format('  "dimensions": %d,\n', self.grid:get_dimensions())
      str = str .. string.format('  "niv": %d,\n', self.grid:get_niv())
      str = str .. string.format('  "njv": %d,\n', self.grid:get_njv())
      str = str .. string.format('  "nkv": %d,\n', self.grid:get_nkv())
      str = str .. string.format('  "nic": %d,\n', self.nic)
      str = str .. string.format('  "njc": %d,\n', self.njc)
      str = str .. string.format('  "nkc": %d,\n', self.nkc)
      local fmt = '  "p%d": {"x":%.18e, "y":%.18e, "z":%.18e},\n'
      for i=0, 3 do
         str = str .. string.format(fmt, i, self.p[i].x, self.p[i].y, self.p[i].z)
      end
      if config.dimensions == 3 then
         for i=4, 7 do
            str = str .. string.format(fmt, i, self.p[i].x, self.p[i].y, self.p[i].z)
         end
      end
   else -- unstructured-grid
      str = str .. string.format('  "dimensions": %d,\n', self.grid:get_dimensions())
      str = str .. string.format('  "nvertices": %d,\n', self.grid:get_nvertices())
      str = str .. string.format('  "ncells": %d,\n', self.grid:get_ncells())
      str = str .. string.format('  "nfaces": %d,\n', self.grid:get_nfaces())
      str = str .. string.format('  "nboundaries": %d,\n', self.grid:get_nboundaries())
   end
   str = str .. '  "bcTags": {\n'
   if self.type == "structured_grid" then
      -- Expect named boundaries
      for k, v in pairs(self.bcTags) do
         str = str .. string.format('    "%s": "%s",\n', k, v)
      end
   else -- unstructured_grid
      -- Expect numbered boundary sets.
      -- Note that Dlang numbering starts at zero.
      for j=0, self.grid:get_nboundaries()-1 do
         str = str .. string.format('    "%d": "%s",\n', j, self.bcTags[j])
      end
   end
   str = str .. '    "dummy_entry_without_trailing_comma": "xxxx"\n'
   str = str .. '  },\n'
   str = str .. string.format('  "ssTag": "%s",\n', self.ssTag)
   str = str .. string.format('  "solidModelTag": "%s",\n', self.solidModelTag)
   str = str .. '  "solidBCTags": {\n'
   -- Only handle structured case presently
   for k, v in pairs(self.solidBCTags) do
      str = str .. string.format('    "%s": "%s",\n', k, v)
   end
   str = str .. '    "dummy_entry_without_trailing_comma": "xxxx"\n'
   str = str .. '  },\n'
   str = str .. string.format('  "gridArrayId": %d\n', self.gridArrayId) -- last item, no comma
   str = str .. '}\n'
   return str
end -- end Grid:tojson()

-------------------------------------------------------------------------
--
-- Structured grids may be connected full face to full face.
--
-- needs storage: connectionList = {}

local function connectGrids(idA, faceA, idB, faceB, orientation,
                            reorient_vector_quantities, RmatrixA, RmatrixB)
   -- in 2D, there is only one orientation for connecting a pair of faces.
   -- Since the user will probably not think of providing it, let's default to 0.
   orientation = orientation or 0
   if reorient_vector_quantities == nil then
      reorient_vector_quantities = false
   end
   RmatrixA = RmatrixA or {1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0}
   RmatrixB = RmatrixB or {1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0}
   if false then -- debug
      local msg = string.format('connectGrids(idA=%d, faceA="%s", idB=%d, faceB="%s", orientation=%d,',
                                idA, faceA, idB, faceB, orientation)
      msg = msg .. string.format(' reorient_vector_quantities=%s, RmatrixA=%s, RmatrixB=%s)',
                                 reorient_vector_quantities, RmatrixA, RmatrixB)
      print(msg)
   end
   local gridA = gridsList[idA+1] -- Note that the id values start at zero.
   local gridB = gridsList[idB+1]
   if gridA.grid:get_type() ~= "structured_grid" or gridB.grid:get_type() ~= "structured_grid" then
      error("connectGrids() Works only for structured grids.", 2)
   end
   connectionList[#connectionList+1] = {
      idA=idA, faceA=faceA, idB=idB, faceB=faceB, orientation=orientation,
      reorient_vector_quantities=reorient_vector_quantities, RmatrixA=RmatrixA, RmatrixB=RmatrixB
   }
end

local json = require 'json'

local function connectionAsJSON(c)
   str = '{'
   str = str .. string.format('"idA": %d, "faceA": "%s", "idB": %d, "faceB": "%s", "orientation": %d,',
                              c.idA, c.faceA, c.idB, c.faceB, c.orientation)
   str = str .. string.format(' "reorient_vector_quantities": %s,', tostring(c.reorient_vector_quantities))
   str = str .. string.format(' "RmatrixA": %s,', json.stringify(c.RmatrixA))
   str = str .. string.format(' "RmatrixB": %s', json.stringify(c.RmatrixB))
   str = str .. '}'
   return str
end


local function identifyGridConnections(includeList, excludeList, tolerance)
   -- Identify grid connections by trying to match corner points.
   -- Parameters (all optional):
   -- includeList: the list of structured grid objects to be included in the search.
   --    If nil, the whole collection is searched.
   -- excludeList: list of pairs of structured grid objects that should not be
   --    included in the search for connections.
   -- tolerance: spatial tolerance for the colocation of vertices
   --
   local myGridList = {}
   if includeList then
      -- The caller has provided a list of grids to bound the search.
      for _,A in ipairs(includeList) do
         if A.grid:get_type() == "structured_grid" then
            myGridList[#myGridList+1] = A
         end
      end
   else
      -- The caller has not provided a list; use the global grids list.
      for _,A in ipairs(gridsList) do
         if A.grid:get_type() == "structured_grid" then
            myGridList[#myGridList+1] = A
         end
      end
   end
   excludeList = excludeList or {}
   -- Put unstructured grid objects into the exclude list because they don't
   -- have a simple topology that can always be matched to a structured grid.
   for _,A in ipairs(myGridList) do
      if A.grid:get_type() == "unstructured_grid" then excludeList[#excludeList+1] = A end
   end
   if false then
      print('DEBUG myGridList=[')
      for _,A in ipairs(myGridList) do print('  ', A.id, ',') end
      print(']')
   end
   if #myGridList == 0 then
      -- There are no structured grids that are to be (potentially) connected.
      return
   end
   tolerance = tolerance or 1.0e-6
   --
   for _,A in ipairs(myGridList) do
      for _,B in ipairs(myGridList) do
         if (A ~= B) and (not isPairInList({A, B}, excludeList)) then
            -- print("Proceed with test for coincident vertices.") -- DEBUG
            local connectionCount = 0
            if config.dimensions == 2 then
	          -- print("2D test A.id=", A.id, " B.id=", B.id) -- DEBUG
               for vtxPairs,connection in pairs(connections2D) do
                  if false then -- debug
                     print("vtxPairs=", tostringVtxPairList(vtxPairs),
                           "connection=", tostringConnection(connection))
                  end
                  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
                     local faceA, faceB, orientation = unpack(connection)
                     connectGrids(A.id, faceA, B.id, faceB, 0)
                     connectionCount = connectionCount + 1
                  end
               end
            else
	             -- print("   3D test")
               for vtxPairs,connection in pairs(connections3D) do
                  if verticesAreCoincident(A, B, vtxPairs, tolerance) then
                     local faceA, faceB, orientation = unpack(connection)
                     connectGrids(A.id, faceA, B.id, faceB, orientation)
                     connectionCount = connectionCount + 1
                  end
	             end
	          end
	          if connectionCount > 0 then
	             -- So we don't double-up on connections.
	             excludeList[#excludeList+1] = {A,B}
	          end
	       end -- if (A ~= B...
      end -- for _,B
   end -- for _,A
end -- identifyGridConnections

return {
   RegisteredGrid = RegisteredGrid,
   connectGrids = connectGrids,
   connectionAsJSON = connectionAsJSON,
   identifyGridConnections = identifyGridConnections
}

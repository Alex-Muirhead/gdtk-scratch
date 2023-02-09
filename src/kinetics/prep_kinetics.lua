#!/usr/bin/env dgd-lua
-- Author: Rowan J. Gollan
-- Remix by: Nick Gibbons
-- Date: 2021-04-08
-- History:
--    2023-02-05 Added handling for V-V rate input
--

local mechanism = require 'mechanism'

userMechs = {}
userCCMechs = {}
mechanisms = {}

local validateMechanism = mechanism.validateMechanism
local addUserMechToTable = mechanism.addUserMechToTable
local mechanismToLuaStr = mechanism.mechanismToLuaStr
local addUserChemMechToTable = mechanism.addUserChemMechToTable

function printHelp()
   print("prep-kinetics --- Prepares an energy exchange kinetics file for Eilmer.")
   print("Usage:")
   print(" > prep-kinetics gmodelFile kinInput output")
   print("")
   print("   gmodelFile  : a gas model file is required as input for context")
   print("   kinInput    : input energy exchange kinetics file in Lua format.")
   print("   output      : output file in format ready for Eilmer4.")
   print("")
   os.exit(1)
end


--%--------------------------------------------
-- Function available to user in input file.
--%--------------------------------------------

function Mechanism(t)
   -- Gather mechanisms, but don't yet validate
   userMechs[#userMechs+1] = t
end

function ChemistryCouplingMechanism(t)
   -- Gather chemistry coupling mechanisms, these don't get validated
   t.type = "C-V"
   userCCMechs[#userCCMechs+1] = t
end

-----------------------------------------------

function buildVerboseLuaFile(fName)
   f = assert(io.open(fName, 'w'))
   f:write("\n")
   f:write("mechanism = {}\n\n")
   for i,mech in ipairs(mechanisms) do
      f:write(mechanismToLuaStr(i, mech))
      f:write("\n")
   end
   f:close()
end


function main()
   local gmodelFile, chemmodelFile, inFname, outFname
   
   if (#arg == 0 or arg[1] == "--help") then
      printHelp()
   end

   if (#arg < 3) then
      print("Not enough arguments or unknown option.")
      print("Exiting program without doing anything.")
      printHelp()
   end

   if (#arg == 3) then
       gmodelFile = arg[1]
       inFname    = arg[2]
       outFname   = arg[3]
   end

   if (#arg == 4) then
       gmodelFile = arg[1]
       chemFile   = arg[2]
       inFname    = arg[3]
       outFname   = arg[4]
   end

   if (#arg > 4) then
      print("Too many arguments.")
      print("Exiting program without doing anything.")
      printHelp()
   end


   outstring = 'Creating kinetics file "' .. outFname .. '" using input file "' .. inFname .. '"...'
   print(outstring)

   -- Execute chemical kinetics file in case we need reaction info
   -- We do this before the gas file because both have a table called
   -- "species", and we want the gas file's version to overwrite the chemFile's one.
   if chemFile ~= nil then
      dofile(chemFile)
   else
      reaction = {} 
   end

   -- Execute gas model file so we can get:
   -- 1. list of species
   -- 2. list of energy modes
   dofile(gmodelFile)

   -- The species table has indices as keys, and species names as values.
   -- Let's augment that with a reverse lookup, names as keys and indices as values.
   -- And we'll use the D-offset for species index (from 0)
   for isp,sp in ipairs(species) do
      species[sp] = isp-1
   end
   -- Do the same for energy_modes
   if energy_modes then
      for imode,mode in ipairs(energy_modes) do
	 energy_modes[mode] = imode-1
      end
   else
      -- For 2-T, we don't require the user to set energy modes explicitly
      -- since we can make that decision. So we'll set it up.
      energy_modes = {}
      for isp,sp in ipairs(species) do
	 if db[sp].type == 'molecule' then
	    energy_modes[sp] = 0
	 end
      end
   end

   -- Load contents from user's file.
   dofile(inFname)

   index = 1
   for i,m in ipairs(userMechs) do
      if validateMechanism(m) then
         index = addUserMechToTable(index, m, mechanisms, species, db)
      else
         print("Error while trying to validate mechanism ", i)
         print(m[1])
         print("Bailing out!")
         os.exit(1)
      end
   end
   for i,m in ipairs(userCCMechs) do
      index = addUserChemMechToTable(index, m, mechanisms, species, db, reaction)
   end
   -- Now write out transformed results
   buildVerboseLuaFile(outFname)
end

main()
   


   
   

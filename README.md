# SteamCracking
Attempt to create rigorous Steam Cracker model

Main files:
  - main.py - Controlls whole model. Run to get simulate Plug-Flow Reactor
  - config.py - Contains initial data for Reactor simulation. Edit to change input data
  - coreobjects.py - Describes main objects in model: 'Species' and 'Stream'
  - chemistry.py - Handles kinetic reactions and plug-flow reactor model. Preforms Plug-Flow Reactor simulation when called from main.py
  - databaes/database_v1.py - Storage of 'Species' and 'Reactions' properies colected in small scal from different sources
  - databases/PRBKV_database.py - Storage of binary interaction parameters for Peng-Robison EOS
  - databases/species_database_ussr.py - Storage of 'Species' properties from soviet source
  - databases/reactions_database_ussr.py - Storage of 'Reactions' properties from soviet source
  - plotter.py - Handles grahps

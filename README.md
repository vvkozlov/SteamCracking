# SteamCracking
Attempt to create rigorous Steam Cracker model

Main files:
  - reactors_main.py - Script that controlls whole model. Run to get simulate Plug-Flow Reactor
  - rctrs_config.py - Script with initial data for Reactor simulation. Edit to change input data
  - rctrs_engine.py - Script with all mathematics required to simulate the Reactor. Preforms Plug-Flow Reactor simulation when called from reactors_main.py
  - database_v2.py - Script containing 'Species' and 'Reaction' models as a sort of database
  - PRBKV_database.py - Script containing binary interaction parameters for Peng-Robison EOS as a sort of database

Other files are legacy from other projects and does not used in calculations

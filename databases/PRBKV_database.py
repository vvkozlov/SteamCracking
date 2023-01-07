import pandas as pd

'''
Binary coefficients for PR EOS database:
- coefficients for ESPR model (Kij = k1 + k2*T + k3/T)
- k2 and k3 are not used yet
'''
'''k1 (PRKBV1) coefficients'''
PRKBV1 = pd.DataFrame()
PRKBV1.loc['Propane', 'Propane'] = 0
PRKBV1.loc['Propene', 'Propene'] = 0
PRKBV1.loc['Ethane', 'Ethane'] = 0
PRKBV1.loc['Ethylene', 'Ethylene'] = 0
PRKBV1.loc['Hydrogen', 'Hydrogen'] = 0
PRKBV1.loc['Propane', 'Propene'] = 7.89980031549931e-003
PRKBV1.loc['Propene', 'Propane'] = 7.89980031549931e-003
PRKBV1.loc['Propane', 'Ethane'] = 1.25795602798462e-003
PRKBV1.loc['Ethane', 'Propane'] = 1.25795602798462e-003
PRKBV1.loc['Propane', 'Hydrogen'] = 0.214200004935265
PRKBV1.loc['Hydrogen', 'Propane'] = 0.214200004935265
PRKBV1.loc['Propane', 'Ethylene'] = 2.67040729522705e-003
PRKBV1.loc['Ethylene', 'Propane'] = 2.67040729522705e-003
PRKBV1.loc['Propene', 'Ethane'] = -1.90000003203750e-003
PRKBV1.loc['Ethane', 'Propene'] = -1.90000003203750e-003
PRKBV1.loc['Propene', 'Hydrogen'] = -0.103600002825260
PRKBV1.loc['Hydrogen', 'Propene'] = -0.103600002825260
PRKBV1.loc['Ethane', 'Ethylene'] = 1.22990002855659e-002
PRKBV1.loc['Ethylene', 'Ethane'] = 1.22990002855659e-002
PRKBV1.loc['Ethane', 'Hydrogen'] = 0.223100006580353
PRKBV1.loc['Hydrogen', 'Ethane'] = 0.223100006580353
PRKBV1.loc['Ethylene', 'Hydrogen'] = 7.40000000223517e-003
PRKBV1.loc['Hydrogen', 'Ethylene'] = 7.40000000223517e-003
PRKBV1.loc['Propene', 'Ethylene'] = 1.59275531768799e-003
PRKBV1.loc['Ethylene', 'Propene'] = 1.59275531768799e-003
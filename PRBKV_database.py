import database_v2 as db
import pandas as pd

'''
Binary coefficients for PR EOS database:
- coefficients for ESPR model (Kij = k1 + k2*T + k3/T)
- k2 and k3 are not used yet
'''
'''k1 (PRKBV1) coefficients'''
PRKBV1 = pd.DataFrame()
PRKBV1.loc[db.PROPANE.name, db.PROPYLENE.name] = 7.89980031549931e-003
PRKBV1.loc[db.PROPYLENE.name, db.PROPANE.name] = 7.89980031549931e-003
PRKBV1.loc[db.PROPANE.name, db.ETHANE.name] = 1.25795602798462e-003
PRKBV1.loc[db.ETHANE.name, db.PROPANE.name] = 1.25795602798462e-003
PRKBV1.loc[db.PROPANE.name, db.H2.name] = 0.214200004935265
PRKBV1.loc[db.H2.name, db.PROPANE.name] = 0.214200004935265
PRKBV1.loc[db.PROPANE.name, db.ETHYLENE.name] = 2.67040729522705e-003
PRKBV1.loc[db.ETHYLENE.name, db.PROPANE.name] = 2.67040729522705e-003
PRKBV1.loc[db.PROPYLENE.name, db.ETHANE.name] = -1.90000003203750e-003
PRKBV1.loc[db.ETHANE.name, db.PROPYLENE.name] = -1.90000003203750e-003
PRKBV1.loc[db.PROPYLENE.name, db.H2.name] = -0.103600002825260
PRKBV1.loc[db.H2.name, db.PROPYLENE.name] = -0.103600002825260
PRKBV1.loc[db.ETHANE.name, db.ETHYLENE.name] = 1.22990002855659e-002
PRKBV1.loc[db.ETHYLENE.name, db.ETHANE.name] = 1.22990002855659e-002
PRKBV1.loc[db.ETHANE.name, db.H2.name] = 0.223100006580353
PRKBV1.loc[db.H2.name, db.ETHANE.name] = 0.223100006580353
PRKBV1.loc[db.ETHYLENE.name, db.H2.name] = 7.40000000223517e-003
PRKBV1.loc[db.H2.name, db.ETHYLENE.name] = 7.40000000223517e-003
PRKBV1.loc[db.PROPYLENE.name, db.ETHYLENE.name] = 1.59275531768799e-003
PRKBV1.loc[db.ETHYLENE.name, db.PROPYLENE.name] = 1.59275531768799e-003
import pandas as pd

# numbering of equation, params, range
mil101_heatCap = [[472.6,221.7,0.800,-0.390,-4.5245,5.634],[80,333,0.683331]]

ad_db=pd.DataFrame(columns=['DA_isotherm','enthalpy','heatCap','ads_speed'])
ad_db.loc['MIL-101']=[[28.67, 4562, 1.204],[19.8593], mil101_heatCap, 35.575*10**(-4)] 
ad_db.loc['Uio-66']=[[6.312,9268, 1.916],[24.8994], 2210,35.575*10**(-4)] 

ad_db = ad_db.T


""" heat transfer coeffiecient is also considered
Mei Yang, PLoS One. 2016; 11(7): e0159602.
Fig. 8
"""

""" MIL-101,
Dubinin Astakhov equation
Simulated in RASPA at 273.15, 298.15, 323.15 K
m = qo * exp(-(A/E)^n)
qo, E, n
# fitting R2 = 0.98959

heat capacity
 Liu, S., et al. J Therm Anal Calorim 129, 509–514 (2017). https://doi.org/10.1007/s10973-017-6168-9
% Molar mass of MIL-101
framework cell mass = 46466.506639999025 [g/mol(frame)]
framework mol = 68 [mol(formula)/mol(frame)]
molar mass (formula) = 46466.506639999025 / 68 = 683.331 [g/mol(formula)] = 0.683331 [kg/mol]

kinetics
Zhing Zhang, et al., Energy Fuels 2011, 25, 2, 835–842, https://doi.org/10.1021/ef101548g
LDF model, experiment at 298-328 K, pressure(0->0.5bar)
"""

""" Uio-66,
Dubinin Astakhov equation
Simulated in RASPA at 273.15, 298.15, 323.15 K
m = qo * exp(-(A/E)^n)
qo, E, n
# fitting R2 =

heat capacity
 Qiang Wang, et al., Jounal of nanomaterial, https://doi.org/10.1155/2019/5154173
 ➡ calculate from Forcite flexible frame geometry optimization (less accurate and maybe too high value)
"""
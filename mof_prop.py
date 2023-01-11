import pandas as pd
from fluidProp import VLEFluid
# numbering of equation, params, range
mil101_heatCap = [[472.6,221.7,0.800,-0.390,-4.5245,5.634],[80,333,0.719392]]

ad_db=pd.DataFrame(columns=['DA_isotherm','enthalpy','heatCap','ads_speed'])
ad_db.loc['MIL-101']=[[28.67, 4562, 1.204],[19.8593], mil101_heatCap, 0.03] 
ad_db.loc['Uio-66']=[[8.243, 7222, 1.666],[24.8994], 658.9, 0.05] 


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
molar mass (formula) = 0.719392 [kg/mol]
http://www.metal-organic-frameworks.eu/pdf/adsorbentien/mil101(cr).pdf  # attention(OH)->(F)

kinetics
Zhing Zhang, et al., Energy Fuels 2011, 25, 2, 835–842, https://doi.org/10.1021/ef101548g
LDF model, experiment at 298-328 K, pressure(0->0.5bar)
35.575*10**(-4)
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

class Uio_66():
    def __init__(self):
        self.name = 'Uio-66'
        self.sorbent_data = ad_db[self.name]
        self.gas = VLEFluid('CO2')
        self.R = 8.314462618 # J/K/mol, gas constant
        self.molar_mass = self.gas.get_molar_mass()  # kg_CO2 / mol
        self.__adsorption_prop()
    
    # freudlich coefficient for the isotherm modeling
    def __adsorption_prop(self):
        self.qo_DA, self.E_DA, self.n_DA = self.sorbent_data.DA_isotherm # mol/kg
        self.heat_coeff, = self.sorbent_data.enthalpy


    # adsorption potential in Polyani adsorption theory
    def ad_pot(self, p, T):
        A = self.R * T * log(self.sat_pressure(T) / p)
        return A

    # equilibrium isotherm with Dubinin Astakhov model
    def eq_loading(self,p,T):
        A = self.ad_pot(p,T)
        loading = self.qo_DA * exp(-(A/self.E_DA)**self.n_DA)
        return loading

    # saturation pressure corresponding to the temperature, over critical pressure is assumed to be constant now
    def sat_pressure(self, T):
        """ above critical point, calculated in param_set.py """
        P = -30874525.840160873 + 444031.88866659196 * T + -2208.9087997569713 * T**2 + 3.821325298211822 * T**3
        return P
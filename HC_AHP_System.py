"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from fluidProp import VLEFluid
from math import sqrt, exp
import log_output

# unit
# energy, J
# mass, kg
# temperature, K
# quantity, mol
# volume, m3
# time, second
# pressure, Pa

"""
Refprop version 10 (not version 9) is necessary to read thermal properties of CO2
"""

class MOF():
    def __init__(self, name):
        self.name = name
        self.gas = VLEFluid('CO2')
        self.R = 8.314462618 # J/K/mol, gas constant
        self.molar_mass = self.gas.get_molar_mass()  # kg_CO2 / mol

        # amount of MIL101 in the tank
        self.MOF_mass = 0.15 # kg

        # adsorption speed coefficient in LDF model
        self.K1_LDF = 0.01  # sec

        # the volume of the tank for the flow of CO2 space
        self.Volume = 0.064 # m3

        # heat capacity of the MOF, which was in the manuscript of Joule paper.
        self.cp_mof = 700.0  # J/kg/K, 
        # ??? Liu, S., et al. J Therm Anal Calorim 129, 509–514 (2017). https://doi.org/10.1007/s10973-017-6168-9

        # input flow of CO2 into the tank
        self.m_out = 6.6e-4 / self.molar_mass # mol/sec = kg/s / (kg/mol)

        # heat transfer coefficient of CO2, 
        self.h_CO2 = 1000.0 # W/m2/K, Mei Yang, PLoS One. 2016; 11(7): e0159602.

        # heat transfer coefficient of water
        self.h_water = 600.0 # W/m2/K

        # surface of the MOF tube
        self.surfaceA = 0.04 # m2

        # HTF mass flow
        self.HTF_flow = 5e-6 # m3/s = 5ml/sec

        # heat capacity of HTF, water
        self.rhoCp_HTF = 997.0 * 4182.0  # J/m3/K = kg/m3 * J/kg/K

        # surface area of the inside of the HTF flow
        self.surfaceA_HTF = self.surfaceA * 0.9

        # heat resistance of the wall of the tube
        # thickness of the tube / heat conductance coefficient of the tube / surface area
        self.R_wall_HTF = 0.005 / 386 / ((self.surfaceA_HTF + self.surfaceA)/2)  

        # set variables in the equation
        self.__set_var()
        # freudlich constant for the isotherm 
        self.__freudlich()
        # set initial condition
        self.__IC()
        # declare dependent variables in the equation
        self.__set_dependant_var()
        

    def __set_var(self):
        # loading to the sorbent
        self.loading = None
        # temperature of the gas
        self.gas_T = None
        # temperature of the MOF
        self.mof_T = None

    
    def __set_dependant_var(self):
        # heat capacity of sorbent incuding the adsorbate
        self.cp_sor = None
        # outflow of CO2 from the tank
        self.m_in = None
        # temperature of the outlet of the HTF
        self.T_HTF_out = None
        # mass of the gas in the tank
        self.m_gas = None

        # calcualte mas of the gas from the temeprature and pressure
        self.m_gas = self.calc_mass_gas()
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()
        log_output.log_any_msg('cp_sor == {}'.format(self.cp_sor))
        # efficiency of the Shell and tube model
        self.epsilonC = self.calc_epsilon_C_heat_exchanger()
        

    # define initial conditions, equilibrium
    def __IC(self):
        # initial pressure (before compressor works)
        self.gas_P_init = 300000.0 # Pa
        # set pressure (constant, after compressor works)
        self.gas_P = 4000000.0 # Pa
        # simulation time (duration)
        self.simulation_time = 1200  # sec

        # inlet temperature of HTF, water
        self.T_HTF_in = 303.15
        # input temperature of CO2 gas
        self.T_in = 303.15  # K

        # initial temperature of gas in the tank
        self.gas_T = self.T_in

        # output temperature of HTF is same as inlet
        self.T_HTF_out = self.T_HTF_in
        # intial MOF temperature is same as gas
        self.mof_T = self.gas_T # K
        # initial loading get equilibrium
        self.loading = self.eq_loading(self.gas_P_init, self.mof_T)

        #initial value of the variables for odeint input
        self.loading_init   = self.loading
        self.gas_T_init     = self.gas_T
        self.mof_T_init     = self.mof_T


    # freudlich coefficient for the isotherm modeling
    def __freudlich(self):
        self.K_fre = 33.8647 # mol/kg
        self.n_fre = 1.4621


    # equilibrium isotherm with freudlich isotherm, p is actual pressure, not partial here
    def eq_loading(self,p,T):
        p_sat = self.sat_pressure(T)
        p_par = p / p_sat
        loading = self.K_fre * p_par**(1/self.n_fre)
        return loading


    " linear driving force model"
    def ads_speed_LDF(self):
        # equilibrium amount
        Qe = self.eq_loading(self.gas_P, self.mof_T)
        dQtdt = self.K1_LDF * (Qe - self.loading)
        return dQtdt # mol/kg/sec

    # saturation pressure corresponding to the temperature, over critical pressure is assumed to be constant now
    def sat_pressure(self, T):
        """ above critical point, calculated in param_set.py """
        P = -30874525.840160873 + 444031.88866659196 * T + -2208.9087997569713 * T**2 + 3.821325298211822 * T**3
        return P


    # calculate the mass of the gas in the tank from the density obtained in REFPROP
    def calc_mass_gas(self):
        v_gas = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T).vol #m3/kg
        m_gas = self.Volume / v_gas
        return m_gas * self.molar_mass # mol = m3/kg * kg/mol

    # calculate the heat capacity of sorbent as a sum of MOF and adsorbate
    def calc_heat_Cap_sor(self):
        #cp_liq = self.gas.calc_VLE_liquid_T(self.mof_T).cp * self.molar_mass
        cp_liq = 8.6979 * self.molar_mass * 1000
        cp_sor = self.cp_mof * self.MOF_mass + self.loading * self.MOF_mass * cp_liq
        return cp_sor # J/K

    # calculate the heat from adsorption amount
    def dHadsdm(self):
        dH = 22 # kJ/mol
        return dH * 1000 # J/mol


    # the efficiency of the heat exchanger, which changes by increasing heat capacity of sorbent, loadings
    # The number of transfer unit of heat exchanger between sorbent and HTF
    # https://jp.mathworks.com/help/physmod/hydro/ref/entuheattransfer.html#:~:text=NTU%20is%20the%20number%20of,or%20finned%2C%20heat%20transfer%20surfaces.
    def calc_epsilon_C_heat_exchanger(self):
        R_overall = 1/(self.surfaceA * self.h_CO2) + self.R_wall_HTF + 1/(self.surfaceA_HTF * self.h_water)
        NTU = 1/(self.rhoCp_HTF*self.HTF_flow)/R_overall

        epsilon = 1 - exp(-NTU)
        log_output.log_any_msg(epsilon)
        return epsilon * self.rhoCp_HTF*self.HTF_flow
    

def test_cp_sor(mof_):
    mof_.calc_heat_Cap_sor()

if __name__=='__main__':
    mof = MOF('MIL-101')
    test_cp_sor(mof)


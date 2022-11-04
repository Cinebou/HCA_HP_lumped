"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from fluidProp import VLEFluid
from math import log, sqrt, exp
from mof_prop import ad_db
# unit
# energy, J
# mass, kg
# temperature, K
# CO2 quantity, mol
# volume, m3
# time, second
# pressure, Pa

"""
Refprop version 10 (not version 9) is necessary to read thermal properties of CO2
"""

class MOF():
    def __init__(self, name):
        self.name = name
        self.sorbent_data = ad_db[name]
        self.gas = VLEFluid('CO2')
        self.R = 8.314462618 # J/K/mol, gas constant
        self.molar_mass = self.gas.get_molar_mass()  # kg_CO2 / mol

        # amount of MIL101 in the tank
        self.MOF_mass = 0.15 # kg

        # adsorption speed coefficient in LDF model
        self.K1_LDF = self.sorbent_data.ads_speed  # sec

        # the volume of the tank for the flow of CO2 space
        self.Volume = 0.064 # m3
    
        # output flow of CO2 into the tank
        self.m_out = 6.6e-4 / self.molar_mass # mol/sec = kg/s / (kg/mol)

        # heat transfer coefficient of CO2, 
        self.h_CO2 = 1000 # W/m2/K, Mei Yang, PLoS One. 2016; 11(7): e0159602.

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
        # constant for this adsorbent
        self.__adsorption_prop()
        # set initial condition
        self.__IC()
        # declare dependent variables in the equation
        self.set_dependant_var_IC()
        

    def __set_var(self):
        # loading to the sorbent
        self.loading = None
        # temperature of the gas
        self.gas_T = None
        # temperature of the MOF
        self.mof_T = None


    # define initial conditions, equilibrium
    def __IC(self):
        # simulation time (duration)
        self.simulation_time = 500  # sec
        # initial temperature of gas in the tank
        self.gas_T = 293.15
        # initial pressure (before compressor works)
        self.gas_P_init = 800000. # Pa
        # set pressure (constant, after compressor works)
        self.gas_P = 2020833.3333333335# Pa
        # time diff
        self.dt = 0.1 # sec, 1/dt should be integral

    
    def set_dependant_var_IC(self):
        # temperature of the gas before compression
        self.T_lower = self.gas_T
        # inlet temperature of HTF, water, is same as initial temperature of the system
        self.T_HTF_in = self.gas_T
        # output temperature of HTF is same as inlet
        self.T_HTF_out = self.T_HTF_in
        # intial MOF temperature is same as gas
        self.mof_T = self.gas_T # K

        # initial loading get equilibrium
        self.loading = self.eq_loading(self.gas_P_init, self.mof_T)
    
        # calcualte mas of the gas from the temeprature and pressure
        self.m_gas = self.calc_mass_gas()
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()
        # efficiency of the Shell and tube model, it is a fixed value
        self.epsilonC = self.calc_epsilon_C_heat_exchanger()
        # input temperature of CO2 gas
        self.T_in, h_after_comp = self.calc_gas_in(self.T_lower, self.gas_P_init, self.gas_P)  # K
        # compressor work
        h_before_comp = self.gas.calc_fluidProp_pT(self.gas_P_init, self.gas_T).h
        self.h_comp = (h_after_comp - h_before_comp)* self.molar_mass # J/mol = J/kg * kg/mol
        
        #initial value of the variables for odeint input
        self.loading_init   = self.loading
        self.gas_T_init     = self.gas_T
        self.mof_T_init     = self.mof_T

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
        cp_liq = self.gas.calc_fluidProp_pT(self.sat_pressure(self.mof_T),self.mof_T).cp * self.molar_mass
        cp_mof = self.heat_capacity_mof(self.mof_T)
        cp_sor = cp_mof * self.MOF_mass + self.loading * self.MOF_mass * cp_liq
        return cp_sor # J/K

    # heat capacity fitting equation from literature
    def heat_capacity_mof(self, temp):
        if self.name == 'MIL-101':
            c1, c2, c3, c4, c5, c6 = self.sorbent_data.heatCap[0]
            a1, a2, a3 =  self.sorbent_data.heatCap[1]
            t = (temp - a2)/a1
            cp = (c1 + c2*t + c3*t**2 + c4*t**3 + c5*t**4 + c6*t**5 ) #/ a3
        if self.name == 'Uio-66':
            cp = self.sorbent_data.heatCap
        return cp


    # calculate the heat from adsorption amount
    def dHadsdm(self):
        dH = self.heat_coeff# kJ/mol
        return dH * 1000 # J/mol


    # the efficiency of the heat exchanger, which changes by increasing heat capacity of sorbent, loadings
    # The number of transfer unit of heat exchanger between sorbent and HTF
    # https://jp.mathworks.com/help/physmod/hydro/ref/entuheattransfer.html#:~:text=NTU%20is%20the%20number%20of,or%20finned%2C%20heat%20transfer%20surfaces.
    def calc_epsilon_C_heat_exchanger(self):
        R_overall = 1/(self.surfaceA * self.h_CO2) + self.R_wall_HTF + 1/(self.surfaceA_HTF * self.h_water)
        NTU = 1/(self.rhoCp_HTF*self.HTF_flow)/R_overall
        epsilon = 1 - exp(-NTU)
        return epsilon * self.rhoCp_HTF*self.HTF_flow

    
    # decide the temperature of the gas just after the addiabatic compression
    def calc_gas_in(self, Tl, pl, ph):
        # entropy of the gas before the compression
        sl = self.gas.calc_fluidProp_pT(pl, Tl).s
        # addiabatic compression means entropy doesn't change
        gas_after_comp = self.gas.calc_fluidProp_ps(ph, sl)
        Th_in = gas_after_comp.T
        h_in = gas_after_comp.h
        return Th_in, h_in


if __name__=='__main__':
    mof = MOF('MIL-101')


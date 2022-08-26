"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from fluidProp import VLEFluid
from math import sqrt, exp

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
        self.molar_mass = self.gas.get_molar_mass() / 1000  # kg_CO2 / mol, 4.40098e-05

        # amount of MIL101 in the tank
        self.MOF_mass = 0.15 # kg

        # adsorption speed coefficient in pseudo second order model
        # this model is for chemisorption basically, I don't know why the author uses this model, it should be refined because this is the most important part in this system
        self.K2 = 4.26    # kg/mol/sec

        # adsorption speed coefficient in LDF model
        self.K1_LDF = 0.01  # sec

        # the volume of the tank for the flow of CO2 space
        self.Volume = 0.005 # m3

        # heat capacity, cp_gas should be refined by refprop
        #self.cp_gas = 36.84   # J/mol/K

        # heat capacity of the MOF, which was in the manuscript of Joule paper.
        self.cp_mof = 700.0  # J/kg/K, 
        # ??? Liu, S., et al. J Therm Anal Calorim 129, 509–514 (2017). https://doi.org/10.1007/s10973-017-6168-9

        # input flow of CO2 into the tank
        self.m_out = 0.01 # mol/sec

        # heat transfer coefficient of CO2, 
        self.h_CO2 = 1000.0 # W/m2/K, Mei Yang, PLoS One. 2016; 11(7): e0159602.

        # heat transfer coefficient of water
        self.h_water = 5000.0 # W/m2/K

        # surface of the MOF tube
        self.surfaceA = 0.2 # m2

        # HTF mass flow
        self.HTF_flow = 5e-6 # m3/s = 5ml/sec

        # heat capacity of HTF, water
        self.rhoCp_HTF = 997.0 * 4182.0  # J/m3/K = kg/m3 * J/kg/K

        # surface area of the inside of the HTF flow
        self.surfaceA_HTF = self.surfaceA * 0.9

        # heat resistance of the wall of the tube
        # thickness of the tube / heat conductance coefficient of the tube / surface area
        self.R_wall_HTF = 0.001 / 0.94 / ((self.surfaceA_HTF + self.surfaceA)/2)  

        # set variables in the equation
        self.__set_var()
        # declare dependent variables in the equation
        self.__set_dependant_var()
        # freudlich constant for the isotherm 
        self.__freudlich()
        # set initial condition
        self.__IC()
        

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
        

    # define initial conditions, equilibrium
    def __IC(self):
        # initial pressure
        self.gas_P_init = 600000.0
        # initial temperature
        self.gas_T = 283.15
        # set pressure (constant)
        self.gas_P = 2500000.0
        # simulation time
        self.simulation_time = 50  # sec

        # inlet temperature of HTF
        self.T_HTF_in = 283.15
        # input temperature of CO2 gas
        self.T_in = 283.15  # K

        # output temperature of HTF is same as inlet
        self.T_HTF_out = self.T_HTF_in
        # intial MOF temperature is same as gas
        self.mof_T = self.gas_T # K
        # initial loading get equilibrium
        self.loading = self.eq_loading(self.gas_P_init, self.mof_T)
        # initial amount of gas in the tank
        self.m_gas = self.calc_mass_gas()

        self.cp_sor = self.calc_heat_Cap_sor()

        #initial value of the variables for odeint input
        self.m_gas_init     = self.m_gas
        self.loading_init   = self.loading
        self.gas_T_init     = self.gas_T
        self.mof_T_init     = self.mof_T
        self.T_HTF_out_init = self.T_HTF_out
        
        # previous value
        # mass of CO2 gas in the tank at previous cycle
        self.m_gas_previous = self.calc_mass_gas()


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
        Qe = self.eq_loading(self.gas_P, self.gas_T)
        dQtdt = self.K1_LDF * (Qe - self.loading)
        return dQtdt # mol/kg/sec


    " pseudo second order model"
    # adsorption speed defined by pseudo-second-order model
    # (Hui Su et al., RCS Adv., 2020, 10, 2198)
    # Qt = k2*Qe^2/(1 + k2*Qe*t) * t
    # t = Qt / (k2*Qe*(Qe - Qt))
    # dQt/dt = k2*Qe^2 / (1 + k2*Qe*t)^2 = k2 * (Qe - Qt)^2
    def ads_speed_second(self):
        # equilibrium amount
        Qe = self.eq_loading(self.gas_P, self.gas_T)
        # loading speed
        dQtdt = self.K2 * (Qe - self.loading)**2
        # return the mass adsorbed in the tank
        return dQtdt  # mol/kg/sec


    # saturation pressure corresponding to the temperature, over critical pressure is assumed to be constant now
    def sat_pressure(self, T):
        f = self.gas.calc_VLE_T(T)
        p = f.p_v
        if p > 0:
            return p
        else:
            #print('pressure over saturation point')
            return self.gas.calc_VLE_T(303).p_v


    # calculate the mass of the gas in the tank from the density obtained in REFPROP
    def calc_mass_gas(self):
        v_gas = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T).vol #m3/kg
        m_gas = self.Volume / v_gas
        return m_gas # kg

    # calculate the heat capacity of sorbent as a sum of MOF and adsorbate
    def calc_heat_Cap_sor(self):
        cp_liq = self.gas.calc_VLE_liquid_T(self.gas_T).cp * self.molar_mass
        cp_sor = self.cp_mof * self.MOF_mass + self.loading * self.MOF_mass * cp_liq
        return cp_sor # J/K

    # calculate the heat from adsorption amount
    def dHadsdm(self):
        dH = 19
        return 19 * 1000 # J/mol


    # the efficiency of the heat exchanger, which changes by increasing heat capacity of sorbent, loadings
    def calc_epsilon_C_heat_exchanger(self):
        # The number of transfer unit of heat exchanger between sorbent and HTF
        # https://jp.mathworks.com/help/physmod/hydro/ref/entuheattransfer.html#:~:text=NTU%20is%20the%20number%20of,or%20finned%2C%20heat%20transfer%20surfaces.
        C_min_HTF = min(self.cp_sor, self.rhoCp_HTF*self.HTF_flow)
        C_rel = min(self.cp_sor, self.rhoCp_HTF*self.HTF_flow) / max(self.cp_sor, self.rhoCp_HTF*self.HTF_flow)
        R_overall = 1/(self.surfaceA * self.h_CO2) + self.R_wall_HTF + 1/(self.surfaceA_HTF * self.h_water)
        NTU = 1/C_min_HTF/R_overall

        epsilon = 2/(1 + C_rel + sqrt(1+C_rel**2)*(1+exp(-NTU*sqrt(1+C_rel**2)))/(1-exp(-NTU*sqrt(1+C_rel**2))))
        return epsilon * C_min_HTF
    


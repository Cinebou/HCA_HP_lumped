"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from fluidProp import VLEFluid
from math import log, exp, pi
from mof_prop import ad_db
from exp_data.exp_sumarize import ResData
import matplotlib.pyplot as plt


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
    def __init__(self, name, res_file, start_time):
        self.name = name
        self.res_file = res_file
        self.start_time = start_time
        self.sorbent_data = ad_db[name]

        self.gas = VLEFluid('CO2')
        self.R = 8.314462618 # J/K/mol, gas constant
        self.cp_water = 4182.0 #  J/kg/K
        self.molar_mass = self.gas.get_molar_mass()  # kg_CO2 / mol

        # amount of MIL101 in the tank
        self.MOF_mass = 0.1349 # kg

        # adsorption speed coefficient in LDF model
        self.K1_LDF = self.sorbent_data.ads_speed  # sec

        # the volume of the tank for the flow of CO2 space
        self.Volume = 0.00368 # m3, calculated from Tanaka san sheet

        # the area of inner surface area of the vessel fo the tank
        # for the calculation of heat loss
        self.A_vessel = 0.11463 # m2
        self.k_water = 1.0
    
        # heat transfer coefficient of CO2, 
        self.h_CO2 = 20 # W/m2/K, Mei Yang, PLoS One. 2016; 11(7): e0159602.

        # surface of the MOF tube
        self.surfaceA = 0.33909 # m2, = num_fin(178) * 0.01905 m * 0.100 m

        # heat capacity of HTF, water
        self.rhoCp_HTF = 997.0 * self.cp_water  # J/m3/K = kg/m3 * J/kg/K

        # surface area of the inside of the HTF flow
        self.surfaceA_HTF = 3.66 * 2*pi * 0.004 # m * 2pi * r
        # mass of water in the pipe
        self.m_water = 3.66 * (0.004**2*pi) * 997.0  # m * m2 * kg/m3 = kg

        # heat resistance of the wall of the tube
        # thickness of the tube / heat conductance coefficient of the tube / surface area
        self.R_wall_HTF = 0.0003 / 386 / ((self.surfaceA_HTF + self.surfaceA)/2)  

        # simulation time (duration)
        self.simulation_time =180 # sec
        # time diff
        self.dt = 0.1 # sec, 1/dt should be integral


        # constant for this adsorbent
        self.__adsorption_prop()
        # set initial condition
        # get condition from experiment 
        if res_file != '':
            self.get_condition_from_exp()
        else:
            self.__IC()
        # declare dependent variables in the equation
        self.set_dependant_var_IC()


    # define initial conditions, equilibrium
    def __IC(self):
        # initial temperature of gas in the tank
        self.gas_T = 21.691 + 273.15  # K
        # initial pressure (before compressor works)
        self.gas_P_init = 0.30908571428571435 * 1e06 # Pa
        # set pressure (constant, after compressor works)
        self.gas_P =  1.8127941176470586 * 1e06 # Pa
        # inlet temperature of HTF, water, is same as initial temperature of the system
        self.T_HTF_in = 15.382857142857139 + 273.15   #self.gas_T
        # initial temperature of HTF
        self.T_HTF = self.T_HTF_in
        # intial MOF temperature 
        self.mof_T = 15.696666666666667 + 273.15 # K
        # experiment setting
        self.T_in = 20.21555555555556 + 273.15 # K
        # temperature of the gas before compression
        self.T_lower = 283.15
        # output flow of CO2 into the tank
        self.m_out = 11.175/60/60 / self.molar_mass # mol/sec = kg/h/3600 / (kg/mol)
        # HTF mass flow
        self.HTF_flow = 0.6278888888888888  /1000/60 # m3/s = l/min / 60
        # initial temp of vessel
        self.vessel_T = 25 + 273.15

    
    def set_dependant_var_IC(self):
        # initial loading get equilibrium
        self.loading = self.eq_loading(self.gas_P_init, self.mof_T)

        # heat transfer coefficient of water
        self.h_water = self.HTrans_water()
        # calcualte mas of the gas from the temeprature and pressure
        self.m_gas = self.calc_mass_gas()
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()
        # efficiency of the Shell and tube model, it is a fixed value
        #self.epsilonC = self.calc_epsilon_C_heat_exchanger()

        """
        # input temperature of CO2 gas
        self.T_in, h_after_comp = self.calc_gas_in(self.T_lower, self.gas_P_init, self.gas_P)  # K
        # compressor work
        h_before_comp = self.gas.calc_fluidProp_pT(self.gas_P_init, self.gas_T).h
        self.h_comp = (h_after_comp - h_before_comp)* self.molar_mass # J/mol = J/kg * kg/mol
        """
        #initial value of the variables for odeint input
        self.loading_init   = self.loading
        self.gas_T_init     = self.gas_T
        self.mof_T_init     = self.mof_T
        self.T_HTF_init     = self.T_HTF

    
    # get the data from literature
    def get_condition_from_exp(self):
        #res_file = './exp_data/res/GL840_01_No2_2022-12-01_17-18-16.csv'; start_sec=3100; exp_span = 450
        #res_file = './exp_data/res/GL840_01_No2_2022-12-06_15-12-36.csv';start_sec=750; exp_span = 450
        time = int(self.simulation_time * 2.2)
        self.res_exp = ResData(self.res_file, self.start_time, exp_span=time, dir_name='')

        # input CO2 temperature
        self.exp_T_in = self.res_exp.d_loc_time("CO2 inlet", self.res_exp.start_t, self.res_exp.finish_t)
        # initial temperature of gas in the tank
        self.gas_T = self.res_exp.T_init_gas + 273.15  # K
        # initial pressure (before compressor works)
        self.gas_P_init = self.res_exp.p_init * 1e06 # Pa
        # initial temperature of water
        self.T_HTF = self.res_exp.T_inlet_HTF + 273.15 # K
        # set pressure (constant, after compressor works)
        self.gas_P =  self.res_exp.p_high * 1e06 # Pa
        # inlet temperature of HTF, water, is same as initial temperature of the system
        self.T_HTF_in = self.res_exp.T_inlet_HTF + 273.15   #self.gas_T
        # intial MOF temperature 
        self.mof_T = self.res_exp.T_init_mof + 273.15 # K
        # experiment setting
        self.T_in = self.mof_T #38.21555555555556 + 273.15 # K
        # output flow of CO2 into the tank
        self.m_out = self.res_exp.CO2_flow /60/60 / self.molar_mass # mol/sec = kg/h/3600 / (kg/mol)
        
        # modify experiment mistake
        self.m_out = self.m_out / 2

        # HTF mass flow
        self.HTF_flow = self.res_exp.water_flow /1000/60 # m3/s = l/min / 60
        # initial temp of vessel
        self.vessel_T = self.res_exp.T_init_vessel + 273.15


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
        v_gas = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T).vol * self.molar_mass #m3/mol =  (m3/kg) * (kg/mol)
        m_gas = self.Volume / v_gas
        return m_gas  # mol 


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
            cp = (c1 + c2*t + c3*t**2 + c4*t**3 + c5*t**4 + c6*t**5 ) / a3
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
        R_overall = self.R_wall_HTF + 1/(self.surfaceA_HTF * self.h_water)
        NTU = 1/(self.rhoCp_HTF*self.HTF_flow)/R_overall
        self.epsilon = 1 - exp(-NTU)
        """
        print(
              "heat cond. copper: ",self.R_wall_HTF/R_overall,\
              "heat trans. HTF  : ", 1/(self.surfaceA_HTF * self.h_water)/R_overall
            )
        """
        return self.epsilon * self.rhoCp_HTF*self.HTF_flow

    
    # decide the temperature of the gas just after the addiabatic compression
    def calc_gas_in(self, Tl, pl, ph):
        # entropy of the gas before the compression
        sl = self.gas.calc_fluidProp_pT(pl, Tl).s
        # addiabatic compression means entropy doesn't change
        gas_after_comp = self.gas.calc_fluidProp_ps(ph, sl)
        Th_in = gas_after_comp.T
        h_in = gas_after_comp.h
        return Th_in, h_in

    # heat transfer coefficeint of water
    def HTrans_water(self):
        d = 0.008 - 0.0003                          # intube diameter
        A = pi * (d/2)**2                           # intube area 
        Re = (self.HTF_flow/A) * d / (1.004*1e-6)   # Reynolds number = flow speed * intube diameter / kinetic viscocity
        Pr = 6.9                                    # water Prandtl number, Microfluidics: Modeling, Mechanics and Mathematics A volume in Micro and Nano Technologies Book • 2017
        St = 0.023 / Pr**(2/3) / Re**(0.2)          # stanton number, 伝熱工学資料　日本機械学会 P221, 液体円管流れの熱伝達計算
        h_water = self.cp_water * (self.HTF_flow*997) / (A) * St # heat transfer = cp m_/A * St
        return h_water * self.k_water


if __name__=='__main__':
    mof = MOF('MIL-101')
    mof.get_condition_from_exp()
    


"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from fluidProp import VLEFluid
from math import log, exp, pi
from exp_data.exp_sumarize import ResData
import matplotlib.pyplot as plt
from DubininAstakov import DA_model

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

class MOF(DA_model):
    def __init__(self, name, res_file, start_time):
        self.name = name
        self.res_file = res_file
        self.start_time = start_time
        super().__init__()

        # discritization number in longtitudnal direction
        self.n_discrete = 5

        self.gas = VLEFluid('CO2')
        self.R = 8.314462618 # J/K/mol, gas constant
        self.cp_water = 4182.0 #  J/kg/K
        self.cp_mof = 658.9 # J/K/kg
        self.cp_HX = 350 # J/K
        self.molar_mass = self.gas.get_molar_mass()  # kg_CO2 / mol

        # amount of MIL101 in the tank
        #self.MOF_mass = 0.1349 # kg
        self.MOF_mass = 0.2608 # kg

        # adsorption speed coefficient in LDF model
        self.K1_LDF = 0.09  # sec

        # the volume of the tank for the flow of CO2 space
        self.Volume = 0.00368 # m3, calculated from Tanaka san sheet

        # the area of inner surface area of the vessel fo the tank
        # for the calculation of heat loss
        self.A_vessel = 0.11463 # m2
        self.k_water = 0.4
    
        # heat transfer coefficient of CO2, 
        self.h_CO2 = 30 # W/m2/K, Mei Yang, PLoS One. 2016; 11(7): e0159602.

        # surface of the MOF tube
        self.surfaceA = 0.33909 # m2, = num_fin(178) * 0.01905 m * 0.100 m

        # heat capacity of HTF, water
        self.rhoCp_HTF = 997.0 * self.cp_water  # J/m3/K = kg/m3 * J/kg/K

        # surface area of the inside of the HTF flow
        self.length_tube = 3.66
        self.d_tube = 0.004 - 0.0003
        self.surfaceA_HTF = self.length_tube * 2*pi * 0.004 # m * 2pi * r
        # mass of water in the pipe
        self.m_water = self.length_tube * (self.d_tube**2*pi) * 997.0  # m * m2 * kg/m3 = kg

        # heat resistance of the wall of the tube
        # thickness of the tube / heat conductance coefficient of the tube / surface area
        self.R_wall_HTF = 0.0003 / 386 / ((self.surfaceA_HTF + self.surfaceA)/2)  

        # simulation time (duration)
        self.simulation_time = 120 # sec
        # time diff
        self.dt = 0.01 # sec, 1/dt should be integral


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
        self.gas_T = 25 + 273.15  # K
        # initial pressure (before compressor works)
        self.gas_P_init = 0.30 * 1e06 # Pa
        # set pressure (constant, after compressor works)
        self.gas_P =  3 * 1e06 # Pa
        # inlet temperature of HTF, water, is same as initial temperature of the system
        self.T_HTF_in = 25 + 273.15   #self.gas_T
        # initial temperature of HTF
        self.T_HTF = self.T_HTF_in
        # intial MOF temperature 
        self.mof_T = 25 + 273.15 # K
        # experiment setting
        self.T_in = 25 + 273.15 # K
        # temperature of the gas before compression
        self.T_lower = 283.15
        # output flow of CO2 into the tank
        self.m_out = 9.7/60/60 / self.molar_mass # mol/sec = kg/h/3600 / (kg/mol)
        # HTF mass flow
        self.HTF_flow = 0.5 /1000/60 # m3/s = l/min / 60
        # initial temp of vessel
        self.vessel_T = 25 + 273.15

    
    def set_dependant_var_IC(self):
        #discritization
        self.dx = self.length_tube / self.n_discrete

        # initial loading get equilibrium
        self.loading = self.eq_loading(self.gas_P_init, self.mof_T)

        # heat transfer coefficient of water
        self.h_water = self.HTrans_water()
        # calcualte mas of the gas from the temeprature and pressure
        self.m_gas = self.calc_mass_gas()
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()

        #initial value of the variables for odeint input
        self.loading_init   = self.loading
        self.gas_T_init     = self.gas_T
        self.mof_T_init     = self.mof_T
        self.T_HTF_init     = [self.T_HTF]*(self.n_discrete+2)

    
    # get the data from literature
    def get_condition_from_exp(self):
        time = int(self.simulation_time*1.2)
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


    " linear driving force model"
    def ads_speed_LDF(self):
        # equilibrium amount
        Qe = self.eq_loading(self.gas_P, self.mof_T)
        dQtdt = self.K1_LDF * (Qe - self.loading)
        return dQtdt # mol/kg/sec


    # calculate the mass of the gas in the tank from the density obtained in REFPROP
    def calc_mass_gas(self):
        v_gas = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T).vol * self.molar_mass #m3/mol =  (m3/kg) * (kg/mol)
        m_gas = self.Volume / v_gas
        return m_gas  # mol 


    # differential adsorption amount with different tempeature
    def dMass_gas_dT(self):
        # m_gas, present
        v_gas = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T).vol * self.molar_mass #m3/mol =  (m3/kg) * (kg/mol)
        m_gas = self.Volume / v_gas
        # m_gas with slightly increased temperature
        dT = 0.1 # K
        v_gas_dv = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T + dT).vol * self.molar_mass #m3/mol =  (m3/kg) * (kg/mol)
        m_gas_dm = self.Volume / v_gas_dv
        return (m_gas_dm - m_gas) / dT


    # calculate the heat capacity of sorbent as a sum of MOF and adsorbate
    def calc_heat_Cap_sor(self):
        cp_ads = self.specific_heat_Cap(self.gas_P, self.mof_T, self.loading) # J/mol/K
        cp_sor = self.cp_mof * self.MOF_mass + self.loading * self.MOF_mass * cp_ads 
        return cp_sor + self.cp_HX # J/K



    # calculate the heat from adsorption amount
    def dHadsdm(self):
        dH = self.current_h_ads(self.loading, self.mof_T)
        return dH * 1000 # J/mol


    # heat transfer coefficeint of water
    def HTrans_water(self):
        A = pi * (self.d_tube)**2                           # intube area 
        Re = (self.HTF_flow/A) * self.d_tube / (1.004*1e-6)   # Reynolds number = flow speed * intube diameter / kinetic viscocity
        Pr = 6.9                                    # water Prandtl number, Microfluidics: Modeling, Mechanics and Mathematics A volume in Micro and Nano Technologies Book • 2017
        St = 0.023 / Pr**(2/3) / Re**(0.2)          # stanton number, 伝熱工学資料　日本機械学会 P221, 液体円管流れの熱伝達計算
        h_water = self.cp_water * (self.HTF_flow*997) / (A) * St # heat transfer = cp m_/A * St
        #print('Re  ', Re)
        #print('Nu   ', h_water*d/0.602)
        #print('v : ',self.HTF_flow/A, ' m/s')
        return h_water * self.k_water


if __name__=='__main__':
    mof = MOF('MIL-101')
    mof.get_condition_from_exp()
    


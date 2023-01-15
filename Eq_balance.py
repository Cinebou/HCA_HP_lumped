from scipy.integrate import ode,cumtrapz, solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from HC_AHP_System import MOF
import log_output 
from os import mkdir, path
from pprint import pprint
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 320
if not path.exists('./Fig'):
    mkdir('./Fig')


def main():
    #mof_type = 'MIL-101'
    mof_type = 'Uio-66'
    res_file = './exp_data/res/2022-12-26/GL840_01_No3_2022-12-26_15-44-03.csv';start_sec=18000
    
    System = balance(mof_type,res_file,start_sec)
    #System = balance(mof_type,'',start_sec)
    totQ = System.solver()

    print('Total hot water of the system is :: ', totQ,' J')
    System.make_gragh(System.plot_data, System.legends)
    #plt.show()


class balance(MOF):
    def __init__(self,name,res_file,start_time):
        super().__init__(name,res_file,start_time)

    # update system with each iteration
    def update_variables(self,t):
        # increase of the loading
        self.dm_loading = self.ads_speed_LDF()
        # increase of the adsorption amount
        self.m_ads = self.dm_loading * self.MOF_mass
        # the mass of CO2 in the tank
        self.m_gas = self.calc_mass_gas()
        # the mass of the outflow gas, record m_out for the calucaltion in energy balance
        #self.m_in = self.m_out + self.m_ads
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()
        # compressor work
        self.compressor_work = 0 
        #get input temperature of CO2 from experiment
        self.T_in = float(self.exp_T_in.iloc[[int(t)]]) + 273.15; assert 270 <= self.T_in <= 400
        
        #boundary condition
        self.T_HTF[0] = self.T_HTF_in

    # mass balance of the adsorbate
    def mass_balance_sor(self,t):
        # write down the log
        #log_output.log_mass_msg(self.m_in, self.m_ads, self.m_out,t)
        return self.dm_loading

    
    # energy balance of the gas in the tank
    def en_balance_gas(self,t):
        # enthalpy of the CO2 flow
        h_in_CO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.T_in).h  * self.molar_mass # J / mol_CO2
        # refprop data of the gas in the tank
        PropCO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T)
        #PropCO2 = self.gas.calc_fluidProp_pT(self.gas_P, 303.15)
        h_out_CO2 = PropCO2.h  * self.molar_mass # J / mol_CO2

        # differential mass of gas in the tank by temperature 
        dm_gas_dT = self.dMass_gas_dT()* self.molar_mass

        # energy flow from the CO2 stream
        en_flow = self.m_out * (h_in_CO2 - h_out_CO2)

        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.mof_T - self.gas_T)

        # expansion heat due to the adsorption
        mRT = self.m_ads * (self.m_ads * h_in_CO2 - self.gas_P * PropCO2.vol * self.molar_mass)

        # heat loss
        Q_loss = self.h_CO2 * self.A_vessel * (self.vessel_T - self.gas_T)

        # temperature development of the gas
        cp_gas = PropCO2.cp * self.molar_mass
        dTgasdt = (en_flow + Htrans + mRT + Q_loss) / (self.m_gas*cp_gas - h_in_CO2*dm_gas_dT)

        # Logger
        #log_output.log_any_msg('{},{},{}, {}'.format(en_flow, Htrans, mRT, Q_loss))
        return dTgasdt
    

    # energy balance of the adsorbent
    def en_balance_sor(self,t):
        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.gas_T - self.mof_T)
        # adsorption heat
        Hads = self.m_ads * self.dHadsdm()
        # heat passed to the HTF
        pass_HTF = self.h_water * self.surfaceA_HTF * (np.mean(self.T_HTF[1:-1]) - self.mof_T)
        
        # temperature development of the sorbent
        dTsordt = (Htrans + Hads + pass_HTF) / self.cp_sor

        # Logger
        #log_output.log_h_sor_msg(dTsordt*self.cp_sor, np.mean(self.T_HTF), self.mof_T, pass_HTF,t)
        #log_output.log_any_msg('{},{}'.format(t,np.mean(self.T_HTF[1:-1])- self.mof_T))
        return dTsordt


    # energy balance of the water
    def en_balance_water(self,t):
        # heat transfer between sorbent and gas
        pass_HTF = self.h_water * self.surfaceA_HTF * (self.mof_T - self.T_HTF[1:-1])/(self.length_tube) # W/m
        pass_HTF = np.hstack([pass_HTF,0.0])
        
        # temperature development of the sorbent
        dTdx = -np.diff(self.T_HTF)/self.dx * self.HTF_flow*997*self.cp_water  # W/△m
        
        dTwdt = (pass_HTF + dTdx) / (self.m_water/self.length_tube * self.cp_water) # W/m
        #log_output.log_any_msg('{}\n{}\n\n'.format(self.T_HTF-self.T_HTF_init, dT))
        return np.hstack([0.0,dTwdt])


    # summarize the time differential equations into ODE system
    def equations(self, t,var):
        # variables
        self.loading   = var[0]
        self.mof_T     = var[1]
        self.gas_T     = var[2]
        self.T_HTF     = var[3:3+self.n_discrete+2]
        #log_output.log_any_msg(('{},'*(len(var)+1)).format(t,*var))

        # update the system
        self.update_variables(t)
        
        # differential
        dmsordt = self.mass_balance_sor(t)
        dTsordt = self.en_balance_sor(t)
        dTgasdt = self.en_balance_gas(t)
        dTwdt   = self.en_balance_water(t)
        return np.hstack([dmsordt, dTsordt, dTgasdt, dTwdt])


    # solver of ODE system
    def solver(self):
        self.t = np.linspace(0,self.simulation_time,int(self.simulation_time/self.dt))
        # initial value in the ODE
        var0 = [self.loading_init, self.mof_T_init, self.gas_T_init] + self.T_HTF_init
        # solver
        sol = solve_ivp(self.equations, [0, self.simulation_time], var0 ,t_eval=self.t,method='LSODA',rtol=1e-3, atol=1e-06,max_step=self.max_time_step)

        # time 
        self.t = sol.t
        
        # answer
        self.loading_list = sol.y[0,:]
        self.mof_T_list   = sol.y[1,:]
        self.gas_T_list   = sol.y[2,:]
        self.T_HTF_list   = sol.y[3+self.n_discrete,:]

        # calc variables
        #T_in_list  = self.find_output_variables()
        Total_obtained_H = self.out_heat_HTF()

        # plot graph
        self.plot_data = [self.loading_list,  self.gas_T_list, self.T_HTF_list, self.mof_T_list]#, T_in_list]
        self.legends = ['loading','temperature of gas in the tank','temperature of HTF out','temperature of sorbent']#, "temperature of input CO2"]
        #self.make_gragh(self.plot_data, self.legends)
        return Total_obtained_H
    

    
    # calculate the amount of heat that HTF obtain
    def out_heat_HTF(self):
        dQdt = self.HTF_flow * self.rhoCp_HTF * (self.T_HTF_list - self.T_HTF_in) # J
        Q_line = cumtrapz(dQdt, self.t, initial = 0)
        tot_Q = Q_line[-1]
        return tot_Q


    # plot the data, temperature and the others are distinguished 
    def make_gragh(self,data,legend):
        self.fig = plt.figure(figsize=(15,5))
        ax_T = self.fig.add_subplot(1,2,1)
        ax = self.fig.add_subplot(1,2,2)

        # plot for each data
        for i in range(len(data)):
            if 'temperature' in legend[i]:
                ax_T.plot(self.t, data[i] - 273.15, label = legend[i].replace('temperature of',''))
                ax_T.set_xlabel(' t sec')
                ax_T.set_ylabel(' T ℃')
                ax_T.legend()
                ax_T.grid(axis='both')
            else:
                ax.plot(self.t, data[i], label = legend[i])
                ax.set_xlabel(' t sec')
                ax.set_ylabel(' M mmol/g')
                ax.grid(axis='both')
                ax.legend()
        plt.savefig('./Fig/'+ self.name + '/results')
        return 0

        # re-calclate non-ODEsystem variables from the results
    def find_output_variables(self):
        W_comp = []
        T_in_list = []
        for i in range(len(self.t)):
            self.loading = self.loading_list[i]
            self.mof_T   = self.mof_T_list[i]
            self.gas_T   = self.gas_T_list[i]
            self.T_HTF   = self.T_HTF_list[i]
            self.update_variables(self.t[i])
            T_in_list.append(self.T_in)
            W_comp.append(self.compressor_work)
        return np.array(T_in_list)

if __name__ == '__main__':
    main()

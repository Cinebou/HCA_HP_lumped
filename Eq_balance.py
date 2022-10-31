from fluidProp import VLEFluid
from scipy.integrate import odeint,cumtrapz
import matplotlib.pyplot as plt
import numpy as np
from HC_AHP_System import MOF
import log_output 
from os import mkdir, path
from pprint import pprint

if not path.exists('./Fig'):
    mkdir('./Fig')


def main():
    mof_type = 'MIL-101'
    #mof_type = 'Uio-66'
    System = balance(mof_type)
    System.solver()


class balance(MOF):
    def __init__(self,name):
        super().__init__(name)

    # update system with each iteration
    def update_variables(self):
        # increase of the loading
        self.dm_loading = self.ads_speed_LDF()
        # increase of the adsorption amount
        self.m_ads = self.dm_loading * self.MOF_mass
        # the mass of the outflow gas, record m_out for the calucaltion in energy balance
        self.m_in = self.m_out + self.m_ads
        #energy balance of the HTF from the efficiency of exchange
        self.T_HTF_out = self.T_HTF_in + self.epsilonC / self.rhoCp_HTF / self.HTF_flow * (self.mof_T - self.T_HTF_in)


    # mass balance of the adsorbate
    def mass_balance_sor(self):
        # write down the log
        log_output.log_mass_msg(self.m_in, self.m_ads, self.m_out)
        return self.dm_loading

    
    # energy balance of the gas in the tank
    def en_balance_gas(self):
        # enthalpy of the CO2 flow
        h_in_CO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.T_in).h  * self.molar_mass # J / mol_CO2
        # refprop data of the gas in the tank
        PropCO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T)
        h_out_CO2 = PropCO2.h  * self.molar_mass # J / mol_CO2

        # energy flow from the CO2 stream
        en_flow = self.m_in * (h_in_CO2 - h_out_CO2)

        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.mof_T - self.gas_T)

        # expansion heat due to the adsorption
        mRT = - self.m_ads * self.gas_P * PropCO2.vol * self.molar_mass

        # temperature development of the gas
        cp_gas = PropCO2.cp * self.molar_mass
        dTgasdt = (en_flow + Htrans + mRT) / (self.m_gas*cp_gas)

        # Logger
        log_output.log_h_gas_msg(dTgasdt*self.m_in*cp_gas, en_flow, Htrans, mRT)
        log_output.log_any_msg('{},{},{}'.format(h_in_CO2, h_out_CO2, self.gas_T ))
        return dTgasdt
    

    # energy balance of the adsorbent
    def en_balance_sor(self):
        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.gas_T - self.mof_T)
        # adsorption heat
        Hads = self.m_ads * self.dHadsdm()
        # heat passed to the HTF
        pass_HTF = self.HTF_flow * self.rhoCp_HTF * (self.T_HTF_in - self.T_HTF_out)
    
        # temperature development of the sorbent
        dTsordt = (Htrans + Hads + pass_HTF) / self.cp_sor

        # Logger
        log_output.log_h_sor_msg(dTsordt*self.cp_sor, Htrans, Hads, pass_HTF)
        return dTsordt


    # summarize the time differential equations into ODE system
    def equations(self, var, t):
        # variables
        self.loading   = var[0]
        self.mof_T     = var[1]
        self.gas_T     = var[2]

        # update the system
        self.update_variables()
        
        # differential
        dmsordt = self.mass_balance_sor()
        dTsordt = self.en_balance_sor()
        dTgasdt = self.en_balance_gas()
        
        print('current time is  : ', t)
        return [dmsordt, dTsordt, dTgasdt]


    # solver of ODE system
    def solver(self):
        # time step
        self.t = np.linspace(0,self.simulation_time,self.simulation_time * 10)
        # initial value in the ODE
        var0 = [self.loading_init, self.mof_T_init, self.gas_T_init]
        # solver
        sol = odeint(self.equations, var0, self.t,atol=1.e-9,rtol=1.e-9)
        
        # answer
        self.loading_list = sol[:,0]
        self.mof_T_list   = sol[:,1]
        self.gas_T_list   = sol[:,2]

        # calc variables
        self.T_HTF_out_list  = self.find_output_variables()
        Total_obtained_H = self.out_heat_HTF()

        # plot graph
        plot_data = [self.loading_list,  self.mof_T_list, self.T_HTF_out_list, self.gas_T_list]
        legends = ['loading','temperature of sorbent','temperature of HTF out', 'temperature of gas in the tank']
        self.make_gragh(plot_data, legends)
        return 0
    
    # re-calclate non-ODEsystem variables from the results
    def find_output_variables(self):
        T_HTF_list = []
        for i in range(len(self.t)):
            self.loading = self.loading_list[i]
            self.mof_T   = self.mof_T_list[i]
            self.update_variables()
            T_HTF_list.append(self.T_HTF_out)
        return np.array(T_HTF_list)

    
    # calculate the amount of heat that HTF obtain
    def out_heat_HTF(self):
        dQdt = self.HTF_flow * self.rhoCp_HTF * (self.T_HTF_out_list - self.T_HTF_in) / 1000 # kJ
        tot_Q = cumtrapz(dQdt, self.t, initial = 0)

        dQdt_trans = self.h_CO2 * self.surfaceA * (self.mof_T_list - self.gas_T) /1000
        tot_H_gas = cumtrapz(dQdt_trans, self.t, initial=0)

        tot_H_ads = self.MOF_mass*(self.loading_list[-1] - self.loading_list[0]) * self.dHadsdm() /1000 # kJ
        
        print('total output heat at time       : ',self.simulation_time,' is   ', tot_Q[-1], 'kJ')
        print('total adsorption heat at time   : ',self.simulation_time,' is   ', tot_H_ads, 'kJ')
        print('total heat transfer CO2 at time : ',self.simulation_time,' is   ', tot_H_gas[-1], 'kJ')
        
        plt.rcParams["font.size"] = 14
        plt.plot(self.t, tot_Q)
        plt.xlabel('t sec')
        plt.ylabel('accumurative heat kJ')
        plt.savefig('./Fig/'+self.name + '/accumurative_heat.png')
        return tot_Q[-1]


    # plot the data, temperature and the others are distinguished 
    def make_gragh(self,data,legend):
        self.fig = plt.figure(figsize=(15,5))
        plt.rcParams["font.size"] = 14
        ax_T = self.fig.add_subplot(1,2,1)
        ax = self.fig.add_subplot(1,2,2)

        # plot for each data
        for i in range(len(data)):
            if 'temperature' in legend[i]:
                ax_T.plot(self.t, data[i], label = legend[i])
                ax_T.set_xlabel(' t sec')
                ax_T.set_ylabel(' T K')
                ax_T.legend()
                ax_T.grid(axis='both')
            else:
                ax.plot(self.t, data[i], label = legend[i])
                ax.set_xlabel(' t sec')
                ax.set_ylabel(' M mmol/g')
                ax.grid(axis='both')
                ax.legend()
        plt.savefig('./Fig/'+ self.name + '/results')
        plt.show()
        return 0


if __name__ == '__main__':
    main()

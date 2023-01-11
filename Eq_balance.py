from scipy.integrate import odeint,cumtrapz
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
    res_file = './exp_data/res/GL840_01_No2_2022-12-06_15-12-36.csv';start_sec=750
    
    System = balance(mof_type,res_file,start_sec)
    scp, totQ, tot_H_ads = System.solver()

    print('SCP of the system is ::  ',scp, '  W')
    print('Total hot water of the system is :: ', totQ,' J')
    print('Total generated adsorption heat is :: ', tot_H_ads,' J')

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
        # the mass of the outflow gas, record m_out for the calucaltion in energy balance
        self.m_in = self.m_out + self.m_ads
        # heat capacity of the sorbent
        self.cp_sor = self.calc_heat_Cap_sor()
        # compressor work
        self.compressor_work = 0 #self.h_comp * self.m_in * self.dt
        #get input temperature of CO2 from experiment
        self.T_in = float(self.exp_T_in.iloc[[int(t)]]) + 273.15; assert 270 <= self.T_in <= 400


    # mass balance of the adsorbate
    def mass_balance_sor(self,t):
        # write down the log
        log_output.log_mass_msg(self.m_in, self.m_ads, self.m_out,t)
        return self.dm_loading

    
    # energy balance of the gas in the tank
    def en_balance_gas(self,t):
        # enthalpy of the CO2 flow
        h_in_CO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.T_in).h  * self.molar_mass # J / mol_CO2
        # refprop data of the gas in the tank
        PropCO2 = self.gas.calc_fluidProp_pT(self.gas_P, self.gas_T)
        #PropCO2 = self.gas.calc_fluidProp_pT(self.gas_P, 303.15)
        h_out_CO2 = PropCO2.h  * self.molar_mass # J / mol_CO2

        # energy flow from the CO2 stream
        en_flow = self.m_in * h_in_CO2 - self.m_out * h_out_CO2
        #en_flow = self.m_in * (h_in_CO2 - h_out_CO2)

        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.mof_T - self.gas_T)

        # expansion heat due to the adsorption
        mRT = - self.m_ads * self.gas_P * PropCO2.vol * self.molar_mass

        # heat loss
        Q_loss = self.h_CO2 * self.A_vessel * (self.vessel_T - self.gas_T)

        # temperature development of the gas
        cp_gas = PropCO2.cp * self.molar_mass
        dTgasdt = (en_flow + Htrans + mRT + Q_loss) / (self.m_gas*cp_gas)

        # Logger
        #log_output.log_h_gas_msg(dTgasdt*self.m_in*cp_gas, en_flow, Htrans, mRT,t)
        #log_output.log_any_msg('{},{},{}, {}'.format(en_flow, Htrans, mRT, Q_loss))
        return dTgasdt
    

    # energy balance of the adsorbent
    def en_balance_sor(self,t):
        # heat transfer between sorbent and gas
        Htrans = self.h_CO2 * self.surfaceA * (self.gas_T - self.mof_T)
        # adsorption heat
        Hads = self.m_ads * self.dHadsdm()
        # heat passed to the HTF
        pass_HTF = self.h_water * self.surfaceA_HTF * (self.T_HTF - self.mof_T)
    
        # temperature development of the sorbent
        dTsordt = (Htrans + Hads + pass_HTF) / self.cp_sor

        # Logger
        #log_output.log_h_sor_msg(dTsordt*self.cp_sor, Htrans, Hads, pass_HTF,t)
        return dTsordt


    # energy balance of the water
    def en_balance_water(self,t):
        # heat transfer between sorbent and gas
        pass_HTF = self.h_water * self.surfaceA_HTF * (self.mof_T - self.T_HTF)
        # energy input and output
        enflow = self.cp_water * (self.HTF_flow*997) * (self.T_HTF_in - self.T_HTF)
        # temperature development of the sorbent
        dTwdt = (pass_HTF + enflow) / (self.m_water*self.cp_water)

        #log_output.log_any_msg('{}, {}'.format(pass_HTF, enflow))
        return dTwdt



    # summarize the time differential equations into ODE system
    def equations(self, var, t):
        # variables
        self.loading   = var[0]
        self.mof_T     = var[1]
        self.gas_T     = var[2]
        self.T_HTF     = var[3]

        # update the system
        self.update_variables(t)
        
        # differential
        dmsordt = self.mass_balance_sor(t)
        dTsordt = self.en_balance_sor(t)
        dTgasdt = self.en_balance_gas(t)
        dTwdt   = self.en_balance_water(t)
        return [dmsordt, dTsordt, dTgasdt, dTwdt]


    # solver of ODE system
    def solver(self):
        # time step
        self.t = np.linspace(0,self.simulation_time,int(self.simulation_time / self.dt))
        
        # initial value in the ODE
        var0 = [self.loading_init, self.mof_T_init, self.gas_T_init, self.T_HTF]
        # solver
        sol = odeint(self.equations, var0, self.t,atol=5.e-6,rtol=5.e-6)
        
        # answer
        self.loading_list = sol[:,0]
        self.mof_T_list   = sol[:,1]
        self.gas_T_list   = sol[:,2]
        self.T_HTF_list   = sol[:,3]

        # calc variables
        #self.T_HTF_list, W_comp_list, T_in_list  = self.find_output_variables()
        Total_obtained_H, tot_H_ads = self.out_heat_HTF()

        # final output
        SCP = Total_obtained_H / self.simulation_time  # kJ / sec
        #COP = Total_obtained_H / sum(W_comp_list)
        
        # plot graph
        self.plot_data = [self.loading_list,  self.gas_T_list, self.T_HTF_list, self.mof_T_list]#, T_in_list]
        self.legends = ['loading','temperature of gas in the tank','temperature of HTF out','temperature of sorbent']#, "temperature of input CO2"]
        #self.make_gragh(self.plot_data, self.legends)
        return SCP, Total_obtained_H, tot_H_ads
    
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
        return np.array(W_comp), np.array(T_in_list)

    
    # calculate the amount of heat that HTF obtain
    def out_heat_HTF(self):
        dQdt = self.HTF_flow * self.rhoCp_HTF * (self.T_HTF_list - self.T_HTF_in) # J
        Q_line = cumtrapz(dQdt, self.t, initial = 0)
        tot_Q = Q_line[-1]
        tot_H_ads = self.MOF_mass*(self.loading_list[-1] - self.loading_list[0]) * self.dHadsdm() # J
        return tot_Q, tot_H_ads

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


if __name__ == '__main__':
    main()

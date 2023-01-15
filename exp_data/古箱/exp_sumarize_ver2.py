import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
from data_index import index_data_ver3
from os import mkdir, path
from math import ceil, sqrt
from fluidProp import VLEFluid


if not path.exists('./Results'):
    mkdir('./Results/')
    mkdir('./Results/Fig/')
plt.rcParams["font.size"] = 24
plt.rcParams["font.family"] = "Times New Roman"


def main():
    #
    file_name ='./2022-12-26/GL840_01_No3_2022-12-26_15-44-03.csv'; start_sec=18000; dir_name = 'expA_300g_eg'
    calc_exp(file_name,start_sec, dir_name)


# file 
def calc_exp(file_name, start_sec, dir_name,exp_span = 300): 
    res1 = ResData(file_name, start_sec, exp_span, dir_name)
    res1.T_graph()
    res1.P_graph()
    res1.work_graph()
    res1.graph_each("${\dot m_{CO2}}$",'[mol/sec]')#, ylim=[0, 0.2])
    res1.graph_each('expected_h_ads', '[W]')# , ylim=[-300, 300])
    res1.graph_each("${\dot v_{water}}$", '[ml/sec]', ylim=[0, 10.5])
    res1.graph_each("MF3-Water2", '[l/min]')
    


class ResData():
    def __init__(self,f_name, start_sec, exp_span, dir_name):
        print(" : ",dir_name.split('/').pop())
        self.f_name = f_name
        self.exp_data = pd.read_csv(f_name, skiprows=88, encoding='Shift-Jis',low_memory=False) # new data
        self.data_index = index_data_ver3()

        # length of the experimental table
        self.tf_all = len(self.exp_data) - 1

        # properties
        self.co2 = VLEFluid('CO2')
        self.cp_water = 4182 # J/kg/K
        self.mass_adsorbent = 0.1349 # kg
        self.cp_mof = self.mass_adsorbent * 658.9 # J/K
        self.cp_HX = 350 # J/K
        self.cp_vessel = 3883 # J/K
   


        # rough start second of the experiment
        self.ts_record = start_sec  # sec
        self.exp_span = exp_span # duration of the experiment    
        self.finish_t = start_sec + exp_span
        # time
        time = self.exp_data.loc[1:, self.data_index['time']]
        self.time_line = [(datetime.datetime.strptime(time[i], '%Y-%m-%d %H:%M:%S') -
                           datetime.datetime.strptime(time[1], '%Y-%m-%d %H:%M:%S')).seconds for i in range(1, len(time)+1)]
        # accurately calculate start time
        self.__start_and_finish()
        # data calibration
        self.__data_calibrate()
        # data arrangement
        self.__data_arrange()
        # experimental condition
        self.__exp_condition()

        # result derectory
        print(dir_name)
        self.file_path = './Results/Fig/' + dir_name+'_'+self.f_name.replace('.csv',"").replace('.CSV',"").split('GL840_01_No').pop()
        #if not path.exists('./Results/Fig/' + dir_name.split('/').pop(0)):
            #mkdir('./Results/Fig/' + dir_name.split('/').pop(0))
        if not path.exists(self.file_path):
            mkdir(self.file_path)
        print('\n\n\n')


    # thermocouple calibration, water calibrated at 2022/12/15, co2 at 2022/12/20
    def __data_calibrate(self):
        self.exp_data["CO2 outlet"] = self.d_loc_time("T26-co2 outlet",0,self.tf_all) * 1.0021 - 0.0199
        self.exp_data["CO2 inlet"] = self.d_loc_time("T25-co2 inlet",0,self.tf_all)  * 1.0038 - 0.1004
        self.exp_data["water inlet"] = self.d_loc_time("T11-tank1 w inlet",0,self.tf_all) * 1.0072 - 0.6292
        self.exp_data["water outlet"] = self.d_loc_time("T12-tank1 w outlet",0,self.tf_all) * 1.0057 - 0.2134


    # data setting arrangement
    def __data_arrange(self):
        self.exp_data["MF1-CO2 flow"] = self.d_loc_time('MF1-CO2 flow_',0,self.tf_all)
        self.exp_data['MF2-Water1'] = self.d_loc_time('MF2-Water1_',0,self.tf_all)
        self.exp_data["${\dot m_{CO2}}$"] = self.exp_data["MF1-CO2 flow"]*1000/44/3600
        self.exp_data["${\dot v_{water}}$"] = self.d_loc_time("MF2-Water1",0,self.tf_all) / 60 * 1000 # ml/sec

        # average temperature in the tank
        self.exp_data['sorbent'] = self.d_average(["T21-tank1 mof", "T22-tank1 mof", "T23-tank1 mof"])

        # average temperature in the tank
        self.exp_data['T_vessel CO2 inlet'] = self.d_average(["T51-tank1 vessel CO2in"])
        self.exp_data['T_vessel bottom'] = self.d_average(["T52-tank1 vessel bottom","T53-tank1 vessel bottom","T54-tank1 vessel bottom","T55-tank1 vessel bottom"])
        self.exp_data['T_vessel middle'] = self.d_average(["T56-tank1 vessel middle","T57-tank1 vessel middle","T58-tank1 vessel middle","T59-tank1 vessel middle"])
        self.exp_data['T_vessel upper'] = self.d_average(["T60-tank1 vessel upper","T61-tank1 vessel upper","T62-tank1 vessel upper","T63-tank1 vessel upper"])
        self.exp_data['T_vessel top'] = self.d_average(["T64-tank1 vessel top"])

        # average vessel temperature of all
        vessel_temps = ["T52-tank1 vessel bottom","T56-tank1 vessel middle","T60-tank1 vessel upper","T64-tank1 vessel top"]
        self.exp_data["vessel"] = self.d_loc_time('T_vessel bottom',0,self.tf_all)*0.07128 \
                                         +self.d_loc_time('T_vessel upper',0,self.tf_all)* 0.32886 + self.d_loc_time('T_vessel top',0,self.tf_all)*0.59986

        # average gas temp
        self.exp_data["gas"] = self.d_average(["T24-tank1 gas low","T29-tank1 gas high"])

        # gas properties
        # high pressure of the tank
        self.p_high = self.d_loc_time("Inlet", self.start_t+20, self.finish_t).mean()
        # initial temperature of the gas
        self.T_init_gas = self.d_loc_time("gas" , self.start_t-1, self.start_t).mean()
        # specific heat capacity of co2
        self.cp_CO2 = self.co2.calc_fluidProp_pT(self.p_high*1e6, self.T_init_gas+273.15).cp # J/kg/K
        self.m_CO2_in_tank = self.co2.calc_fluidProp_pT(self.p_high*1e6, self.T_init_gas+273.15).vol

        # water output
        self.exp_data['dT_water1'] = self.d_loc_time('water outlet',0,self.tf_all) - self.d_loc_time('water inlet',0,self.tf_all)
        self.exp_data['$Q_{water}$'] = self.exp_data['dT_water1'] * self.d_loc_time('MF2-Water1',0,self.tf_all) * self.cp_water/60 # W 

        # CO2 flow
        self.exp_data['dT_co2_flow'] = self.d_loc_time("CO2 inlet",0,self.tf_all) - self.d_loc_time("CO2 outlet",0,self.tf_all)
        self.exp_data['$Q_{CO2}$'] = self.exp_data['dT_co2_flow'] * self.d_loc_time('MF1-CO2 flow',0,self.tf_all) * self.cp_CO2/3600 # W, cp_co2 = 1.0 kJ/kg/h

        # mof heat cap, temperature difference by time
        self.exp_data['Q_loss mof'] = self.diff_times('sorbent') * (self.cp_mof + self.cp_HX)
        self.exp_data['Q_loss vessel'] = (self.diff_times("T_vessel bottom")*0.07128 +self.diff_times('T_vessel upper')* 0.32886 \
                                          + self.diff_times('T_vessel top')*0.59986) * self.cp_vessel # J/sec
        self.exp_data['Q_gas_change'] = self.diff_times("gas") * self.cp_CO2     # W = K/s * kg * J/K/kg

        """
        Whole tank vessel = 113.36 m3,
        bottom sphere = 8.03 e-5 m3, cylinder = 37.28 e-5 m3, flange = 68.00 e-5 m3
        ratio [bottom, cylinder, flange] = [0.07128, 0.32886, 0.59986]
        """

        # expected heat of adsorption from heat balance including heat losses at gas, vessel, HX     
        self.exp_data['expected_h_ads'] = self.exp_data['$Q_{water}$'] + self.exp_data['Q_loss mof'] + self.exp_data['Q_loss vessel']\
                                         +self.exp_data['Q_gas_change'] - self.exp_data['$Q_{CO2}$']

        # uncertainity analysis
        self.uncertains()

        

    # findout the starting point of the experiment
    def __start_and_finish(self):
        p_tank = self.d_loc_time("Inlet", self.ts_record, self.finish_t)
        self.start_t = self.ts_record
        if not self.start_t == 0:
            for i in range(3, len(p_tank)):
                if p_tank.iloc[i] > p_tank.iloc[i-1] * 1.04 + 0.1: # if there is a sudden increase of pressure
                    self.start_t = self.ts_record + i-1  # start moment of experiment is decided
                    self.ts_record = self.start_t - 30
                    break
        self.finish_t = min(self.start_t + self.exp_span, self.tf_all)
        self.time_line = np.array(self.time_line[self.ts_record:self.finish_t]) - self.start_t


    # data extract
    def d_loc(self, data_name):
        if data_name in self.data_index:
            name = self.data_index[data_name]
        else:
            name = data_name
        data = self.exp_data.loc[1:, name]\
            .replace('BURNOUT', None).replace('+++++++', None).replace('-------', None).replace('*******',None).astype(float)
        return data[self.ts_record : self.finish_t]


    # data extract
    def d_loc_time(self, data_name, ts, tf):
        if data_name in self.data_index:
            name = self.data_index[data_name]
        else:
            name = data_name
        data = self.exp_data.loc[1:, name]\
            .replace('BURNOUT', None).replace('+++++++', None).replace('-------', None).replace('*******',None).astype(float)
        return data[ts: tf]


    # take average over several data set
    def d_average(self, data_names):
        d_ave = np.empty_like(self.d_loc_time(data_names[0],0,self.tf_all))
        for data in data_names:
            d_ave += self.d_loc_time(data,0, self.tf_all)
        d_ave /= len(data_names)
        return d_ave


    # moving average over time
    def moving_average(self,data):
        n_conv = 15 # seconds average
        b = np.ones(n_conv)/n_conv
        m_time_ave = self.valid_convolve(data, n_conv)
        return np.append(0,m_time_ave)


    # diff time
    def diff_times(self, data_name):
        diff = np.diff(self.d_loc_time(data_name,0,self.tf_all))
        return np.concatenate([[0,0], diff])


    # combined standard uncertainity
    def __combined_uncertain(self,m_dot, t1, t2, sigma_m, sigma_t1, sigma_t2):
        m_variance = sigma_m**2 / self.exp_data[m_dot]**2
        T_variance = (sigma_t1**2 + sigma_t2**2) / (self.exp_data[t1] - self.exp_data[t2])**2
        combined_std = (m_variance + T_variance)**0.5
        return combined_std


    def uncertains(self):
        # standard deviation of thermo couples, calibrated at 2022/12/15
        self.u_T11 = 0.08089 # K
        self.u_T12 = 0.06945 # K
        self.u_T25 = 0.099259 # K
        self.u_T26 = 0.09641  # K
        # catalog accuracy of flow meter, assuming rectangle distribution
        self.u_m_co2 = 1.2 /sqrt(3) # kg/h
        self.u_m_water = 0.040 /sqrt(3) # l/min, 定格5 l/minの繰り返し精度0.8 % F.S. = 0.040 l/min, 30sec出力
        # expansion coefficient for uncertainity to include 95 % of possibility
        expansion_coeff = 2
        self.exp_data['Q_CO2_uncertain'] = self.exp_data['$Q_{CO2}$'] *expansion_coeff* self.__combined_uncertain('MF1-CO2 flow','CO2 inlet', 'CO2 outlet',self.u_m_co2, self.u_T25, self.u_T26)
        self.exp_data['Q_w_uncertain'] = self.exp_data['$Q_{water}$'] *expansion_coeff* self.__combined_uncertain('MF2-Water1','water inlet', 'water outlet',self.u_m_water, self.u_T11, self.u_T12)


    # graph. work
    def work_graph(self):
        fig_work = plt.figure(figsize=(10.0, 6.0))
        ax = fig_work.add_subplot(111)
        ax.set_position([0.12,0.15,0.82,0.8])
        d_names = ["$Q_{water}$",
                   "$Q_{CO2}$",
                   #'Q_loss vessel',
                   #'Q_loss mof',
                   #'Q_gas_change'
                  ]
        for name in d_names:
            data = self.d_loc(name)
            ax.plot(self.time_line, data, label=name, linewidth=5)
            #print(name, ' total  : ', sum(self.d_loc_time(name,self.start_t,self.finish_t)), ' J')
            ax.set_xlim(0,len(data))

        ax.set_xlabel('time [sec]')
        ax.set_ylabel('ΔQ [W]')
        ax.legend(bbox_to_anchor=(0.7, 0.95), loc='upper left', borderaxespad=0)
        ax.grid()
        fig_work.savefig(self.file_path +  '/work.png')


    # temperature in teh tank
    def T_graph(self):
        fig = plt.figure(figsize=(12.0, 6.0))
        ax = fig.add_subplot(111)
        ax.set_position([0.1,0.15,0.6,0.8])
        d_names = ["CO2 inlet",
                   "CO2 outlet",
                   "water outlet",
                   #"gas",
                   "vessel",
                   #"T2-compressor outlet",
                   #"T7-WHX1 comp",
                   #"T8-WHX1 expV",
                   #"water inlet",
                   #"T21-tank1 mof",
                   #"T22-tank1 mof",
                   #"T23-tank1 mof",
                   #"T24-tank1 gas low",
                   #"T29-tank1 gas high",
                   'sorbent', 
                   #"T-room",
                   #'T_vessel CO2 inlet',
                   #'T_vessel bottom',
                   #'T_vessel middle',
                   #'T_vessel upper',
                   #"T_vessel top"
                  ]
        
        max_temp = 30
        for name in d_names:
            data = self.d_loc(name)
            #max_temp = max(max_temp, max(data))
            name = name.replace('CO2', 'CO$_2$')
            ax.plot(self.time_line, data, label=name, linewidth=3)

        ax.set_xlabel('time [sec]')
        ax.set_ylabel('temperature [℃]')
        #ax.set_ylim(22, max_temp+4)
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
        ax.grid()
        fig.savefig(self.file_path + '/temp.png')


    # Pressure in the tank
    def P_graph(self):
        fig = plt.figure(figsize=(9.0, 6.0))
        ax = fig.add_subplot(111)
        ax.set_position([0.16,0.15,0.8,0.8])
        d_names = [#"PS1-compressor inlet",
                   #"PS2-compressor outlet",
                   "Inlet",
                   "Outlet",
                   #"PS7-WHX1 comp",
                   #"PS8-WHX1 expV"
                  ]
        for name in d_names:
            data = self.d_loc(name)
            ax.plot(self.time_line, data, label=name, linewidth=3)  
        ax.set_xlabel('time [sec]')
        ax.set_ylabel('pressure [MPa]')
        ax.grid(axis="x")
        #ax.minorticks_on()
        ax.grid(which = "both", axis="y")
        ax.set_ylim(0, max(data)+1.0)
        ax.legend(bbox_to_anchor=(0.65, 0.22), loc='upper left', borderaxespad=0)
        fig.savefig(self.file_path + '/pressure.png')


    # calculate experimental condition
    def __exp_condition(self):
        #print("experiment start sec is : ", self.start_t,', ', self.exp_data.loc[self.start_t, self.data_index['time']])
        #print("experiment finish sec is : ", self.finish_t,', ', self.exp_data.loc[self.finish_t, self.data_index['time']])

        p_init = self.d_loc_time("Inlet", self.start_t-60, self.start_t-10).mean()
        T_inlet_HTF = self.d_loc_time("water inlet", self.start_t, self.finish_t).mean()
        CO2_flow = self.d_loc_time("MF1-CO2 flow", self.start_t+100, self.finish_t).mean()
        water_flow = self.d_loc_time("MF2-Water1",  self.start_t, self.finish_t).mean()
        T_init_mof = self.d_loc_time('sorbent', self.start_t-20, self.start_t-10).mean()
        self.T_init_gas = self.d_loc_time("gas" , self.start_t-1, self.start_t).mean()
        self.T_fin_gas = self.d_loc_time("gas" , self.finish_t-1, self.finish_t).mean()
        T_fin_mof = self.d_loc_time('sorbent', self.finish_t-5, self.finish_t).mean()
        T_init_vessel = self.d_loc_time("vessel", self.start_t-5, self.start_t-1).mean()
        T_fin_vessel = self.d_loc_time("vessel", self.finish_t-1, self.finish_t).mean()
        T_room = self.d_loc("T-room").mean()
        Q_loss_vessel = (T_fin_vessel - T_init_vessel) * self.cp_vessel
        Q_HX = (T_fin_mof - T_init_mof) * self.cp_mof
        self.en_CO2 = sum(self.d_loc_time('$Q_{CO2}$',self.start_t,self.finish_t))
        self.en_water = sum(self.d_loc_time('$Q_{water}$',self.start_t,self.finish_t))
        Q_residual = self.en_water + Q_HX - self.en_CO2 + Q_loss_vessel 
        # uncertainity analysis
        Uc_Q_CO2 = sum(self.d_loc_time('Q_CO2_uncertain',self.start_t,self.finish_t))
        Uc_water_CO2 = sum(self.d_loc_time('Q_w_uncertain',self.start_t,self.finish_t))

        # calculated adsorption amount with DA model
        ads_amount_init = self.uio66.eq_loading(p_init*1e6, T_init_mof+273.15)
        ads_amount_fin = self.uio66.eq_loading(self.p_high*1e6, T_fin_mof+273.15)
        d_ads = ads_amount_fin - ads_amount_init
        h_ads = d_ads * self.mass_adsorbent * self.uio66.heat_coeff * 1000 # J

        ###output
        print('Initial pressure in the tank  : ',p_init, ' MPa')
        print('High pressure in the tank     : ',self.p_high, ' MPa')
        print('temperature of HTF at the inlet of the tank  : ', T_inlet_HTF, ' C')
        print('Frow rate of CO2 at constant flow  : ', CO2_flow, ' kg/h')
        print('Frow rate of water in tank 1  : ', water_flow, ' l/min')
        print('initial temperature of the mof  : ', T_init_mof, ' C')
        print('initial temperature of the gas in the tank  : ', self.T_init_gas, ' C')
        print('final temperature of the gas in the tank  : ', self.T_fin_gas, ' C')
        print('final temperature of mof  : ', T_fin_mof, ' C')
        print('initial temperature of the vessel  : ', T_init_vessel, ' C')
        print('final temperature of the vessel  : ', T_fin_vessel, ' C')
        print("temperatrue of the room  : ", T_room, ' C')
        print('Heat loss to the vessel is :', Q_loss_vessel, ' J')
        print('Heat loss to the HX : ', Q_HX, ' J')
        print('Input heat CO2 : ', self.en_CO2, ' J')
        print('output heat water : ', self.en_water, ' J')
        print('Residual heat : ', Q_residual, ' J')
        print('Uncertainity of CO2 is : ', Uc_Q_CO2, ' J')
        print('Uncertainity of water is : ', Uc_water_CO2, ' J')

        print('expected adsorption amount initial : ',ads_amount_init, ' mol/kg' )
        print('expected adsorption amount final   : ',ads_amount_fin , ' mol/kg' )
        print('expected adsorption amount diff    : ',d_ads          , ' mol/kg' )
        print('expected adsorption heat           : ',h_ads          , ' J' )
      

    # output each data
    def graph_each(self, data_name, data_unit, ylim=[None,None]):
        data = self.d_loc(data_name)
        fig_ = plt.figure()
        ax = fig_.add_subplot(111)
        ax.plot(self.time_line, data, linewidth=3)  
        ax.set_xlabel('time [sec]')
        ax.set_ylabel(data_name+'  '+ data_unit)
        ax.set_position([0.25,0.15,0.7,0.8])
        ax.grid()
        ax.set_ylim(ylim[0], ylim[1])
        fig_.savefig(self.file_path + '/' + data_name.replace('\\','') + '.png')


    # remove noise, low pass filter
    def low_pass(self, data):
        fc = 1.0  # カットオフ周波数

        N =  len(data); dt = 1.0  #number of samples; #sampling rate
        t = np.arange(0, N*dt, dt); freq = np.linspace(0, 1.0/dt, N) # 時間軸 周波数軸
        F = np.fft.fft(data)# 高速フーリエ変換（周波数信号に変換）
        F = F/(N/2)  # 正規化 + 交流成分2倍
        F[0] = F[0]/2
        F2 = F.copy()  # 配列Fをコピー
        F2[(freq > fc)] = 0  # ローパスフィル処理（カットオフ周波数を超える帯域の周波数信号を0にする）
        f2 = np.fft.ifft(F2) # 高速逆フーリエ変換（時間信号に戻す）
        f2 = np.real(f2*N) # 振幅を元のスケールに戻す
        return np.append(0,f2)

    # xxに対してsize個での移動平均を取る
    def valid_convolve(self, xx, size):
        b = np.ones(size)/size
        xx_mean = np.convolve(xx, b, mode="same")
        n_conv = ceil(size/2)
        # 補正部分
        xx_mean[0] *= size/n_conv
        for i in range(1, n_conv):
            xx_mean[i] *= size/(i+n_conv)
            xx_mean[-i] *= size/(i + n_conv - (size % 2)) 
	    # size%2は奇数偶数での違いに対応するため
        return xx_mean


if __name__ == '__main__':
    main()

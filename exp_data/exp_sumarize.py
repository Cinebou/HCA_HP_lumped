import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
from data_index import index_data, index_data_previous
from os import mkdir, path
import seaborn as sns
#from fluidProp import VLEFluid
from pprint import pprint
if not path.exists('./Results'):
    mkdir('./Results/')
    mkdir('./Results/Fig/')
    mkdir('./Results/temp_data/')
plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "Times New Roman"


def main():
    file_name = './res/GL840_01_No6_2022-11-10_16-25-01.CSV'; start_sec=630
    #file_name = './res/GL840_01_no5_2022-11-10_16-20-31.CSV'; start_sec=4600
    
    exp_span = 424
    res1 = ResData(file_name, start_sec, exp_span)
    #res1.enthalpy_temp("T3-tank1 comp")

    
    res1.T_graph()
    #res1.P_graph()
    #res1.work_graph()
    #res1.graph_each("MF1-CO2 flow",'[kg/h]')
    #res1.graph_each('d_inout_en', '[W]')
    #res1.graph_each("MF2-Water1", '[l/min]', ylim=[0, 0.8])
    #res1.graph_each("MF3-Water2", '[l/min]')
    #plt.show()
    
    


class ResData():
    def __init__(self,f_name, start_sec, exp_span):
        self.f_name = f_name
        self.exp_data = pd.read_csv(f_name, skiprows=58, encoding='Shift-Jis') # new data
        self.data_index = index_data()
        #self.exp_data = pd.read_csv(f_name, skiprows=56, encoding='Shift-Jis')  # older data
        #self.data_index = index_data_previous()

        # length of the experimental table
        self.tf_all = len(self.exp_data) - 1
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

        # gas properties
        #self.gas = VLEFluid('CO2')
        # data arrangement
        self.__data_arrange()
        # experimental condition
        self.__exp_condition()

        # result derectory
        self.file_path  = './Results/Fig/' + self.f_name.replace('.CSV',"").split('GL840_01_').pop()
        if not path.exists(self.file_path):
            mkdir(self.file_path)
        
    
    # data setting arrangement
    def __data_arrange(self):
        # average temperature in the tank
        self.exp_data['T-tank1 ave'] = (self.d_loc_time("T21-tank1 mof",0,self.tf_all) + self.d_loc_time("T22-tank1 mof",0,self.tf_all) \
            + self.d_loc_time("T23-tank1 mof",0,self.tf_all) + self.d_loc_time("T24-tank1 mof",0,self.tf_all)) / 4
        #self.exp_data['T-tank2 ave'] = (self.d_loc_time("T25-tank2 mof",0,self.tf_all) + self.d_loc_time("T26-tank2 mof",0,self.tf_all) \
        #    + self.d_loc_time("T27-tank2 mof",0,self.tf_all) + self.d_loc_time("T28-tank2 mof",0,self.tf_all)) / 4

        # water output
        self.exp_data['dT_water1'] = self.d_loc_time('T12-tank1 water outlet',0,self.tf_all) - self.d_loc_time('T11-tank1 water inlet',0,self.tf_all)
        self.exp_data['Q_water_out'] = self.exp_data['dT_water1'] * self.d_loc_time('MF2-Water1',0,self.tf_all) * 4182/60 # W 

        # CO2 flow
        self.exp_data['dT_co2_flow'] = self.d_loc_time("T3-tank1 comp",0,self.tf_all) - self.d_loc_time("T4-tank1 expV",0,self.tf_all)
        self.exp_data['dh_co2_flow'] = self.exp_data['dT_co2_flow'] * self.d_loc_time('MF1-CO2 flow',0,self.tf_all) * 1000.0/3600 # W, cp_co2 = 1.0 kJ/kg/h

        self.exp_data['d_inout_en'] = self.exp_data['dh_co2_flow'] - self.exp_data['Q_water_out']
        # mof heat cap, temperature difference by time
        cp_HX = 350 # J/K
        self.exp_data['dhdt_mof'] = self.diff_times('T-tank1 ave') * (0.098 * 550 + cp_HX)
        # W = K/s * kg * J/K/kg    


    # findout the starting point of the experiment
    def __start_and_finish(self):
        p_tank = self.d_loc_time("PS3-tank1 comp", self.ts_record, self.finish_t)
        self.start_t = 0
        for i in range(3, len(p_tank)):
            if p_tank.iloc[i] > p_tank.iloc[i-3] * 1.5 + 0.15: # if there is a sudden increase of pressure
                self.start_t = self.ts_record + i-3  # start moment of experiment is decided
                break
        self.finish_t = self.start_t + self.exp_span
        print("experiment start sec is : ", self.start_t,', ', self.exp_data.loc[self.start_t, self.data_index['time']])
        print("experiment finish sec is : ", self.finish_t,', ', self.exp_data.loc[self.finish_t, self.data_index['time']])

        self.time_line = np.array(self.time_line[self.ts_record:self.finish_t]) - self.start_t


    # refprop heat capacity
    def enthalpy_temp(self,temp_name):
        h_list = []
        T_list = self.d_loc_time("T3-tank1 comp",0,self.tf_all)
        P_list = self.d_loc_time("PS3-tank1 comp",0,self.tf_all)
        len(T_list)
        len(P_list)
        for index, rowT in T_list.iteritems():
            h = self.gas.calc_fluidProp_pT(P_list.iloc[index] ,T_list.iloc[index]).h
            #print(h)
        return 0


    # data extract
    def d_loc(self, data_name):
        if data_name in self.data_index:
            name = self.data_index[data_name]
        else:
            name = data_name
        data = self.exp_data.loc[1:, name]\
            .replace('BURNOUT', None).replace('+++++++', None).replace('-------', None).astype(float)
        return data[self.ts_record : self.finish_t]


    # data extract
    def d_loc_time(self, data_name, ts, tf):
        if data_name in self.data_index:
            name = self.data_index[data_name]
        else:
            name = data_name
        data = self.exp_data.loc[1:, name]\
            .replace('BURNOUT', None).replace('+++++++', None).replace('-------', None).astype(float)
        return data[ts: tf]


    # diff time
    def diff_times(self, data_name):
        dx_list = np.array([0,0,0])
        for t in range(1,self.tf_all-1):
            dx = self.d_loc_time(data_name,1,self.tf_all).iloc[t] - self.d_loc_time(data_name,0,self.tf_all-1).iloc[t]
            dx_list = np.append(dx_list,dx)
        return dx_list


    # graph. work
    def work_graph(self):
        fig_work = plt.figure(figsize=(10.0, 6.0))
        ax = fig_work.add_subplot(111)
        ax.set_position([0.1,0.1,0.8,0.8])
        d_names = ["Q_water_out",
                   "dh_co2_flow",
                   #'dhdt_mof'
                  ]
        for name in d_names:
            data = self.d_loc(name)
            ax.plot(self.time_line, data, label=name)
            print(name, ' total  : ', sum(self.d_loc_time(name,self.start_t,self.finish_t)), ' J')

        ax.set_xlabel('time [sec]')
        ax.set_ylabel('ΔQ [W]')
        ax.legend(bbox_to_anchor=(0.68, 0.94), loc='upper left', borderaxespad=0)
        ax.grid()
        fig_work.savefig(self.file_path +  '/work.png')


    # temperature in teh tank
    def T_graph(self):
        fig = plt.figure(figsize=(10.0, 6.0))
        ax = fig.add_subplot(111)
        ax.set_position([0.1,0.1,0.55,0.8])
        d_names = [#"T2-compressor outlet",
                   "T3-tank1 comp",
                   "T4-tank1 expV",
                   #"T7-WHX1 comp",
                   #"T8-WHX1 expV",
                   #"T11-tank1 water inlet",
                   "T12-tank1 water outlet",
                   #"T21-tank1 mof",
                   #"T22-tank1 mof",
                   #"T23-tank1 mof",
                   #"T24-tank1 mof",
                   #"T25-tank2 mof",
                   #"T26-tank2 mof",
                   #"T27-tank2 mof",
                   #"T28-tank2 mof",
                   'T-tank1 ave',
                   "T-room"  
                  ]
        
        for name in d_names:
            data = self.d_loc(name)
            ax.plot(self.time_line, data, label=name)

        ax.set_xlabel('time [sec]')
        ax.set_ylabel('temperature [℃]')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        ax.grid()
        fig.savefig(self.file_path + '/temp.png')


    # Pressure in the tank
    def P_graph(self):
        fig = plt.figure(figsize=(10.0, 6.0))
        ax = fig.add_subplot(111)
        ax.set_position([0.1,0.1,0.8,0.8])
        d_names = [#"PS1-compressor inlet",
                   #"PS2-compressor outlet",
                   "PS3-tank1 comp",
                   "PS4-tank1 expV",
                   #"PS7-WHX1 comp",
                   #"PS8-WHX1 expV"
                  ]
        for name in d_names:
            data = self.d_loc(name)
            ax.plot(self.time_line, data, label=name)  
        ax.set_xlabel('time [sec]')
        ax.set_ylabel('pressure [MPa]')
        ax.legend(bbox_to_anchor=(0.65, 0.22), loc='upper left', borderaxespad=0)
        ax.grid()
        fig.savefig('./Results/Fig/'+self.f_name.replace('.CSV',"").split('GL840_01_').pop() + '/pres.png')


    # calculate experimental condition
    def __exp_condition(self):
        p_init = self.d_loc_time("PS3-tank1 comp", self.start_t-60, self.start_t-10).mean()
        print('Initial pressure in the tank  : ',p_init, ' MPa')
        p_high = self.d_loc_time("PS3-tank1 comp", self.start_t+10, self.finish_t).mean()
        print('High pressure in the tank     : ',p_high, ' MPa')
        T_inlet_CO2 = self.d_loc_time("T3-tank1 comp", self.start_t, self.finish_t).mean()
        print('temperature of CO2 at the inlet of the tank  : ', T_inlet_CO2, ' ℃')
        T_inlet_HTF = self.d_loc_time("T11-tank1 water inlet", self.start_t, self.finish_t).mean()
        print('temperature of HTF at the inlet of the tank  : ', T_inlet_HTF, ' ℃')
        CO2_flow = self.d_loc_time("MF1-CO2 flow", self.start_t, self.finish_t).mean()
        print('Frow rate of CO2 before compression  : ', CO2_flow, ' kg/h')
        water_flow = self.d_loc_time("MF2-Water1",  self.start_t, self.finish_t).mean()
        print('Frow rate of water in tank 1  : ', water_flow, ' l/min')
        T_init_mof = self.d_loc_time('T-tank1 ave', self.start_t-60, self.start_t-10).mean()
        print('initial temperature of the mof  : ', T_init_mof, ' ℃')
        T_init_gas = self.d_loc_time("T3-tank1 comp", self.start_t-10, self.start_t).mean()
        print('initial temperature of the gas in the tank  : ', T_init_gas, ' ℃')
        T_fin_tank = self.d_loc_time('T-tank1 ave', self.finish_t-2, self.finish_t).mean()
        print('final temperature of the tank  : ', T_fin_tank, ' ℃')
        T_room = self.d_loc("T-room").mean()
        print("temperatrue of the room  : ", T_room, ' ℃')
        Q_loss_equil = self.d_loc_time('d_inout_en', self.finish_t - 100, self.finish_t).mean()
        print('The average heat loss after reaching euilibrium is :', Q_loss_equil, ' W')



    # output each data
    def graph_each(self, data_name, data_unit, ylim=[None,None]):
        data = self.d_loc(data_name)
        fig_ = plt.figure()
        ax = fig_.add_subplot(111)
        ax.plot(self.time_line, data)  
        ax.set_xlabel('time [sec]')
        ax.set_ylabel(data_name+'  '+ data_unit)
        ax.set_position([0.15,0.15,0.8,0.8])
        ax.grid()
        ax.set_ylim(ylim[0], ylim[1])
        fig_.savefig(self.file_path + '/' + data_name + '.png')


if __name__ == '__main__':
    main()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import mkdir, path
from exp_sumarize import ResData
if not path.exists('./Results'):
    mkdir('./Results/')
    mkdir('./Results/Fig/')
    mkdir('./Results/compare')
plt.rcParams["font.size"] = 23
plt.rcParams["font.family"] = "Times New Roman"

def main():
    exp_span = 240
    files, times = file_list()
    data_all(exp_span, files, times)


def data_all(span, files, times,dir, legends = ['0.3 MPa','0.5 MPa','0.8 MPa']):
    fig_P = plt.figure(figsize=(9,6))
    ax_P = fig_P.add_subplot(111)
    fig_T = plt.figure(figsize=(9,6))
    ax_T = fig_T.add_subplot(111)
    fig_flow = plt.figure(figsize=(9,6))
    ax_f = fig_flow.add_subplot(111)
    fig_flow_co2 = plt.figure(figsize=(9,6))
    ax_f_co2 = fig_flow_co2.add_subplot(111)
    fig_W = plt.figure(figsize=(9,6))
    ax_w = fig_W.add_subplot(111)
    d_all = []

    for i in range(len(files)):
        data = data_collect(files[i],times[i],span)
        ax_P = graph_P(ax_P,data['timeline'], data['pressure'], legends[i])
        ax_f = graph_flow(ax_f,data['timeline'], data["water flow"], legends[i],'${\dot v_{w}}$ [ml/sec]',ylim=[0,12])
        ax_f_co2 = graph_flow(ax_f_co2,data['timeline'], data["CO2 flow"], legends[i],'${\dot m_{CO2}}$ [mol/sec]',ylim=[0,0.10])
        ax_T = graph_T(ax_T,data['timeline'], data['water temperature'], legends[i])
        ax_w = graph_W(ax_w,data['timeline'], data['Q_water'], legends[i])#+' $Q_{water}$')
        #ax_w = graph_W(ax_w,data['timeline'], data['Q_co2'], legends[i]+' $Q_{CO2}$')
        d_all.append(data['res'])

    if not path.exists('./Results/compare/'+dir):
        mkdir('./Results/compare/'+dir)
    fig_P.savefig('./Results/compare/'+dir+'pressure.png')
    fig_T.savefig('./Results/compare/'+dir+'temperature.png')
    fig_flow.savefig('./Results/compare/'+dir+'flowWater.png')
    fig_flow_co2.savefig('./Results/compare/'+dir+'flowCO2.png')
    fig_W.savefig('./Results/compare/'+dir+'work.png')
    d_all = pd.DataFrame(d_all)
    d_all.to_csv('./Results/compare/'+dir+'data.csv', header=None)


def file_list():
    # experiment A
    file = [
    # low pressure 0.3MPa -> 1.9 MPa
    './2022-12-06/GL840_01_No2_2022-12-06_15-12-36_flow_mod.csv',
    # middle pressure 0.5MPa -> 1.9 MPa
    './2022-12-06/GL840_01_No2_2022-12-06_10-02-01_flow_mod.csv',
    # high pressure 0.8 MPa -> 1.9 MPa
    './2022-12-06/GL840_01_No2_2022-12-06_15-12-36_flow_mod.csv']
    time = [750,5360,2500]
    return file, time


def data_collect(file,time,span):
    res = ResData(file,time,exp_span=span,dir_name='')
    data = {
        'timeline' : res.time_line,
        'pressure' : res.d_loc('Inlet'),
        'water temperature' : res.d_loc('water outlet'),
        "CO2 inlet" : res.d_loc("CO2 inlet"),
        "CO2 outlet" : res.d_loc("CO2 outlet"),
        "CO2 flow"  : res.d_loc("${\dot m_{CO2}}$"),
        "gas in tank" : res.d_loc('gas'),
        "water flow"  : res.d_loc('${\dot v_{water}}$'),
        "Q_water" : res.d_loc('$Q_{water}$'),
        "Q_co2"   : res.d_loc("$Q_{CO2}$"),
        "res"   : res.res_list
        }
    return data


def graph_P(ax,time,data,legend):
    ax.plot(time, data, label = legend, linewidth=2.8)
    ax.set_position([0.12,0.15,0.58,0.8])
    ax.set_xlabel('time [sec]')
    ax.set_ylabel('${P}$ [MPa]')    
    #ax.set_ylim(0,np.max(data)+0.5)
    ax.legend(bbox_to_anchor=(1.03, 1.0), loc='upper left', borderaxespad=0)
    ax.grid(axis='both')
    return ax

def graph_W(ax,time,data,legend):
    ax.plot(time, data, label = legend, linewidth=2.8)
    ax.set_position([0.12,0.15,0.58,0.8])
    ax.set_xlabel('time [sec]')
    ax.set_ylabel('${Q}$ [W]')    
    ax.legend(bbox_to_anchor=(1.03, 1.0), loc='upper left', borderaxespad=0)
    ax.grid(axis='both')
    return ax

def graph_T(ax,time,data,legend):
    ax.plot(time, data, label = legend, linewidth=2.7)
    ax.set_position([0.12,0.15,0.58,0.8])
    ax.set_xlabel('time [sec]')
    ax.set_ylabel('${T}$ [℃]')    
    ax.legend(bbox_to_anchor=(1.03, 1.0), loc='upper left', borderaxespad=0)
    ax.grid()
    return ax

def graph_flow(ax,time,data,legend,label,ylim=[0,0.13]):
    ax.plot(time, data, label = legend, linewidth=2.7)
    ax.set_position([0.14,0.15,0.58,0.8])
    ax.set_xlabel('time [sec]')
    ax.set_ylabel(label)
    ax.legend(bbox_to_anchor=(1.005, 1.0), loc='upper left', borderaxespad=0)
    ax.grid()
    ax.set_ylim(ylim[0],ylim[1])
    return ax


if __name__=='__main__':
    main()
import multiprocessing
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from Eq_balance import balance
import numpy as np
import time
plt.rcParams["font.size"] = 20
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 320

def main():
    start_time = time.time()
    file_list =[
        './exp_data/res/2022-12-26/GL840_01_No3_2022-12-26_15-44-03.csv'
    ]*4

    #start_list = [17100,18000,18950,19600,20300,21300,22300,23100]
    start_list = [18000,18950,19600,20300]
    params = [62.69737115,  0.07982275,  0.41252353]

    params = np.array([10,   0.07855777,  0.38171265])
    params = fitting_param(params, file_list,start_list)
    compare_all(params, file_list,start_list)
    print(' 計算時間　',(time.time() - start_time)/3600,' 時間')


# curve fitting
def fitting_param(corr0, file_, time_):
    bound = ([5, 0.003, 0.1],[100, 0.3,3])
    xscales = [5,0.01, 0.02]
    res_perf = least_squares(lsq_fit_multi,corr0, bounds=bound, verbose = 1,method='trf',x_scale=xscales,args=(file_,time_))
    fitted_params = res_perf.x
    print(fitted_params)
    return fitted_params


# give params and calc simulation results
def param_set(system, corr):
    system.h_CO2 = corr[0]
    system.K1_LDF = corr[1]
    system.k_water = corr[2]

    system.simulation_time = 120
    system.set_dependant_var_IC()
    return system


# run lsq for many cases
def lsq_fit(corr,res_list, time_list):
    diff_all = []
    for i in range(len(res_list)):
        diff_all.extend(diff_sim_exp(corr, res_list[i],time_list[i]))
    return diff_all


# run lsq for many cases with multi process calculation
def lsq_fit_multi(corr,res_list, time_list):
    proc = []
    proc_answer = []
    # prepare the multi process
    for i in range(len(res_list)):
        get_rev,send_rev  = multiprocessing.Pipe(False)
        t = multiprocessing.Process(target=diff_sim_exp,args=(corr, res_list[i],time_list[i],send_rev))
        proc_answer.append(get_rev)
        t.start()
        proc.append(t)
    
    #　計算実行
    for i in range(len(proc)):
        proc[i].join()
    
    #　結果の取得
    diff_all = []
    for i in range(len(proc)):
        diff_all.extend(proc_answer[i].recv())
    return diff_all


# calculate the difference between experiment and simulation with given params
def diff_sim_exp(corr,res_exp,start_sec,send_rev):
    mof_type = 'Uio-66'
    system = balance(mof_type,res_exp,start_sec)

    # set parameters and solve it
    system = param_set(system, corr)
    system.solver()

    # experimental data
    exp_temp_sor = system.res_exp.d_loc('sorbent')
    #exp_temp_gas = system.res_exp.d_loc("CO2 outlet")
    exp_temp_water = system.res_exp.d_loc("water outlet")
    
    # comparison of exp and sim
    dT_sor = []; dT_water = [];dT_gas = []
    for t in range(int(max(system.t))):
        sim_index = np.where(system.t.astype(int)==t)[0][0]
        exp_index = np.where(system.res_exp.time_line==t)[0][0]

        #dT_gas.append(relative_diff(system.gas_T_list[sim_index], exp_temp_gas.iloc[exp_index]))
        dT_sor.append(relative_diff(system.mof_T_list[sim_index], exp_temp_sor.iloc[exp_index]))
        dT_water.append(relative_diff(system.T_HTF_list[sim_index], exp_temp_water.iloc[exp_index]))

    # ratio of importanec of parameters
    #dT_gas = np.array(dT_gas) * 1.0
    #dT_sor = np.array(dT_sor) * 1.0
    dT_water = np.array(dT_water) * 1.0
    # return the difference 
    #diff_sim = np.concatenate((dT_water, dT_sor), axis = 0)
    diff_sim = dT_water
    send_rev.send(diff_sim)
    return diff_sim



# return relative difference
def relative_diff(data_sim, data_exp):
    data_sim = float(data_sim) # K
    data_exp = float(data_exp) + 273.15 # modify celcius to kelvin
    error = abs((data_sim - data_exp))
    return error


# make a graph for each case
def compare_all(params, res_list,time_list):
    for i in range(len(res_list)):
        tot_Q_w = compare_temp(params, res_list[i],time_list[i],i+1)
        print('num : ',i+1,'  ', tot_Q_w, '  J ')


# compare the temperatures from simulation adn experiment
def compare_temp(params, file,start,num):
    mof_type = 'Uio-66'

    System = balance(mof_type,file,start)
    param_set(System, params)
    tot_Q = System.solver()

    fig = plt.figure(figsize=(8,5))
    ax_T = fig.add_subplot(1,1,1)

    # simulation plot
    sim_legends = ['water(sim)']
    sim_temp_data = [System.T_HTF_list-273.15]
    for i in range(len(sim_temp_data)):
        ax_T.plot(System.t, sim_temp_data[i], label = sim_legends[i],linewidth = 2)


    # experimental data
    exp_temp_sor = System.res_exp.d_loc('sorbent')
    #exp_temp_gas = System.res_exp.d_loc("CO2 outlet")
    exp_temp_water = System.res_exp.d_loc("water outlet")
    exp_legends = ['water(exp)']
    exp_temp_data = [exp_temp_water]
    markers = ['s']
    for i in range(len(exp_temp_data)):
        ax_T.scatter(System.res_exp.time_line, exp_temp_data[i], label = exp_legends[i], marker=markers[i], s=15, edgecolors='k')
    ax_T.set_xlabel(' t sec')
    ax_T.set_ylabel(' T ℃')
    ax_T.set_xlim(-1,System.t[-1])
    ax_T.legend()
    ax_T.grid(axis='both')
    fig.savefig('./Fig/temp_compare_'+str(start)+'.png')
    return tot_Q


if __name__=='__main__':
    main()
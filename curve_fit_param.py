
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from Eq_balance import balance
import numpy as np
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 320

def main():
    """
    file_list = [
        './exp_data/res/GL840_01_No2_2022-12-06_15-12-36.csv', # low,750
        './exp_data/res/GL840_01_No2_2022-12-06_10-02-01.csv', # middle,5360,
        './exp_data/res/GL840_01_No2_2022-12-06_15-12-36.csv'  # high, 2500
    ]"""

    file_list =[
        './exp_data/res/GL840_01_No3_2022-12-26_15-44-03.csv'
    ]*8

    start_list = [17100,18000,18950,19600,20300,21300,22300,23100]
    params = np.array([10,  0.32])
    params = fitting_param(params, file_list,start_list)
    compare_all(params, file_list,start_list)


# curve fitting
def fitting_param(corr0, file_, time_):
    bound = ([5, 0.003],[100, 0.5])
    x_scales = ([2, 0.01])
    res_perf = least_squares(lsq_fit,corr0, bounds=bound, verbose = 1, x_scale=x_scales, method='trf',args=(file_,time_))
    fitted_params = res_perf.x
    print(fitted_params)
    return fitted_params


# give params and calc simulation results
def param_set(system, corr):
    system.h_CO2 = corr[0]
    system.K1_LDF = corr[1]
    #system.k_water = corr[2]
    system.set_dependant_var_IC()
    return system


# run lsq for many cases
def lsq_fit(corr,res_list, time_list):
    diff_all = []
    for i in range(len(res_list)):
        diff_all.extend(diff_sim_exp(corr, res_list[i],time_list[i]))
    return diff_all


# calculate the difference between experiment and simulation with given params
def diff_sim_exp(corr,res_exp,start_sec):
    mof_type = 'Uio-66'
    system = balance(mof_type,res_exp,start_sec)

    # set parameters and solve it
    system = param_set(system, corr)
    system.solver()

    # experimental data
    exp_temp_sor = system.res_exp.d_loc('sorbent')
    exp_temp_gas = system.res_exp.d_loc("gas")
    exp_temp_water = system.res_exp.d_loc("water outlet")
    
    # comparison of exp and sim
    dT_gas = []; dT_sor = []; dT_water = []
    for t in range(int(max(system.t))):
        sim_index = np.where(system.t.astype(int)==t)[0][0]
        exp_index = np.where(system.res_exp.time_line==t)[0][0]

        dT_gas.append(relative_diff(system.gas_T_list[sim_index], exp_temp_gas.iloc[exp_index]))
        dT_sor.append(relative_diff(system.mof_T_list[sim_index], exp_temp_sor.iloc[exp_index]))
        dT_water.append(relative_diff(system.T_HTF_list[sim_index], exp_temp_water.iloc[exp_index]))

    # ratio of importanec of parameters
    dT_gas = np.array(dT_gas) * 0.7
    dT_sor = np.array(dT_sor) * 1.0
    dT_water = np.array(dT_water) * 2.0
    # return the difference 
    diff_sim = np.concatenate((dT_gas, dT_sor, dT_water), axis = 0)
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
    scp, tot_Q, tot_H_ads = System.solver()

    fig = plt.figure(figsize=(8,5))
    ax_T = fig.add_subplot(1,1,1)

    # simulation plot
    sim_legends = ['gas (sim)', 'sor(sim)', 'water(sim)']
    sim_temp_data = [System.gas_T_list-273.15, System.mof_T_list-273.15, System.T_HTF_list-273.15]
    for i in range(len(sim_temp_data)):
        ax_T.plot(System.t, sim_temp_data[i], label = sim_legends[i],linewidth = 2)


    # experimental data
    exp_temp_sor = System.res_exp.d_loc('sorbent')
    exp_temp_gas = System.res_exp.d_loc("gas")
    exp_temp_water = System.res_exp.d_loc("water outlet")
    exp_legends = ['gas (exp)', 'sor(exp)', 'water(exp)']
    exp_temp_data = [exp_temp_gas, exp_temp_sor, exp_temp_water]
    for i in range(len(exp_temp_data)):
        ax_T.scatter(System.res_exp.time_line, exp_temp_data[i], label = exp_legends[i], s=15, edgecolors='k')
    ax_T.set_xlabel(' t sec')
    ax_T.set_ylabel(' T ℃')
    ax_T.set_xlim(-1,System.t[-1])
    ax_T.legend()
    ax_T.grid(axis='both')
    fig.savefig('./Fig/temp_compare_'+str(num)+'.png')
    return tot_Q


if __name__=='__main__':
    main()
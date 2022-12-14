import pandas as pd
from time import time
import numpy as np
from HC_AHP_System import MOF
from Eq_balance import balance 
from os import mkdir, path

if not path.exists('./Results'):
    mkdir('./Results')

def main():
    mof_type = 'MIL-101'
    #mof_type = 'Uio-66'
    res_file = './Results/303_2_30.csv'
    multi_solver(mof_type, res_file)


# solve the system for several cases
def multi_solver(name, file):
    # length of the simulation time
    simulation_time_list = np.linspace(200, 1800,num=50) # sec
    # lower side temperature
    #Tl_list = np.linspace(278.15, 313.15, num=11)
    # higher side pressure
    #Ph_list = np.linspace(21e5, 51e5, num=11)
    # lower side temperature
    #Pl_list = np.linspace(1e5, 12e5, num=13)

    """one temp pres
    """
    Tl_list = np.linspace(303.15, 303.15, num=1)
    Ph_list = np.linspace(30e5, 30e5, num=1)
    Pl_list = np.linspace(2e5, 2e5, num=1)

    # set result format
    df_index = pd.DataFrame(['simulation_time','gas_T','gas_P_init','gas_P','SCP','COP','totQ']).T
    df_index.to_csv(file, header=False, index=False)

    # calc for many cases one by one
    for t in simulation_time_list:
        for Tl in Tl_list:
            for Pl in Pl_list:
                for Ph in Ph_list:
                    params = {
                        'simulation_time' : t,
                        'gas_T'           : Tl,
                        'gas_P_init'      : Pl,
                        'gas_P'           : Ph,   
                        }
                    res = one_solver(name,params)
                    params.update(res)
                    df_res = pd.DataFrame(params.values()).T
                    df_res.to_csv(file, mode='a',index=False, header=False)


# solve the system for one case
def one_solver(name,param_dict):
    # initialize the system
    case = balance(name)

    # set variables 
    print(param_dict)
    case.simulation_time = param_dict['simulation_time']
    case.gas_T = param_dict['gas_T']
    case.gas_P_init = param_dict['gas_P_init']
    case.gas_P = param_dict['gas_P']

    # reset the variables
    case.set_dependant_var_IC()

    # solve the system
    scp, cop, totQ = case.solver()
    res_dict = {
        'SCP' : scp,
        'COP' : cop,
        'totQ': totQ
    }

    return res_dict


if __name__=='__main__':
    main()


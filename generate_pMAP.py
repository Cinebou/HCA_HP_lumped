import pandas as pd
from time import time
import matplotlib.pyplot as plt
import numpy as np
from HC_AHP_System import MOF
from Eq_balance import balance 
from os import mkdir, path

if not path.exists('./Results'):
    mkdir('./Results')

def main():
    start_time = time()
    mof_type = 'MIL-101'
    #mof_type = 'Uio-66'
    time_len_solver(mof_type)


# solve the system for several cases
def time_len_solver(name):
    simulation_time_list = np.linspace(500, 1000,50) # sec
    for t in simulation_time_list:
        params = {
            'simulation_time' : t
        }
        res = one_solver(name,params)
        params.update(res)
        df_res = pd.DataFrame(params.values()).T
        df_res.to_csv('./Results/time_len_solver.csv', mode='a',index=False, header=False)


# solve the system for one case
def one_solver(name,param_dict):
    # initialize the system
    case = balance(name)

    # set variables 
    #case.gas_P_init = param_dict['gas_P_init']
    #case.gas_P = param_dict['gas_P']
    case.simulation_time = param_dict['simulation_time']

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

g= {
  "brand": "Ford",
  "model": "Mustang",
  "year": 1964
}
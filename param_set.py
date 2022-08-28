import sys
sys.path.append('../')
from fluidProp import VLEFluid

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
from pprint import pprint


def main():
    #param_sat_pressure()
    param_heat_ad()


def param_heat_ad():
    data = pd.read_csv('./Parameter_fix/heat.csv')
    amount, heat = data["adsorption amount"], data["experiment"]*1000
    popt, pcov = curve_fit(func_exp, amount, heat, p0=[17,1,-1],  maxfev=5000) # poptは最適推定値、pcovは共分散
    
    amount_2 = np.linspace(0, 30, 30)
    estimated_heat = func_exp(amount_2,popt[0],popt[1],popt[2])
    fig_sat = plt.figure()
    ax = fig_sat.add_subplot(111)
    ax.plot(amount, heat, label = 'refprop')
    ax.plot(amount_2, estimated_heat, label='curve fit')
    ax.set_xlabel('adsorption amount mol/kg')
    ax.set_ylabel('heat J/mol')
    ax.legend()
    print('Heat of adsorption curve is \n H = a + b * exp(c * M)'.format(popt[0],popt[1],popt[2]))
    plt.savefig('./Fig/heat of adsorpion.png')
    plt.show()


def param_sat_pressure():
    Temp, pres = sat_curve()
    popt, pcov = curve_fit(func3,Temp,pres) # poptは最適推定値、pcovは共分散

    Temp_2 = np.linspace(220, 400, 180)
    estimated_Psat = func3(Temp_2,popt[0],popt[1],popt[2], popt[3])
    fig_sat = plt.figure()
    ax = fig_sat.add_subplot(111)
    ax.plot(Temp, pres, label = 'refprop')
    ax.plot(Temp_2,estimated_Psat, label='curve fit')
    ax.set_xlabel('temperature K')
    ax.set_ylabel('pressure Pa')
    ax.legend()
    print('saturation curve of CO2 is \n P = {} + {} * T + {} * T^2 + {} * T^3'.format(popt[0],popt[1],popt[2], popt[3]))
    plt.savefig('./Fig/satureation_pressure.png')
    plt.show()
    return 0

# saturation vapor-liquid curve of the CO2
def sat_curve():
    gas = VLEFluid('CO2')
    p_sat_list = []
    temp_list = np.linspace(220, 303, 200)
    for temp in temp_list:
        p_sat = gas.calc_VLE_T(temp).p_v
        p_sat_list.append(p_sat)
    return np.array(temp_list), np.array(p_sat_list)


def func3(X, a, b, c, d): # 3次式近似
    Y = a + b * X + c * X ** 2 + d * X ** 3
    return Y

def func_exp(X, a, b, c): # 3次式近似
    Y = a + b * np.exp(c * X)
    return Y

if __name__ == "__main__": 
    main()
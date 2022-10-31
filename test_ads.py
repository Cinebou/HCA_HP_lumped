from HC_AHP_System import MOF
from fluidProp import VLEFluid
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from pprint import pprint
from math import exp, log 
from HC_AHP_System import MOF

def main():
    cp_calc()
    #isotherm_example()
    #relaxation_curve()
    #test_refprop()
    #sat_curve()
    #DubininAstakhov()


# calculation of heat capacity of mil101
def cp_calc():
    mof = MOF('MIL-101')
    temp_line = np.linspace(200,400,200)
    cp_line = [mof.heat_capacity_mof(t) for t in temp_line]
    plt.plot(temp_line, cp_line)
    plt.show()

# saturation vapor-liquid curve of the CO2
def sat_curve():
    gas = VLEFluid('CO2')
    p_sat_list = []
    temp_list = np.linspace(220, 303, 200)
    print(temp_list)
    for temp in temp_list:
        p_sat = gas.calc_VLE_T(temp).p_v
        p_sat_list.append(p_sat)

    fig_sat = plt.figure()
    ax = fig_sat.add_subplot(111)
    ax.plot(temp_list, p_sat_list)
    ax.set_xlabel('temperature K')
    ax.set_ylabel('pressure Pa')
    plt.savefig('./Fig/satureation_pressure.png')
    plt.show()

def DubininAstakhov():
    gas = VLEFluid('CO2')
    p_list = np.linspace(1, 1000000, 20)
    
    # DA parameters
    qo = 1.08 * 1000/44 # mol(CO2)/kg(MOF)
    E = 4400
    n = 1.1
    temp1 = 273.15
    ad_amount1 = [qo*exp(-(Ad_pot(p, temp1)/E)**n) for p in p_list]
    plt.plot(p_list,ad_amount1)

    temp2 = 298.15
    ad_amount2 = [qo*exp(-(Ad_pot(p, temp2)/E)**n) for p in p_list]
    plt.plot(p_list,ad_amount2)

    temp3 = 323.15
    ad_amount3 = [qo*exp(-(Ad_pot(p, temp3)/E)**n) for p in p_list]
    plt.plot(p_list,ad_amount3)
    plt.show()

def Ad_pot(pressure, temp):
    gas = VLEFluid('CO2')
    mof = MOF('MIL-101')
    sat_p = mof.sat_pressure(temp)
    A = mof.R * temp * log(sat_p / pressure)
    return A

# calculate isotherm 
def isotherm_example():
    test_tank = MOF("MIL-101")
    T = 273.15
    p_po = np.linspace(0.1,1,100)

    amount = []
    for p in p_po:
        p_Pa = p * test_tank.sat_pressure(T)
        amount.append(test_tank.eq_loading(p_Pa, T))

    fig_iso = plt.figure()
    ax_iso = fig_iso.add_subplot(1,1,1)
    ax_iso.plot(p_po * test_tank.sat_pressure(T)/100000, amount)
    ax_iso.set_xlabel('p bar')
    ax_iso.set_ylabel('M mol/kg')
    plt.show()
    plt.savefig('./Fig/isotherm.jpg')
    

# calculate reluxation curve
def relaxation_curve():
    test_tank = MOF("MIL101")
    init_P = test_tank.gas_P_init

    init_amount = test_tank.eq_loading(init_P, test_tank.gas_T)
    test_tank.loading = init_amount
    print('initial adsorption amount      :  ', init_amount)
    print('equilibrium adsorption amount  :  ', test_tank.eq_loading(test_tank.gas_P, test_tank.gas_T))

    max_t = 500 # sec
    t = np.linspace(0, max_t, max_t)
    sol = odeint(m, init_amount, t, args=(test_tank,))

    fig_rel = plt.figure()
    ax_rel = fig_rel.add_subplot(1,1,1)
    ax_rel.plot(t, sol[:,0])
    ax_rel.set_xlabel('t sec')
    ax_rel.set_ylabel('M mol/kg')
    plt.savefig('./Fig/relaxation.jpg')


def m(var, t, tank):
    tank.loading = var
    dmdt = tank.ads_speed_LDF() 
    return dmdt

def test_refprop():
    gas = VLEFluid('CO2')
    T = 283.15
    P = 2500000
    h = gas.calc_fluidProp_pT(P, T).h 
    print("enthalpy  J/kg : ",h)
    print("enthalpy  J/mol: ",h* gas.get_molar_mass())

    print('molar mass of CO2', gas.get_molar_mass())

    print("cp of the gas  : ", gas.calc_VLE_liquid_T(T).cp)
    print("volume of the gas : ", gas.calc_fluidProp_pT(P,T).vol)
    return 0

if __name__ == "__main__": 
    main()
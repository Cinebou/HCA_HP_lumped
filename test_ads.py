from HC_AHP_System import MOF
from fluidProp import VLEFluid
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint
from pprint import pprint
from math import exp, log 
from multiprocessing import Pool
plt.rcParams["font.size"] = 24
plt.rcParams["font.family"] = "Times New Roman"

def main():
    #mass_flow()
    #heat_tot()
    #cp_calc()
    #isotherm_example()
    #relaxation_curve()
    #test_refprop()
    sat_curve()
    #DubininAstakhov()
    #entropy_example()
    #multi_process()
    #latent_heat()
    #accum_heat_example()

def mass_flow():
    data_mass = pd.read_csv('./Log/log_mass.csv',names=["m_in", "m_ads","m_out", "time"])
    plt.plot(data_mass["time"], data_mass["m_in"] ,label="m_in")
    #plt.plot(data_mass["time"], data_mass["m_in"] / data_mass["m_out"] ,label="in out rratio")
    plt.plot(data_mass["time"], data_mass["m_out"],label="m_out")
    plt.xlabel("time [sec]")
    plt.ylabel("$\dot m_{CO2}$ [mol/sec]")
    plt.legend()
    plt.show()

def heat_tot():
    test_tank = MOF("Uio-66",res_file='',start_time=0)
    Ts = 21.6
    Tf = 20.4
    ps =  0.7802e06
    pf =  1.935809090e06
    p_po = np.linspace(0.001,1,100)
    amount_s = []
    amount_f = []
    for p in p_po:
        p_Pa_s = p * test_tank.sat_pressure(Ts + 273.15)
        amount_s.append(test_tank.eq_loading(p_Pa_s, Ts + 273.15))
        p_Pa_f = p * test_tank.sat_pressure(Tf + 273.15)
        amount_f.append(test_tank.eq_loading(p_Pa_f, Tf + 273.15))

    fig_iso = plt.figure(figsize=(7.0, 4.5))
    ax_iso = fig_iso.add_subplot(1,1,1)
    ax_iso.plot(p_po * test_tank.sat_pressure(Ts + 273.15)/1000000, amount_s, label = str(Ts) + ' ℃')
    ax_iso.plot(p_po * test_tank.sat_pressure(Tf + 273.15)/1000000, amount_f, label = str(Tf) + ' ℃')
    
    
    ax_iso.scatter(ps/1e06, test_tank.eq_loading(ps, Ts + 273.15), s=60)
    ax_iso.scatter(pf/1e06, test_tank.eq_loading(pf, Tf + 273.15), s=60)
    """
    ps = 0.7802e06
    ax_iso.scatter(ps/1e06, test_tank.eq_loading(ps, Ts + 273.15), s=60,label='0.8 MPa')
    ps = 1.961e06
    ax_iso.scatter(ps/1e06, test_tank.eq_loading(ps, Ts + 273.15), s=60,c='b')"""
    dm = test_tank.eq_loading(pf, Tf + 273.15)-test_tank.eq_loading(ps, Ts + 273.15)
    print(' dm = ',dm,' mol/kg')
    print('H ads. = ', dm * 0.1349* 24.8994 * 1000, ' J')
    ax_iso.set_xlabel('p MPa')
    ax_iso.set_ylabel('M mol/kg')
    ax_iso.set_xlim(0, 2.7)
    ax_iso.set_ylim(0,9)
    ax_iso.set_position([0.1,0.16, 0.85,0.8])
    ax_iso.legend()
    ax_iso.grid()
    plt.savefig('./Fig/isotherm.jpg')
    


def latent_heat():
    gas = VLEFluid('CO2')
    T_list = np.linspace(263, 400,num=20000)
    h_fg_list = np.array([gas.calc_fg_latent(T) for T in T_list])
    plt.scatter(T_list, h_fg_list/1000)
    plt.ylim(0 ,15)
    plt.show()

def sample_func(initial_num):
    for i in range(100):
        initial_num += 1

def multi_process():
    initial_num_list = list(range(1000000))
    with Pool(processes=4) as p:
        p.map(func=sample_func, iterable=initial_num_list)

def test():
    mof = MOF('MIL-101')
    print(mof.eq_loading(300000,303.15))

# entropy reading
def entropy_example():
    gas = VLEFluid('CO2')
    #h = gas.calc_fluidProp_ps(1500000, 2356.7).h
    mof = MOF('MIL-101')
    t,h = mof.calc_gas_in(303.15, 300000, 4000000)
    print(t,'\n',h)

# calculation of heat capacity of mil101
def cp_calc():
    mof = MOF('MIL-101')
    temp_line = np.linspace(200,400,200)
    cp_line = [mof.heat_capacity_mof(t) for t in temp_line]

    print(mof.heat_capacity_mof(298.15))
    plt.plot(temp_line, cp_line)
    #plt.xlim(240, 440)
    #plt.ylim(200, 800)
    plt.ylabel('Cp [J/kg/K]')
    plt.xlabel('temperature [K]')
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
    test_tank = MOF("Uio-66",res_file='',start_time=0)
    T = 303
    p_po = np.linspace(0.001,0.5,100)
    amount_a = []
    for p in p_po:
        p_Pa = p * test_tank.sat_pressure(T)
        amount_a.append(test_tank.eq_loading(p_Pa, T))

    fig_iso = plt.figure(figsize=(9,6))
    ax_iso = fig_iso.add_subplot(1,1,1)
    ax_iso.plot(p_po * test_tank.sat_pressure(T)/1000000, amount_a, label='a')
    ax_iso.set_xlabel('${p}$ MPa')
    ax_iso.set_ylabel('${M}$ mmol/g')
    ax_iso.set_xlim(0,3.2)
    ax_iso.set_position([0.1,0.15,0.8,0.8])
    ax_iso.grid()
    plt.savefig('./Fig/isotherm.jpg')

# calculate isotherm 
def accum_heat_example():
    test_tank = MOF("Uio-66",res_file='',start_time=0)
    T = 303.15
    p_po = np.linspace(0.001,0.5,100)
    amount_a = []
    for p in p_po:
        p_Pa = p * test_tank.sat_pressure(T)
        amount_a.append(test_tank.eq_loading(p_Pa, T))

    fig_iso = plt.figure()
    ax_iso = fig_iso.add_subplot(1,1,1)
    p = p_po * test_tank.sat_pressure(T)/1000000
    h =  np.array(amount_a) * 24.8994
    ax_iso.plot(p,h, label='a')
    data = pd.DataFrame([p,h])
    data.to_csv('./Results/accum_heat.csv')
    ax_iso.set_xlabel('p MPa')
    ax_iso.set_ylabel('M kJ/kg')
    #plt.show()
    plt.savefig('./Fig/accum_heat.jpg')

# calculate reluxation curve
def relaxation_curve():
    test_tank = MOF("MIL101",res_file='',start_time=0)
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
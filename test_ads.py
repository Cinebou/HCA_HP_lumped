from HC_AHP_System import MOF
from fluidProp import VLEFluid
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from pprint import pprint

def main():
    #isotherm_example()
    #relaxation_curve()
    #plt.show()
    test_refprop()



# calculate isotherm 
def isotherm_example():
    test_tank = MOF("MIL101")
    T = 273.15
    p_po = np.linspace(0,1,100)

    amount = []
    for p in p_po:
        p_Pa = p * test_tank.sat_pressure(T)
        amount.append(test_tank.eq_loading(p_Pa, T))

    fig_iso = plt.figure()
    ax_iso = fig_iso.add_subplot(1,1,1)
    ax_iso.plot(p_po * test_tank.sat_pressure(T)/100000, amount)
    ax_iso.set_xlabel('p bar')
    ax_iso.set_ylabel('M mol/kg')
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
from HC_AHP_System import MOF
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from pprint import pprint

def main():
    moisture_content_Xeq()
    #curve_sat()


# saturation pressure
def curve_sat():
    mof = MOF('TMPS-1.5')
    Td = 300
    #print('saturation pressure at ', Td, ' K is ',mof.sat_pressure(Td) )
    Tlin = np.linspace(273, 350, 100)
    Psat = []
    for T in Tlin:
        Psat.append(mof.sat_pressure(T))

    plt.plot(Tlin, Psat)
    plt.show()


# moisture content of the TMPS-15
def moisture_content():
    mof = MOF('TMPS-1.5')
    Plin = np.linspace(0, 0.85, 100)

    # adsorption
    Td = 0 + 273.15
    Cw0_ad = [mof.isotherm(Td, p, hys='ad') for p in Plin*mof.sat_pressure(Td)]
    Td = 25 + 273.15
    Cw25_ad = [mof.isotherm(Td, p, hys='ad') for p in Plin*mof.sat_pressure(Td)]
    Td = 50 + 273.15
    Cw50_ad = [mof.isotherm(Td, p, hys='ad') for p in Plin*mof.sat_pressure(Td)]

    # desorption
    Td = 0 + 273.15
    Cw0_de = [mof.isotherm(Td, p, hys='de') for p in Plin*mof.sat_pressure(Td)]
    Td = 25 + 273.15
    Cw25_de = [mof.isotherm(Td, p, hys='de') for p in Plin*mof.sat_pressure(Td)]
    Td = 50 + 273.15
    Cw50_de = [mof.isotherm(Td, p, hys='de') for p in Plin*mof.sat_pressure(Td)]

    # plot graph
    fig_iso = plt.figure()
    ax = fig_iso.add_subplot(1,1,1)
    ax.plot(Cw0_ad,  Plin*mof.sat_pressure(0 + 273.15), label = "ad 0℃")
    ax.plot(Cw25_ad, Plin*mof.sat_pressure(25 + 273.15), label = "ad 25℃")
    ax.plot(Cw50_ad, Plin*mof.sat_pressure(50 + 273.15), label = "ad 50℃")
    ax.plot(Cw0_de,  Plin*mof.sat_pressure(0 + 273.15), label = "de 0℃")
    ax.plot(Cw25_de, Plin*mof.sat_pressure(25 + 273.15),  label = "de 25℃")
    ax.plot(Cw50_de, Plin*mof.sat_pressure(50 + 273.15),  label = "de 50℃")
    ax.legend()
    ax.set_xlabel('pressure Pa')
    ax.set_ylabel('Cw kg/kg')
    plt.savefig('./Fig/adsorption amount.png')
    plt.show()


# moisture content of the TMPS-15
def moisture_content_Xeq():
    mof = MOF('TMPS-1.5')
    Xeq = np.linspace(0, 5e-3, 100)
    Td = 0 + 273.15
    # adsorption
    Cw0_ad = [mof.isotherm_Xeq(Td, x, hys='ad') for x in Xeq]

    # desorption
    Cw0_de = [mof.isotherm_Xeq(Td, x, hys='de') for x in Xeq]

    # plot graph
    fig_iso = plt.figure()
    ax = fig_iso.add_subplot(1,1,1)
    ax.plot(Cw0_ad,  Xeq * 1000, label = "ad " + str(Td - 273.15))
    ax.plot(Cw0_de,  Xeq * 1000, label = "de " + str(Td - 273.15))

    ax.legend()
    ax.set_ylabel('Xeq')
    ax.set_xlabel('Cw kg/kg')
    plt.savefig('./Fig/adsorption amount_Xeq.png')
    plt.show()

if __name__ == "__main__": 
    main()
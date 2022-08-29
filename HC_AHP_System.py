"""
Author Hibiki Kimura,
2022/08/23, lumped parameter batch mode MOF heat pump model

"""
from math import sqrt, exp, log
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint

# unit
# energy, J
# mass, kg
# temperature, K
# volume, m3
# time, second
# pressure, Pa

"""
Refprop version 10 (not version 9) is necessary to read thermal properties of CO2
"""

class MOF():
    def __init__(self, name):
        self.name = name

    # water saturation pressure, IAPWS-IF97
    # http://twt.mpei.ac.ru/mcs/worksheets/iapws/IAPWS-IF97-Region4.xmcd    
    def sat_pressure(self, Td):
        n1 = 1167.0521452767; n2 = -724213.16703206; n3 = -17.073846940092; n4 = 12020.82470247; n5 =-3232555.0322333
        n6 = 14.91510861353; n7 = -4823.2657361591; n8 = 405113.40542057; n9 = -0.23855557567849; n10 = 650.17534844798

        theta = Td + n9 / (Td - n10)
        A = theta**2 + n1 * theta + n2
        B = n3*theta**2 + n4 * theta + n5
        C = n6*theta**2 + n7*theta + n8
        Ps = (2*C / (-B + sqrt(B**2-4*A*C)))**4
        return Ps * 1000000  # Pa


    # equilibrium adsorption / desorption equation
    def isotherm(self, Td, P, hys = 'ad'):
        # relative pressure
        Ps = self.sat_pressure(Td)
        Rp = P/Ps
        
        # parameters 
        n0, n1, n2, n3 = self.read_isotherm(Td, Rp, hys)
        # Moisure content, kg_water/kg_desiccant
        Cw = n0 + n1 * Rp + n2 * Rp**2 + n3 * Rp**3
        return Cw

 
    # equilibrium adsorption / desorption equation from the air humidity ratio
    def isotherm_Xeq(self, Td, Xeq, hys = 'ad'):
        # relative pressure
        Ps = self.sat_pressure(Td)
        Rp = Xeq* 101325.0/(Xeq+1)/Ps
        print('Pa  : ', Rp * Ps)
        
        # parameters 
        n0, n1, n2, n3 = self.read_isotherm(Td, Rp, hys)
        # Moisure content, kg_water/kg_desiccant
        Cw = n0 + n1 * Rp + n2 * Rp**2 + n3 * Rp**3
        return Cw

    # calculate the slope of the isotherm
    def isotherm_diff(self, Td, Xeq, hys = 'ad'):
        # relative pressure
        Ps = self.sat_pressure(Td)
        Rp = Xeq* 101325.0/(Xeq+1.0)/Ps
        print('Pa  : ', Rp * Ps)
        
        # parameters 
        n0, n1, n2, n3 = self.read_isotherm(Td, Rp, hys)
        # Moisure content, kg_water/kg_desiccant
        dpdx = 1/(Xeq+1)**2 * 101325.0
        dCwdX = (n1 + 2 * n2 * Rp + 3 * n3 * Rp**2) * dpdx
        return dCwdX

    
    # read parameters from file
    def read_isotherm(self, Td, Rp, hys):
        params = pd.read_csv('./parameter_fix/isotherm_'+hys+'.csv')
        # temperature
        if -12.5 <= Td -273.15 < 12.5:
            Trep = 0
        elif 12.5 <= Td - 273.15 < 37.5:
            Trep = 25
        elif 37.5 <= Td -273.15 <= 62.5:
            Trep = 50
        else:
            assert (Td >= -12.5 + 273.15), "Temperature is out of range in the isotherm, too cold"
            assert (Td <= 62.5 + 273.15), "Temperature is out of range in the isotherm, too hot"
        
        # pressure
        Rp_range = pd.Series(params[params['temp']==Trep]['p/po'])
        if 0 <= Rp < Rp_range.iloc[0]:
            Prep = Rp_range.iloc[0]
        elif Rp_range.iloc[0] <= Rp < Rp_range.iloc[1]:
            Prep = Rp_range.iloc[1]
        elif Rp_range.iloc[1] <= Rp <= Rp_range.iloc[2]:
            Prep = Rp_range.iloc[2]
        else:
            assert (Rp >= 0), "partial pressure is out of range in the isotherm, too low presssure"
            assert (Rp <= Rp_range.iloc[2]), "partial pressure is out of range in the isotherm, too high pressure "

        # take out parameter
        n0 = params[(params['temp']==Trep) & (params['p/po']==Prep)]['n0']
        n1 = params[(params['temp']==Trep) & (params['p/po']==Prep)]['n1']
        n2 = params[(params['temp']==Trep) & (params['p/po']==Prep)]['n2']
        n3 = params[(params['temp']==Trep) & (params['p/po']==Prep)]['n3']
        return float(n0), float(n1), float(n2), float(n3)

    
    # calculate the coefficient in the LDF equation
    def LDF_const(self, Td, p, key = 'ad'):
        # gas phase diffusion
        rp = 1.8e-9 / 2
        D_ao = 1.735e-9 * Td**1.685 / p
        D_ak = 97 * rp * sqrt(Td/18.015) 
        egv = 0.3
        tgv = 3.0
        Dp = 1/(1/D_ao + 1/D_ak*tgv/egv)

        # surface diffusion
        # https://www.sciencedirect.com/science/article/pii/S0038092X97000807
        D_o = 1.6e-6
        ks = -0.974e-3
        Qst = 2257 # kJ/kg, approximatry same as vaporization heat
        D_s = D_o * exp(ks * Qst/Td)

        # diffusivity in the desiccant felt
        Xeq = p/(101325-p)
        dCdX = self.isotherm_diff(Td, Xeq, hys = key)

        rho_d = 281 # kg/m3
        rho_air = 1.29 # kg/m3
        Dd = Dp + D_s*rho_d/rho_air*dCdX

        # the average mass transfer coefficient inside the desiccantfelt (kd)
        thickness = 0.15e-3
        kd = Dd*rho_air/thickness
        return kd

def test():
    mof = MOF('TMPS-1.5')
    
if __name__ == '__main__':
    test()
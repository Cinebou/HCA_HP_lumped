# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 09:21:41 2020

@author: gibelhaus
"""
from pylab import *
import Fluid_RP10_2 as fl

class VLEFluid:
    """Fluid
    
    Class containing the states of the fluid and 
    corresponding functions
    """
    def __init__(self,name, h=None, h_l=None,h_v=None,p_v=None, cp = None, vol = None):
        self.fluid = name
        fl.init_RP(self.fluid)
        self.h = h
        self.h_l = h_l
        self.h_v = h_v
        self.p_v = p_v
        self.vol = vol
        self.cp = cp


    """
    def __init__(self,name,computeVLEAdditionalProperties=False):
        # Constructor
        Instantiate the fluid
        # name: String defining name of the working pair
        self.fluidProp = prop.VLEFluid(name,computeVLEAdditionalProperties=computeVLEAdditionalProperties)
    """
    # the properties of saturated vapor calculated from temperature
    def calc_VLE_T(self,T,unit_T='K'):
        """Calculates vapor liquid equilibrium
            
        T:      Array_like temperature in K
        unit_T: String defining unit of temperature input 'K' or 'C'
        
        returns: NumPy array   
        """

        if unit_T=='K':
            T = T
        elif unit_T=='C':
            T = T + 273.15
        else:
            pass

        q = 1 #fluid quality is saturated vapor
        prop_list = ["T",'p','v','u','h','s','BETA','q']
        VLE = fl.zs(["T","q"],[T,q],prop_list,self.fluid,Eh="default")
        VLEFluid.h_v = VLE[4]
        VLEFluid.p_v = VLE[1]
        return VLEFluid

    # the properties of saturated liquid calculated from temperature
    def calc_VLE_liquid_T(self,T,unit_T='K'):
        """Calculate the enthalpy of saturated liquid"""
        if unit_T=='K':
            T = T
        elif unit_T=='C':
            T = T + 273.15
        else:
            pass

        q = 0 #fluid quality is saturated liquid
        prop_list = ["T",'p','cp','q']
        VLE = fl.zs(["T","q"],[T,q],prop_list,self.fluid,Eh="default")
        VLEFluid.cp = VLE[2]  # J/kg/K
        return VLEFluid


    # the parameter of saturated vapor calcurated from pressure 
    def calc_VLE_p(self,p):
        """Calculates vapor liquid equilibrium
            
        p:      Array_like pressure in Pa
        
        returns: NumPy array   
        """
        q = 1 #fluid quality is saturated vapor
        prop_list = ["T",'p','v','u','h','s','q']
        VLE = fl.zs(["p","q"],[p, q],prop_list,self.fluid,Eh="default")
        VLEFluid.h_v = VLE[4]
        return VLEFluid

    # the parameter of unsaturated vapor calculated from temperature and pressure
    def calc_fluidProp_pT(self,p,T):
        """Calculated fluid properties
        """
        VLE = fl.zs(["T", "p"],[T, p],["T",'p','v','h','cp','q'],self.fluid,Eh="default")
        VLEFluid.h = VLE[3]
        VLEFluid.vol = VLE[2]
        VLEFluid.cp = VLE[4]
        return VLEFluid

    # search molar mass of the fluid
    def get_molar_mass(self):
        return fl.getInfo(self.fluid)[0]
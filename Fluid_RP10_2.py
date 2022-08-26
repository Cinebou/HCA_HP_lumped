# -*- coding: utf-8 -*-

import os; os.environ['RPPREFIX'] = r'c:/Program Files (x86)/REFPROP'
from numpy import *
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import time




Gr0=["T",'p','v','u','h','s','BETA','q','dvis','kvis','lam','Pr',"cp0","cp","M"]

dTk=273.15
fluid_old="default"

################## Einheitensystem ###################################

#default : K, Pa, J, W, m^3, kg
#CBar    : C, bar, kJ, W, m^3, kg
#CKPa    : C, kPa, kJ, W, m^3, kg
#mol     : K,Pa, J, W, m^3, mol


def init_RP(name,Eh="default"):
    global M,Units,qflag, RP, Out0, z0, iMass
   
  
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    RP.SETPATHdll(os.environ['RPPREFIX'])
    #RP.SETFLUIDSdll(name)
    
    if ".MIX" in name:
        z0=RP.SETMIXTUREdll(name)[0]  # mixed fluid
    else:
        RP.SETFLUIDSdll(name)  # pure fluid
        z0=[1.0]

    if Eh=="mol":
        Units = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
        Out0="T;P;V;E;H;S;BETA"
        qflag=0
        iMass = 0
    else:
        Units = RP.GETENUMdll(0, "MASS BASE SI").iEnum
        Out0="T;P;V;E;H;S;BETA"
        qflag=0
        iMass = 0
    M=RP.REFPROPdll("","","MM",Units,0,0,0,0,[1.0]).Output[0] #Molmasse kg/mol
    return






def zs(Gr,In,Out,name,Eh="default"):
    
    global fluid_old  # at first, fluid_old == "default"
    """
    Gr: List containing two strings of symbols of the state that shall be inserted, e.g., ["T","s"] -> Input: temperature and spec. entropy
    In: List of values of the state variables defined in Gr
    Out: List containing strings of symbols of the state variables that shall be computed, e.g., ["T","lam"] -> Output: temperature and heat conductivity
    Name: String of the fluid name as defined in the list below
    Eh: String concerning the uni system
    
    Supported input combinations of Gr
         ["T","p"]  temperature, pressure   
         ["T","q"]  temperature, steam quality  (q=0:saturated liquid quality, q=1:saturated vapor quality)
         ["T","v"]  temperature, spec. volume
         ["p","v"]  pressure, spec. volume
         ["p","q"]  pressure, steam quality
         ["p","h"]  pressure, spec. enthalpy
         ["p","s"]  pressure, spec. entropy
    
    Supported outputs (Out) 
        T    temperature                        
        p    pressure                             
        v    spec. volume                        
        u    spec. internal energy                
        h    spec. enthalpy                      
        s    spec. entrop                         
        q    steam quality                       
        dvis dynamic viscosity                   
        kvis kinematic viscosity                   
        lam  heat conductivity                     
        Pr   Prandtl-number
        cp0  isobaric ideal gas heat capacity     
        cp   isobaric heat capacity                
        M    molar mass            
        BETA   isothermal expansion                
        
    Units for in- and output, defined by Eh
    
       Eh= "default"  "CBar"   "CKPa"   "mol"
                      ¦        ¦        ¦
        T     K        C        C        K
        p     Pa       bar      kPa      Pa
        v     m3/kg    m3/kg    m3/kg    m3/mol
        u     J/kg     kJ/kg    kJ/kg    J/mol
        h     J/kg     kJ/kg    kJ/kg    J/mol
        s     J/kg/K   kJ/kg/K  kJ/kg/K  J/mol/K
        q     kg/kg    kg/kg    kg/kg    mol/mol
        dvis  Pa s     Pa s     Pa s     Pa s
        kvis  m2/s     m2/s     m2/s     m2/s
        lam   W/m/K    kW/m/K   kW/m/K   W/m/K
        Pr   
        cp0   J/kg/K   kJ/kg/K  kJ/kg/K  J/mol/K
        cp    J/kg/K   kJ/kg/K  kJ/kg/K  J/mol/K
        M     kg/mol   kg/mol   kg/mol   kg/mol
        BETA    /K       /K       /K       /K
    """
    
    if name!=fluid_old:
        init_RP(name,Eh)   # setting of Out0, iUnits and fluid_name
        fluid_old=name 
    
    
    z=z0
    
    # T und p

    if Gr[0]=="T" and Gr[1]=="p":
        if Eh=="CBar":
            T=In[0]+dTk
            p=In[1]*1e5
        elif Eh=="CKPa":
            T=In[0]+dTk
            p=In[1]*1e3
        else:
            T=In[0]
            p=In[1]
        X=RP.REFPROPdll("", "TP", Out0, Units,iMass, qflag, T, p,z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # p und T
    if Gr[0]=="p" and Gr[1]=="T":
        if Eh=="CBar":
            T=In[1]+dTk
            p=In[0]*1e5
        elif Eh=="CKPa":
            T=In[1]+dTk
            p=In[0]*1e3
        else:
            T=In[1]
            p=In[0]
        X=RP.REFPROPdll("", "TP", Out0, Units, iMass, qflag, T, p, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
 ###########################################################       
    # T and q
    if Gr[0]=="T" and Gr[1]=="q":
        if Eh=="CBar":
            T=In[0]+dTk
            q=In[1]
        elif Eh=="CKPa":
            T=In[0]+dTk
            q=In[1]
        else:
            T=In[0]
            q=In[1]
        X=RP.REFPROPdll( "","TQ", Out0, Units,iMass, qflag, T, q, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # q and T
    if Gr[0]=="q" and Gr[1]=="T":
        if Eh=="CBar":
            T=In[1]+dTk
            q=In[0]
        elif Eh=="CKPa":
            T=In[1]+dTk
            q=In[0]
        else:
            T=In[1]
            q=In[0]
        X=RP.REFPROPdll("", "TQ", Out0, Units,iMass, qflag, T, q, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
##############################################################
    # p and q
    if Gr[0]=="p" and Gr[1]=="q":
        if Eh=="CBar":
            p=In[0]*1e5
            q=In[1]
        elif Eh=="CKPa":
            p=In[0]*1e3
            q=In[1]
        else:
            p=In[0]
            q=In[1]
        X=RP.REFPROPdll("", "PQ", Out0, Units,iMass, qflag, p, q, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q

    # q and p
    if Gr[0]=="q" and Gr[1]=="p":
        if Eh=="CBar":
            p=In[1]*1e5
            q=In[0]
        elif Eh=="CKPa":
            p=In[1]*1e3
            q=In[0]
        else:
            p=In[1]
            q=In[0]
        X=RP.REFPROPdll("", "PQ", Out0, Units,iMass, qflag, p, q, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
##############################################################      
    # p and h
    if Gr[0]=="p" and Gr[1]=="h":
        if Eh=="CBar":
            p=In[0]*1e5
            h=In[1]*1000.
        elif Eh=="CKPa":
            p=In[0]*1e3
            h=In[1]*1000.
        else:
            p=In[0]
            h=In[1]
        X=RP.REFPROPdll("", "PH", Out0, Units,iMass, qflag, p, h, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
    
    # h and p
    if Gr[0]=="h" and Gr[1]=="p":
        if Eh=="CBar":
            p=In[1]*1e5
            h=In[0]*1000.
        elif Eh=="CKPa":
            p=In[1]*1e3
            h=In[0]*1000.
        else:
            p=In[1]
            h=In[0]
        X=RP.REFPROPdll("", "PH", Out0, Units,iMass, qflag, p, h, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
##############################################################
    # p and s
    if Gr[0]=="p" and Gr[1]=="s":
        if Eh=="CBar":
            p=In[0]*1e5
            s=In[1]*1000.
        elif Eh=="CKPa":
            p=In[0]*1e3
            s=In[1]*1000.
        else:
            p=In[0]
            s=In[1]
        X=RP.REFPROPdll("", "PS", Out0, Units,iMass, qflag, p, s,z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
    
    # s and p
    if Gr[0]=="h" and Gr[1]=="p":
        if Eh=="CBar":
            p=In[1]*1e5
            s=In[0]*1000.
        elif Eh=="CKPa":
            p=In[1]*1e3
            s=In[0]*1000.
        else:
            p=In[1]
            s=In[0]
        X=RP.REFPROPdll("", "PH", Out0, Units,iMass, qflag, p, s, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
##############################################################

    # T und v

    if Gr[0]=="T" and Gr[1]=="v":
        if Eh=="CBar":
            T=In[0]+dTk
            v=In[1]
        elif Eh=="CKPa":
            T=In[0]+dTk
            v=In[1]
        else:
            T=In[0]
            v=In[1]
        X=RP.REFPROPdll("", "TD", Out0, Units,iMass, qflag, T, 1./v,z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # v und T
    if Gr[0]=="v" and Gr[1]=="T":
        if Eh=="CBar":
            T=In[1]+dTk
            v=In[0]
        elif Eh=="CKPa":
            T=In[1]+dTk
            v=In[0]
        else:
            T=In[1]
            v=In[0]
        X=RP.REFPROPdll("", "TD", Out0, Units,iMass, qflag, T, 1./v, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
##############################################################
    # p und v

    if Gr[0]=="p" and Gr[1]=="v":
        if Eh=="CBar":
            p=In[0]*1e5
            v=In[1]
        elif Eh=="CKPa":
            p=In[0]#1e3
            v=In[1]
        else:
            p=In[0]
            v=In[1]
        X=RP.REFPROPdll("", "pD", Out0, Units,iMass, qflag, p, 1./v, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # v und p
    if Gr[0]=="v" and Gr[1]=="T":
        if Eh=="CBar":
            p=In[1]*1e5
            v=In[0]
        elif Eh=="CKPa":
            p=In[1]*1e3
            v=In[0]
        else:
            p=In[1]
            v=In[0]
        X=RP.REFPROPdll("", "pD", Out0, Units,iMass, qflag, p, 1./v,z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q

####################################################################
    # u und v

    if Gr[0]=="u" and Gr[1]=="v":
        if Eh=="CBar":
            u=In[0]*1000
            v=In[1]
        elif Eh=="CKPa":
            u=In[0]*1000.
            v=In[1]
        else:
            u=In[0]
            v=In[1]
        X=RP.REFPROPdll("", "ED", Out0, Units,iMass, qflag, u, 1./v, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # v und p
    if Gr[0]=="v" and Gr[1]=="T":
        if Eh=="CBar":
            u=In[1]*1000
            v=In[0]
        elif Eh=="CKPa":
            u=In[1]*1000
            v=In[0]
        else:
            u=In[1]
            v=In[0]
        X=RP.REFPROPdll("", "ED", Out0, Units,iMass, qflag, p, 1./v, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
####################################################################
    # h und s

    if Gr[0]=="h" and Gr[1]=="s":
        if Eh=="CBar":
            h=In[0]*1000.
            s=In[1]*1000.
        elif Eh=="CKPa":
            h=In[0]*1000.
            s=In[1]*1000.
        else:
            h=In[0]
            s=In[1]
        X=RP.REFPROPdll("", "HS", Out0, Units,iMass, qflag, h, s, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q
        
    # s und h
    if Gr[0]=="s" and Gr[1]=="h":
        if Eh=="CBar":
            h=In[1]*1000.
            s=In[0]*1000.
        elif Eh=="CKPa":
            h=In[1]*1000.
            s=In[0]*1000.
        else:
            h=In[1]
            s=In[0]
        X=RP.REFPROPdll("", "ED", Out0, Units,iMass, qflag, h, s, z)
        T,p,v,u,h,s,BETA=X.Output[0:7]
        q=X.q  
        
    
    
    # Quality output correction
    if q>1.:
        q=5
    elif q<0.:
        q=-5

    # Calculation of transport properties if required
    if ("dvis" in Out or "kvis" in Out or "lam" in Out or "Pr" in Out or "cp0" in Out or "cp" in Out) and (q==5 or q==-5):
        if Eh=="CBar":
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TP", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, p, z).Output[0:5]
        elif Eh=="CKPa":
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TP", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, p,z).Output[0:5]
        else:
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TP", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, p, z).Output[0:5]
        kvis=dvis*v
    elif q==1. or q==0.:
        if Eh=="CBar":
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TQ", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, q, z).Output[0:5]
        elif Eh=="CKPa":
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TQ", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, q,z).Output[0:5]
        else:
            dvis, lam, Pr,cp0,cp=RP.REFPROPdll("", "TQ", "ETA;TCX;PRANDTL;CP0;CP", Units,iMass, qflag, T, q, z).Output[0:5]
        kvis=dvis*v
        
    else:
        dvis=kvis=Pr=lam=cp0=cp=nan
    

    if Eh=="mol":
        kvis=kvis/M
    
     # Output unit correction
    if Eh=="CBar":
            T=T-dTk
            p=p*1e-5
            u=u/1000.
            h=h/1000.
            s=s/1000.
            cp=cp/1000.
            cp0=cp0/1000.
            lam=lam/1000.
     
    elif Eh=="CKPa":
            T=T-dTk
            p=p*1e-3
            u=u/1000.
            h=h/1000.
            s=s/1000.
            cp=cp/1000.
            cp0=cp0/1000.
            lam=lam/1000
    
    
    

    z=[T,p,v,u,h,s,BETA,q,dvis,kvis,lam,Pr,cp0,cp,M]
    erg=zeros(len(Out))
    for k in range(len(Out)):
        erg[k]=z[Gr0.index(Out[k])]
    
    
    return erg

def getInfo(name):

    z=[1.0]
    if ".MIX" in name:
        z=RP.SETMIXTUREdll(name)[0]
    else:
        
        z=[1.0]
    Units = RP.GETENUMdll(0, "MASS BASE SI").iEnum

    M,Tc,Tmin,Tmax,pc=RP.REFPROPdll("","","M,TC,TMIN,TMAX,PC",Units,0,0,0,0,z).Output[0:5] #Molmasse kg/mol

    return [M,Tc-273.15,Tmin-273.15,Tmax-273.15,pc*1e-3]
    
def getNames():
    names=["13BUTADIENE","1BUTYNE","1PENTENE","22DIMETHYLBUTANE","23DIMETHYLBUTANE","3METHYLPENTANE",
       "acetone","acetylene","ammonia","argon","benzene","butane","1BUTENE",
       "CO2","CO","COS","chlorine","chlorobenzene","C2BUTENE","cyclobutene","CYCLOHEX","CYCLOPEN","CYCLOPRO",
       "D4","D5","D6","DEA","decane","D2","R150","DEE","DMC",
       "DME","C22","C12","ethane","ethanol","EGLYCOL","EBENZENE","ethylene","ETHYLENEOXIDE","fluorine","D2O","helium","heptane","C16","hexane","HYDROGEN","HCL",
       "H2S","ISOBUTAN","IBUTENE","IHEXANE","IOCTANE","IPENTANE","krypton","MD2M","MD3M","MD4M","MDM","MEA","methane","methanol","MLINOLEA","MLINOLEN","MOLEATE",
       "MPALMITA","MSTEARAT","C1CC6","MM","MXYLENE","neon","NEOPENTN","nitrogen","NF3","N2O","nonane","NOVEC649","octane","ORTHOHYD","oxygen","OXYLENE","PARAHYD",
       "pentane","C4F10","C6F14","C5F12","propadiene","propane","C3CC6","PROPYLEN","PROPYLENEOXIDE","propyne","PXYLENE","R11","R1123","R113","R114","R115","R116",
       "R12","R1216","R1224YDZ","R123","R1233ZDE","R1234yf","R1234ZEE","R1234ZEZ","R124","R1243zf","R125","R13","R1336MZZZ","R134a","CF3I","R14","R141b","R142b",
       "R143a","R152a","R161","R21","R218","R22","R227ea","R23","R236ea","R236fa","R245ca","R245fa","R32","R365mfc","R40","R41","RC318","RE143a","RE245cb2",
       "RE245fa2","RE347MCC","SO2","SF6","toluene","T2BUTENE","C11","VINYLCHLORIDE","water","xenon"]
    return names    


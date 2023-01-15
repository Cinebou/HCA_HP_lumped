import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fluidProp import VLEFluid
from math import exp, log
plt.rcParams["font.size"] = 23
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 100

def main():
    mof = DA_model()
    #mof.v_ads_plot()
    #mof.isotherm_plot()
    #mof.Qst_plot()
    print(mof.current_h_ads(15,298.15))
    


class DA_model():
    def __init__(self):
        self.co = 0.380 * 1000 * 1e-6# m3/kg_MOF = cm3/g * 1000 * 1e-6, Belsorp experiment
        # co2 体積、0.000638978　m3/kg

        # isotherm params
        self.n = 1.42040287e+00
        self.E = 5.32025217e+03 # J/mol
        # thermal expansion
        self.alpha = 2.97307623e-03
        # saturation over critical point coefficient
        self.k_over_crit = 5.92587884e+00

        self.co2_prop = VLEFluid('CO2')
        self.__coeff()


    def __coeff(self):
        # gas const
        self.R = 8.314 # J/K/mol
        # critical point
        self.T_crit = 30 + 273.16  # critical temperature, 30℃
        
        # boiling point data, wikipedia
        self.T_boiling = 194.65 # K
        self.p_boiling = 1e5 # Pa
        self.v_boiling = 1/(1.565*1e3) # m3/kg_CO2, Icarus Volume 294, 15 September 2017, Pages 201-208
        self.v_gas_boiling = 0.36071 # m3/kg



    # accumurative heat
    def accum_H_ads(self, T):
        pres_point = 100
        sat_P = self.sat_pressure(T)
        p_line = np.linspace(1e3, max(sat_P-1e6,4.2e6), pres_point)

        H_ads_list = [0]
        for i in range(1, len(p_line)):
            dm_ads = self.eq_loading(p_line[i],T) - self.eq_loading(p_line[i-1],T)
            Qst = self.isosteric_heat((p_line[i]+p_line[i-1])/2, T)
            H_ads_list.append(H_ads_list[-1] + Qst*dm_ads)
        return np.array(p_line), np.array(H_ads_list)


    # return the accumurative heat at specific temperature and pressure
    def accum_H_pT(self,p,T):
        p_l, H_l = self.accum_H_ads(T)
        index = self.idx_of_the_nearest(p_l,p)
        return H_l[index]


    # isotherm DA model
    def eq_loading(self, p, T):
        sat_P = self.sat_pressure(T)
        theta = exp(-((self.R * T) / self.E * log(sat_P/p))**self.n)
        _, v_ads = self.vol_gas_adsobate(p, T)
        w = self.co/v_ads * theta * 1000/44 # mol / kg_MOF
        return w


    # saturation pressure
    def sat_pressure(self, T):
        if 194 < T < 196:
            return 0.10397 * 1e6 # Pa
        elif T < self.T_crit:
            return self.co2_prop.calc_VLE_T(T).p_v
        else:
            p_crit = self.co2_prop.calc_VLE_T(self.T_crit).p_v
            ps = (T/self.T_crit)**self.k_over_crit * p_crit
            return ps


    # 飽和蒸気圧の温度微分  
    def TdlnPs_aT(self, T):
        if T >= self.T_crit:
            return self.k_over_crit
        else:
            dT = 0.1 # K
            Ps_T = self.sat_pressure(T)
            Ps_T_dT = self.sat_pressure(T+dT)
            dPsdT = (Ps_T_dT - Ps_T) / dT
            #print('飽和蒸気圧での傾き',dPsdT)
            return dPsdT * T / Ps_T


    # volume of gas and adsorbate
    def vol_gas_adsobate(self, p, T):
        # volume of gas and adsorbate
        v_gas = self.co2_prop.calc_fluidProp_pT(p,T).vol # m3/kg
        
        v_ads = self.v_boiling * exp(self.alpha * (T - self.T_boiling)) # m3/kg_CO2
        return v_gas, v_ads


    #微分吸着熱の計算, Azaharのモデル,  kJ/mol_co2で返す
    def isosteric_heat(self, p, T):
        sat_P = self.sat_pressure(T)
        v_gas, v_ads = self.vol_gas_adsobate(p, T)
        
        dlnP_dT = self.TdlnPs_aT(T)
        dlnP_dT += log(sat_P/p)
        dlnP_dT += (self.alpha/self.n) * (self.E/self.R)**self.n * (T * log(sat_P/p))**(1-self.n)
        Qst = p * (v_gas - v_ads) * dlnP_dT  # J/kg
        return Qst * 0.001 * 44/1000 # kJ/mol

    
    def current_h_ads(self, loading, T):
        p_line = np.linspace(3e4, self.sat_pressure(T)*0.98,15)
        m_line = np.array([self.eq_loading(p,T) for p in p_line])
        qst_line = np.array([self.isosteric_heat(p,T) for p in p_line])

        # index of correspoding heat of adsorption
        index = self.idx_of_the_nearest(m_line, loading)
        return qst_line[index]


    
    #　吸着相の比熱、液体比熱近似すると臨界点で発散するよ
    # Kazi Afzalur Rahman et al. / Procedia Engineering 56 ( 2013 ) 118 – 125
    def specific_heat_Cap(self,p,T, loading):
        loading = loading * 44/1000
        # 被吸着率
        _, v_ads = self.vol_gas_adsobate(p, T)
        theta = loading * v_ads / self.co
        #ガス相比熱
        cp_gas = self.co2_prop.calc_fluidProp_pT(p,T).cp * 44/1000# J/mol/K
        term2 = (self.alpha**2*(1-self.n)/self.n**2*self.E*T) * (log(1/theta))**(1/self.n-2) # J/mol/K
        #気液変化の凝縮熱の項は圧力一定なのでなし
        #term3 = -self.TdlnPs_aT(T) * self.R # J/mol/K
        cp_ads = cp_gas + term2 
        #print(cp_ads, '\t', term3)
        return cp_ads # J/mol/K


    # plot accum heat
    def accum_plot(self):
        T = 303.15
        p_line, Q_line = self.accum_H_ads(T)
        plt.plot(p_line, Q_line)
        plt.xlim(0,4.1e6)
        plt.ylim(0,max(Q_line)+30)
        plt.show()


    def v_ads_plot(self):
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111)
        T_line = np.linspace(195, 320, 50)
        v_ads = [self.v_boiling * exp(self.alpha * (T - self.T_boiling))*1e3 for T in T_line]# m3/kg
        ax.plot(T_line, v_ads,label='$v_\mathrm{ads}$')

        ax.legend()
        ax.set_ylim(0,2)
        ax.set_xlabel('$T$ [K]')
        ax.set_ylabel('$v$ [$\mathrm{cm}^3/\mathrm{g}$]')
        fig.savefig('./Fig/v_ads_liq.png')
        

    # plot Qst
    def isotherm_plot(self):
        fig = plt.figure(figsize=(10,7.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.14,0.14,0.8,0.8])
        T = 288
        sat_P = self.sat_pressure(T)
        p_line = np.linspace(1e2, sat_P*0.99,100)
        w_line = np.array([self.eq_loading(p,T) for p in p_line])
        ax.plot(p_line/1e6, w_line,label='288K',linewidth=4)

        T = 303
        sat_P = self.sat_pressure(T)
        p_line = np.linspace(1e2, sat_P*0.99,100)
        w_line = np.array([self.eq_loading(p,T) for p in p_line])
        ax.plot(p_line/1e6, w_line,label='303K',linewidth=4)

        T = 313
        sat_P = self.sat_pressure(T)
        p_line = np.linspace(1e2, sat_P*0.99,100)
        w_line = np.array([self.eq_loading(p,T) for p in p_line])
        ax.plot(p_line/1e6, w_line,label='313K',linewidth=4)
        ax.set_xlim(0,4.1)
        ax.set_ylim(0,12)
        ax.legend()
        ax.set_xlabel('$p$ [MPa]')
        ax.set_ylabel('$M$ [mmol/g]')
        fig.savefig('./Fig/eq_loading.png')


    # plot Qst
    def Qst_plot(self):
        fig = plt.figure(figsize=(10,7.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.14,0.14,0.8,0.8])
        T = 288
        p_line = np.linspace(3e4, self.sat_pressure(T)*0.98,60)
        m_line = np.array([self.eq_loading(p,T) for p in p_line])
        qst_line = np.array([self.isosteric_heat(p,T) for p in p_line])
        ax.scatter(m_line, qst_line,label='288K',s=65)

        T = 303
        p_line = np.linspace(3e4, self.sat_pressure(T)*0.98,60)
        m_line = np.array([self.eq_loading(p,T) for p in p_line])
        qst_line = np.array([self.isosteric_heat(p,T) for p in p_line])
        ax.scatter(m_line, qst_line,label='303K',s=65)

        T = 313
        p_line = np.linspace(3e4, self.sat_pressure(T)*0.98,60)
        m_line = np.array([self.eq_loading(p,T) for p in p_line])
        qst_line = np.array([self.isosteric_heat(p,T) for p in p_line])
        ax.scatter(m_line, qst_line,label='313K',s=65)
        ax.set_ylim(0,32)
        ax.set_xlim(0,9.5)
        ax.set_ylabel("$Q_\mathrm{ads}$ [kJ/mol]")
        ax.set_xlabel("$m_\mathrm{ads}$ [mmol/g]")
        ax.legend(loc='lower left', borderaxespad=0.5)
        fig.savefig('./Fig/Q_ads_heat.png')


    # plot saturation pressure 
    def p_plot(self):
        T_line = np.linspace(240,323.15,100)
        p_line = np.array([self.sat_pressure(t) for t in T_line])
        plt.plot(T_line, p_line)
        pT = pd.DataFrame([T_line] + [p_line]).T
        pT.to_csv('./sat_P/sat_P.csv',header=None,index=None)
        plt.show()


    # volume of gas and adsorbate
    def v_plot(self):
        T = 298
        p_line = np.linspace(1e5,3e6,100)
        v_line = np.array([self.vol_gas_adsobate(p, T) for p in p_line]).T
        plt.plot(p_line,v_line[1])
        plt.show()

    # plot Qst
    def sat_P_plot(self):
        fig = plt.figure(figsize=(10,7.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.14,0.14,0.8,0.8])
        T_line = np.linspace(263, 323,100)
        sat_P_line = np.array([self.sat_pressure(T) for T in T_line])

        T_line_ref = np.linspace(263, 303,100)
        sat_P_ref = np.array([self.co2_prop.calc_VLE_T(T).p_v for T in T_line_ref])

        ax.plot(T_line, sat_P_line,label='approx.',linewidth=3.5)
        ax.scatter(T_line_ref, sat_P_ref, label = 'vapor liquid',s=40,color = 'red')
        ax.legend(loc='lower right')
        ax.set_ylabel('$p$ [MPa]')
        ax.set_xlabel('$T$ [K]')
        
        T = 313
        print('temp\t:\t',T)
        print('pres\t:\t',self.sat_pressure(T))
        fig.savefig('./Fig/sat_P_approx.png')
    

    # find index near to the experiment point
    def idx_of_the_nearest(self, data, value):
        idx = np.argmin(np.abs(np.array(data) - value))
        return idx



if __name__ == '__main__':
    main()

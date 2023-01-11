"""
Created on Tue 01 Feb

@author: Hibiki Kimura
"""

import logging
from os import mkdir, path, remove


""" 
The results or Logger can be summarized in Log/ folder
Each Logging can be turned off by 'lg.setlevel(logging.WARN)'
"""

if not path.exists('./Log'):
    mkdir('./Log')

class Log:
        """ Define the Logging settings, log file, format, Log on/off
        """
        def __init__(self):
                log_list = ["Log/log_mass.csv", "Log/log_heat_gas.csv","Log/log_any.csv","Log/log_heat_sor.csv"]
                for f in log_list:
                        if path.exists(f):
                                remove(f)

                # Mass CO2 flow
                lg_m = logging.getLogger('mass_balance')
                handler_m = logging.FileHandler(filename = "Log/log_mass.csv",mode = 'a')
                handler_m.setFormatter(logging.Formatter("%(message)s"))
                lg_m.setLevel(logging.DEBUG) # log on
                #lg_m.setLevel(logging.WARN)   # log off
                lg_m.addHandler(handler_m)

                # Heat of gas
                lg_h = logging.getLogger('heat_balance_gas')
                handler_h = logging.FileHandler(filename = "Log/log_heat_gas.csv",mode = 'a')
                handler_h.setFormatter(logging.Formatter("%(message)s"))
                lg_h.setLevel(logging.DEBUG) # log on
                #lg_h.setLevel(logging.WARN)   # log off
                lg_h.addHandler(handler_h)

                # Heat of sorbent
                lg_hs = logging.getLogger('heat_balance_sor')
                handler_hs = logging.FileHandler(filename = "Log/log_heat_sor.csv",mode = 'a')
                handler_hs.setFormatter(logging.Formatter("%(message)s"))
                lg_hs.setLevel(logging.DEBUG) # log on
                #lg_hs.setLevel(logging.WARN)   # log off
                lg_hs.addHandler(handler_hs)

                # Any
                lg_any = logging.getLogger('Any')
                handler_any = logging.FileHandler(filename = "Log/log_any.csv",mode = 'a')
                handler_any.setFormatter(logging.Formatter("%(message)s"))
                lg_any.setLevel(logging.DEBUG) # log on
                #lg_any.setLevel(logging.WARN)   # log off
                lg_any.addHandler(handler_any)

                """set the logger"""
                Log.__instance = self     


        def log_set_m(self,msg):
                log = logging.getLogger('mass_balance')
                log.debug(msg)

        def log_set_h_gas(self,msg):
                log = logging.getLogger('heat_balance_gas')
                log.debug(msg)

        def log_set_h_sor(self,msg):
                log = logging.getLogger('heat_balance_sor')
                log.debug(msg)

        def log_set_any(self,msg):
                log = logging.getLogger('Any')
                log.debug(msg)


""" creating instance of Logger """
l = Log()

""" output any message to "Log/log_mass.csv" file """
def log_mass_msg(m_in, m_ads, m_out,t):
        msg = '{},{},{},{}'.format( m_in,m_ads, m_out,t)
        l.log_set_m(msg)
        return 0

""" output any message to "Log/log_heat_gas.csv" file """
def log_h_gas_msg(mcpdTgasdt, en_flow, Htrans, mRT, t):
        msg = '{},{},{},{},{}'.format(mcpdTgasdt, en_flow, Htrans, mRT,t)
        l.log_set_h_gas(msg)
        return 0

""" output any message to "Log/log_heat_sor.csv" file """
def log_h_sor_msg(dTdt, Htrans,Hads,  passHTF,t):
        msg = '{},{},{},{},{}'.format(dTdt, Htrans, Hads, passHTF,t)
        l.log_set_h_sor(msg)
        return 0


""" output any message to "Log/log_any.csv" file """
def log_any_msg(msg):
        l.log_set_any(msg)
        return 0

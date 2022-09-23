# Requirement
REFPROP version 10
'c:/Program Files (x86)/REFPROP'
### Not version 9 !!!


# How to run simulation
in the derectory
```
python Eq_balance.py
```
The graph of the temperature histograph can be found in './Fig'

# initial condition
Initial condition of the simulation can be modified in '__IC()' of 'HC_AHP_System.py', line 106 ~ 118.

# Needs to be refined
1. The heat capacity of copper tube is not included in this model.
2. The heat capacity of the adsorbed CO2 should be depend on the temperature, not constant.
3. Heat of adsorption should be defferent with increasing adsorption amount.
4. Adsorption speed should be dependent on temperature and pressure. It can be divided into surface diffusion and knudsen diffusion and internal particle diffusion.
5. Heat transfer coefficient can be estimated from Reynold's number and Prank number. 
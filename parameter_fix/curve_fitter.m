de_data = readmatrix('desorption0C.csv');
p_de = de_data(:,1);
cw_de = de_data(:,2);

ad_data = readmatrix('adsorption0C.csv');
p_ad = ad_data(:,1);
cw_ad = ad_data(:,2);

sat_data = readmatrix('sat_pressure_water.csv');
T_sat = sat_data(1:100,1);
P_sat = sat_data(1:100,2).*1e6;
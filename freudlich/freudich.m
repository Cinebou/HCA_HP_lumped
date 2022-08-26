isotherm = readmatrix('isotherm.csv');
p_sat = 64.3;
p = isotherm(:,1);
loading = isotherm(:,2);

p_par = p./p_sat;
reg = polyfit(log(p_par),log(loading),1);

n = 1/reg(1);
k = exp(reg(2));

p_point = 0.01:0.01:1;
m_fre = k .* p_point.^(1/n);

figure(1)
plot(p_point, m_fre)
hold on
scatter(p_par, loading)
hold off
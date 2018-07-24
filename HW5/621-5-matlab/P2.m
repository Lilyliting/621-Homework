%% 
S0=100;K=100;V0=0.010201;kappa=6.21;theta=0.019;sigma=0.61;rho=-0.7;
n=300;m=1e6;Tm=1;r=0.0319;
[Price, Bias, RMSE, Time] = MC_Heston(1, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)
[Price, Bias, RMSE, Time] = MC_Heston(2, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)
[Price, Bias, RMSE, Time] = MC_Heston(3, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)
[Price, Bias, RMSE, Time] = MC_Heston(4, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)
[Price, Bias, RMSE, Time] = MC_Heston(5, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)

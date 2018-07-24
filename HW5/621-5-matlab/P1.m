%% 
iscall=1;S0=100;K=100;Tm=1;r=0.06;sigma=0.2;div=0.03;beta=-1;n=300;m=1e6;
[Price, SD, SE, Time] = Monte_Carlo(iscall, S0, K, Tm, r, sigma, div, n, m)
%% 
[Price, SD, SE, Time] = Monte_Carlo_AVR(iscall, S0, K, Tm, r, sigma, div, n, m)
%% 
[Price, SD, SE, Time] = Monte_Carlo_DC(iscall, S0, K, Tm, r, sigma, div, beta, n, m)
%% 
[Price, SD, SE, Time] = Monte_Carlo_AVRDC(iscall, S0, K, Tm, r, sigma, div, beta, n, m)
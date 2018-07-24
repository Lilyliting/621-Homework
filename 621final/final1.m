%% a
S0=100; K=100; Tm=5; r=0.03; sigma=0.3;
[Price_G] = Black_Scholes_Asian(S0, K, Tm, r, sigma)

%% b,c
% Arithmetic Asian Call options
type=1; n=252*5; m=1e6; cl=0.99;
[Price_A_sim, lower_bound, upper_bound, Time] = Monte_Carlo_Asian(type, S0, K, Tm, r, sigma, n, m, cl)

% Geometric Asian Call options
type=2;
[Price_G_sim, lower_bound, upper_bound, Time] = Monte_Carlo_Asian(type, S0, K, Tm, r, sigma, n, m, cl)

%% d
m=1e4;
[Price_A, Price_G, b_star] = Monte_Carlo_Asian_duo(S0, K, Tm, r, sigma, n, m)
m=1e5;
[Price_A, Price_G, b_star] = Monte_Carlo_Asian_duo(S0, K, Tm, r, sigma, n, m)
m=5e5;
[Price_A, Price_G, b_star] = Monte_Carlo_Asian_duo(S0, K, Tm, r, sigma, n, m)
m=1e6;
[Price_A, Price_G, b_star] = Monte_Carlo_Asian_duo(S0, K, Tm, r, sigma, n, m)

%% e
error_G = Price_G - Price_G_sim

%% f
Price_A_star = Price_A_sim - b_star*error_G;

S0=100; K=100; Tm=5; r=0.03; sigma=0.3;
[Price_G] = Black_Scholes_Asian(S0, K, Tm, r, sigma)
type=2; n=252*5; m=1e6; cl=0.99;
[Price_G_sim, lower_bound, upper_bound, Time] = Monte_Carlo_Asian(type, S0, K, Tm, r, sigma, n, m, cl)
error_G = Price_G - Price_G_sim
function [Price, SD, SE, Time]=Monte_Carlo_AVR(iscall, S0, K, Tm, r, sigma, div, n, m)
% Monte Carlo valuation with Antithetic Variance reduction
% iscall=1: european call option; iscall=0: put option
% S0: Initial stock price
% K: Strike price
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility
% delta: Continuous dividend rate
% n: Time steps
% m: Quantity of trials
tic;
dt = Tm/n;
S0 = ones(1,m).*S0;
ST1 = Path_to_T(1, S0, r, sigma, div, dt, n, m);
ST2 = Path_to_T(-1, S0, r, sigma, div, dt, n, m);

if iscall==1
    profits = (max(ST1-K,0)+max(ST2-K,0))/2*exp(-r*Tm);
else
    profits = (max(K-ST1,0)+max(K-ST2,0))/2*exp(-r*Tm);
end

Price = mean(profits);
SD = std(profits);
SE = SD/sqrt(m);
Time=toc;
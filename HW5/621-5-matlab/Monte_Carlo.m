function [Price, SD, SE, Time] = Monte_Carlo(iscall, S0, K, Tm, r, sigma, div, n, m)
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
S_T = Path_to_T(1, S0, r, sigma, div, dt, n, m); % Stock prices at maturity

if iscall==1
    profits = max(S_T-K,0)*exp(-r*Tm);
else
    profits = max(K-S_T,0)*exp(-r*Tm);
end

Price = mean(profits);
SD = std(profits);
SE = SD/sqrt(m);
Time=toc;
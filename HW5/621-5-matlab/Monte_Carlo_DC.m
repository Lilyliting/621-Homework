function [Price, SD, SE, Time]=Monte_Carlo_DC(iscall, S0, K, Tm, r, sigma, div, beta, n, m)
% Monte Carlo valuation with delta control
% iscall=1: european call option; iscall=0: put option
% S0: Initial stock price
% K: Strike price
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility
% div: Continuous dividend rate
% n: Time steps
% m: Quantity of trials
tic;

dt = Tm/n;
nudt = (r-div-sigma^2/2)*dt;
sigsdt = sigma*sqrt(dt);
erddt = exp((r-div)*dt);

S0 = ones(1,m).*S0;

St = S0;
cv = zeros(1,m);


for i=1:n
    t=(i-1)*dt;
    d1 = (log(St/K)+(r-div+sigma^2/2)*(Tm-t))/sigma/sqrt(Tm-t);
    if iscall==1 % Black-Scholes Formulas for Delta
        delta = exp(-div*(Tm-t))*normpdf(d1,0,1);
    else
        delta = exp(-div*(Tm-t))*(normpdf(d1,0,1)-1);
    end
    epsilon = normrnd(0,1,1,m);
    Stn = St.*exp(nudt+sigsdt.*epsilon); % evolve the stock price
    cv = cv+delta.*(Stn-St*erddt);
    St = Stn;
end

if iscall==1
    profits = max(St-K,0)+beta*cv;
else
    profits = max(K-St,0)+beta*cv;
end

Price = mean(profits)*exp(-r*Tm);
SD = std(profits);
SE = SD/sqrt(m);
Time=toc;
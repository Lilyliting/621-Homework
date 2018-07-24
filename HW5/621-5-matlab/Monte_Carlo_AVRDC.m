function [Price, SD, SE, Time]=Monte_Carlo_AVRDC(iscall, S0, K, Tm, r, sigma, div, beta, n, m)
% Monte Carlo valuation with Antithetic Variance reduction and delta control
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

St1 = S0;
St2 = S0;
cv1 = zeros(1,m);
cv2 = zeros(1,m);

for i=1:n
    t=(i-1)*dt;
    d1 = (log(St1/K)+(r-div+sigma^2/2)*(Tm-t))/sigma/sqrt(Tm-t);
    d2 = (log(St2/K)+(r-div+sigma^2/2)*(Tm-t))/sigma/sqrt(Tm-t);
    if iscall==1 % Black-Scholes Formulas for Delta
        delta1 = exp(-div*(Tm-t))*normpdf(d1,0,1);
        delta2 = exp(-div*(Tm-t))*normpdf(d2,0,1);
    else
        delta1 = exp(-div*(Tm-t))*(normpdf(d1,0,1)-1);
        delta2 = exp(-div*(Tm-t))*(normpdf(d2,0,1)-1);
    end

    epsilon = normrnd(0,1,1,m);
    Stn1 = St1.*exp(nudt+sigsdt.*epsilon); % evolve the stock price
    Stn2 = St2.*exp(nudt-sigsdt.*epsilon);
    cv1 = cv1+delta1.*(Stn1-St1*erddt);
    cv2 = cv2+delta2.*(Stn2-St2*erddt);
    St1 = Stn1;
    St2 = Stn2;
end

if iscall==1
    profits = (max(St1-K,0)+beta*cv1 + max(St2-K,0)+beta*cv2)/2;
else
    profits = (max(K-St1,0)+beta*cv1 + max(K-St2,0)+beta*cv2)/2;
end

Price = mean(profits)*exp(-r*Tm);
SD = std(profits);
SE = SD/sqrt(m);
Time=toc;
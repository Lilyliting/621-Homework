function [Price] = Black_Scholes_Asian(S0, K, Tm, r, sigma)
% geometric Asian Call option in the Black-Scholes model
% S0: Initial stock price
% K: Strike price
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility

N = Tm*252;
sigmahat = sigma*sqrt((2*N+1)/6/(N+1));
rho = (r-sigma^2/2+sigmahat^2)/2;
d1 = (log(S0/K) + (rho+sigmahat^2/2)*Tm)/sqrt(Tm)/sigmahat;
d2 = (log(S0/K) + (rho-sigmahat^2/2)*Tm)/sqrt(Tm)/sigmahat;
Price = exp(-r*Tm)*(S0*exp(rho*Tm)*normcdf(d1)-K*normcdf(d2));

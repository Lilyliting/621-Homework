function [Price]=Basket_option(iscall,S0,K,Tm,dt,A,sigma,MU,m)
% iscall=1: european call option; iscall=0: put option
% S0: Initial stock prices
% K: Strike price
% Tm: Time to maturity
% dt: Time intervel
% A: Correlation matrix
% sigma: Volatilities
% MU: Risk free rates
% m: Quantity of trials

St = Correlated_BM(S0,Tm,dt,A,sigma,MU,m);

Ut = mean(St,3); % a simple average basket
UT = Ut(end,:);
Price = mean(max((-1)^iscall*(K-UT),0));
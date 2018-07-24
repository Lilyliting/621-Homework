function [S_T]=Path_to_T(symbol, S0, r, sigma, div, dt, n, m)
% symbol=1 means that the symbol of standard normal item is positive
% symbol=-1 means negative

lnS_T = log(S0);
nudt = (r-div-sigma^2/2)*dt;
sigsdt = sigma*sqrt(dt);

for i=1:n
    epsilon = symbol*normrnd(0,1,1,m);
    dS = nudt+sigsdt.*epsilon; % evolve the stock price
    lnS_T = lnS_T+dS;
end
S_T = exp(lnS_T);

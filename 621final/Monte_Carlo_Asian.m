function [Price, lower_bound, upper_bound, Time] = Monte_Carlo_Asian(type, S0, K, Tm, r, sigma, n, m, cl)
% Type=1: arithmetic Asian call options; Type=2: geometric Asian Call
% options
% S0: Initial stock price
% K: Strike price
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility
% n: Time steps
% m: Quantity of trials
% cl: Confidence level 

tic;
dt = Tm/n;
S0 = ones(1,m).*S0;

if type==1
    sumST = S0;
else
    productST = S0.^(1/(n+1));
end

lnS_T = log(S0);
nudt = (r-sigma^2/2)*dt;
sigsdt = sigma*sqrt(dt);

for i=1:n
    epsilon = normrnd(0,1,1,m);
    dS = nudt+sigsdt.*epsilon; % evolve the stock price
    lnS_T = lnS_T+dS;
    if type==1
        sumST = sumST + exp(lnS_T);
    else
        productST = productST.*exp(lnS_T).^(1/(n+1));
    end
end


if type==1
    profits = max(sumST/(n+1) - K, 0);
else
    profits = max(productST - K, 0);
end

Price = mean(profits)*exp(-r*Tm);
SD = std(profits);
SE = SD/sqrt(m);

z=-norminv((1-0.95)/2);
lower_bound = Price - z*SE;
upper_bound = Price + z*SE;

Time=toc;



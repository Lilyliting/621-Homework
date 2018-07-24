function [Price_A, Price_G, b_star] = Monte_Carlo_Asian_duo(S0, K, Tm, r, sigma, n, m)
% Pricing arithmetic and geometric Asian Call with same random variables
% options
% S0: Initial stock price
% K: Strike price
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility
% delta: Continuous dividend rate
% n: Time steps
% m: Quantity of trials
dt = Tm/n;
S0 = ones(1,m).*S0;

sumST = S0;
productST = S0.^(1/(n+1));

lnS_T = log(S0);
nudt = (r-sigma^2/2)*dt;
sigsdt = sigma*sqrt(dt);

for i=1:n
    epsilon = normrnd(0,1,1,m);
    dS = nudt+sigsdt.*epsilon; % evolve the stock price
    lnS_T = lnS_T+dS;
    sumST = sumST + exp(lnS_T);
    productST = productST.*exp(lnS_T).^(1/(n+1));
end

Yi = max(sumST/(n+1) - K, 0);
Xi = max(productST - K, 0);

Price_A = mean(Yi)*exp(-r*Tm);
Price_G = mean(Xi)*exp(-r*Tm);

Xbar = mean(Xi); Ybar = mean(Yi);
b_star = sum((Xi-Xbar).*(Yi-Ybar))/sum((Xi-Xbar).^2);


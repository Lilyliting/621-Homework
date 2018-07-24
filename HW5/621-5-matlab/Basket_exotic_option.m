function [Price]=Basket_exotic_option(S0,K,B,Tm,dt,A,sigma,MU,m)
% For european call option
% S0: Initial stock prices
% K: Strike price
% B: Barrier
% Tm: Time to maturity
% dt: Time intervel
% A: Correlation matrix
% sigma: Volatilities
% MU: Risk free rates
% m: Quantity of trials

St = Correlated_BM(S0,Tm,dt,A,sigma,MU,m);
UT = zeros(1,m); % Define option price matrix
for i=1:m
    pathi = St(:,i,:);
    pathi = reshape(pathi,(Tm/dt+1)*1,3);
    
    if max(pathi(:,2)) > B % condition (i)
        UT(i) = max(pathi(end,2)-K,0);
        continue;
    elseif max(pathi(:,2)) > max(pathi(:,3)) % condition (ii)
        UT(i) = max(pathi(end,2)-K,0)^2;
        continue;
    elseif mean(pathi(:,2)) > mean(pathi(:,3)) % condition (iii)
        UT(i) = max(mean(pathi(:,2))-K,0);
        continue;
    else % condition (iv)
        ST = mean(pathi(end,:));
        UT(i) = max(ST-K,0);
    end
end

Price = mean(UT)*exp(-mean(MU)*Tm);
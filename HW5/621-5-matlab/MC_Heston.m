function [Price, Bias, RMSE, Time]=MC_Heston(Scheme, S0, K, V0, Tm, r, sigma, kappa, theta, rho, n, m)
% For european call option
% S0: Initial stock price
% K: Strike price
% V0: Initial volatility
% Tm: Time to maturity
% r: Interest rate
% sigma: Volatility variance
% kappa: speed of mean-reversion of the variance
% theta: the long-term average variance
% rho: Correlation coefficient between two Brownian motions
% n: Time steps
% m: Quantity of trials

tic;
dt = Tm/n;
St = ones(1,m).*S0;
lnSt = log(St);
Vt1 = ones(1,m).*V0;
dWv = normrnd(0,1,n,m)*sqrt(dt);
dZ = normrnd(0,1,n,m)*sqrt(dt);
dWs = rho*dWv+sqrt(1-rho^2)*dZ;
    
for i=1:n
    if Scheme==1 || Scheme==4 || Scheme==5
        Vt2 = max(Vt1,0); % Effective variance
    else
        Vt2 = abs(Vt1);
    end
    
    lnSt = lnSt+(r-Vt2/2)*dt+sqrt(Vt2).*dWs(i,:);
    if Scheme==1
        f1=max(Vt1,0);
        f2=f1;
        f3=f1;
    elseif Scheme==2
        f1=abs(Vt1);
        f2=f1;
        f3=f1;
    elseif Scheme==3
        f1=Vt1;
        f2=f1;
        f3=abs(Vt1);
    elseif Scheme==4
        f1=Vt1;
        f2=f1;
        f3=max(Vt1,0);
    else
        f1=Vt1;
        f2=max(Vt1,0);
        f3=f2;
    end

    Vt1 = f1+kappa*(theta-f2)*dt+sigma.*f3.^(1/2).*dWv(i,:);
    


end
St = exp(lnSt);

Calls = max(St-K, 0);
Price = mean(Calls)*exp(-r*Tm);
Bias = abs(Price-6.8061);
RMSE = rms(Calls-6.8061);

Time = toc;
function [St]=Correlated_BM(S0,Tm,dt,A,sigma,MU,m)
% S0: Initial stock prices
% Tm: Time to maturity
% dt: Time intervel
% A: Correlation matrix
% sigma: Volatilities
% MU: Risk free rates
% m: Quantity of trials

L = chol(A,'lower'); % Cholesky factorization

n = Tm/dt;
nudt = (MU-0.5*sigma.^2)'*dt*ones(1,m);
St = zeros(n+1,m,3);
St(1,:,1)=S0(1);
St(1,:,2)=S0(2);
St(1,:,3)=S0(3);
lnSt = log(St);

for i=1:n
    Zt = normrnd(0,1,3,m)*sqrt(dt);
    Wt = L*Zt;
    for j=1:3
        lnSt(i+1,:,j)=lnSt(i,:,j)+nudt(j,:)+Wt(j,:).*sigma(j);
    end
end

St = exp(lnSt);

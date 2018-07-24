%% a

A=[1,0.5,0.2;0.5,1,-0.4;0.2,-0.4,1];
L = chol(A,'lower');

%% b

S0 = [100,101,98];
MU = [0.03,0.06,0.02];
sigma = [0.05,0.2,0.15];

Tm = 100/365;
m = 1000;
dt = 1/365;

St = Correlated_BM(S0,Tm,dt,A,sigma,MU,m);
Onepath = St(:,1,:);
Onepath = reshape(Onepath,(Tm/dt+1)*1,3);
plot(Onepath)
legend('Stock1','Stock2','Stock3')

%% c
iscall=1;K=100;m=1e4;
Basket_option(iscall,S0,K,Tm,dt,A,sigma,MU,m)

%% d
B=104;
Basket_exotic_option(S0,K,B,Tm,dt,A,sigma,MU,m)

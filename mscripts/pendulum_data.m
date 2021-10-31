clear all;
close all;

%test
A=[0 1; 0 0];
B=[0;1];
C=[1 1];

rho=0.1;
[Gl,Gn,P]=disobs(A,B,C);
Tc=[1 0; 1 1];
Ac=Tc*A*inv(Tc);
%[Ac,Bc,Cc,Tc] = obsfor(A,B,C);
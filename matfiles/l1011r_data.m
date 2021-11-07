% load the data from mat-file
clear all;
clc;

data=load('l1011r.mat');
A=data.A;
B=data.B;
C=data.C;

[n,m]=size(B);
[p,n]=size(C);

%% Chap. 5.7.1 Model-Reference System Using only Outputs

% design the swithcing surface
% lambda = [-0.05 -2 1.5 -1.5 1.5];
% nocomp = 1;
% specent=[  1 NaN NaN NaN NaN;
%          NaN   1   0 NaN NaN;
%          NaN   0   1 NaN NaN;
%          NaN NaN NaN   1 NaN;
%          0 NaN NaN NaN   1];
% 
% %[S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
% [S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
% S=inv(S*B)*S;

%S=rpp(A,B,lambda(1:n-m));

% Phi=diag([-5 -5]);

% [Lx,P,Lam]=contl(A,B,S,Phi);

% Design diagonal weighting matrix for the state vector
% Q=diag([5 1 1 5 5]); % compare with Eq. 4.139

% Designs a hyperplane S by minimising the linear quadratic cost function
% [S,E]=lqcf(A,B,Q);
% L=-inv(S*B)*(S*Am-Phi*S);

%[Acal,Dcal,Ccal,Lbar,D2,TL,r]=canon(A,B,C,-5);

Lx=data.Lx; % compare with Eq. 5.119

%[Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,-5,[-0.4470 -2.5061 -2.2372],[-5 -5]);

% Model-Reference matrices
Am=A+B*Lx;
eig(Am)

Cm=data.Cm;

Lr=-inv(Cm*inv(A+B*Lx)*B);

Bm=B*Lr;

% Design diagonal weighting matrix for the state vector
Q=diag([5 1 1 5 5]);

% Designs a hyperplane S by minimising the linear quadratic cost function
[S,E]=lqcf(A,B,Q);
L=-inv(S*B)*(S*Am-Phi*S);

% Design parameter for switching function 
rho=1;
delta=0.001;

% Initial state vector
x0=zeros(n);

%% compare calculated valued with original
fields = fieldnames(data);
for i=1:length(fields)
    if abs(data.(fields{i})-eval(fields{i})) > 1e-2
        error([fields{i},' is not equal'])
    end
end

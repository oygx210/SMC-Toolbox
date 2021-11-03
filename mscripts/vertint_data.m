clear all;

load vertint.mat;

%% Chap. 4.5.1 Model-Reference Design
lambda = [-5.6 4.2 -1];
nocomp=1;
specent=[  0   0   1;
           1 NaN   0;
         NaN   1   0;
         NaN NaN NaN;
         NaN NaN NaN];

% lambda = [-5.6 4.2 -1 -20 -20];
% specent=[  0   0   1 NaN NaN;
%            1 NaN   0 NaN NaN;
%          NaN   1   0 NaN NaN;
%          NaN NaN NaN   1 NaN;
%          NaN NaN NaN NaN   1];

[S,V]=dea(A,B,specent,lambda,nocomp);
S=inv(S*B)*S;

Phi=-20*eye(2);
[F,P,Lam]=contl(A,B,S,Phi);

H=[1 0 1 0 0;
   1 0 0 0 0];
G=-inv(H*inv(A+B*F)*B);

Am=A+B*F;
Bm=B*G;
Q=diag([10 5 5 20 20]);
[S,E]=lqcf(A,B,Q);
L=-inv(S*B)*(S*Am-Phi*S);

%[Af,Bf,Cf,Tcan,r]=outfor(A,B,C);

%
rho=1;
delta=0.001;
[n,m]=size(B);
[p,n]=size(C);

x0=zeros(n);

%% Chap. 4.5.2 Integral Action Based Design
[Aa,Ba]=intac(A,B,C);

nocomp=1;
specent=[NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN;
           0   0   1   0   1;
           1 NaN   0   1   0;
         NaN   1   0 NaN   0;
         NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN];

  
lambda = [-5.6 4.2 -1 -0.4 -0.7];

[S,V]=dea(Aa,Ba,specent,lambda,nocomp);
S=inv(S*Ba)*S;

Phi=-20*eye(2);

[L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

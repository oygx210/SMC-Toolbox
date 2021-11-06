% load the data from mat-file
clear all;
clc;

data=load('vertint.mat');
A=data.A;
B=data.B;
C=data.C;

[n,m]=size(B);
[p,n]=size(C);

%% Chap. 4.5.1 Model-Reference Design

% design the swithcing surface
lambda = [-5.6 4.2 -1 -20 -20];
nocomp = 1;
specent=[  0   0   1 NaN NaN;
           1 NaN   0 NaN NaN;
         NaN   1   0 NaN NaN;
         NaN NaN NaN   1 NaN;
         NaN NaN NaN NaN   1];

[S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
S=inv(S*B)*S;

Phi=diag(lambda(n-m+1:end));
[F,P,Lam]=contl(A,B,S,Phi);

H=[1 0 1 0 0;
   1 0 0 0 0]; % compare with Eq. 4.131
G=-inv(H*inv(A+B*F)*B); % compare with Eq. 4.133

% Model-Reference matrices
Am=A+B*F;
Bm=B*G;

% Design diagonal weighting matrix for the state vector
Q=diag([10 5 5 20 20]); % compare with Eq. 4.139

% Designs a hyperplane S by minimising the linear quadratic cost function
[S,E]=lqcf(A,B,Q); % compare with Eq. 4.140
L=-inv(S*B)*(S*Am-Phi*S);  % compare with Eq. 4.141

% Design parameter for switching function 
rho=1;
delta=0.001;

% Initial state vector
x0=zeros(n);

%% Chap. 4.5.2 Integral Action Based Design
lambda = [-5.6 4.2 -1 -0.4 -0.7];
specent=[NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN;
           0   0   1   0   1;
           1 NaN   0   1   0;
         NaN   1   0 NaN   0;
         NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN];

[bigA,bigB]=intac(A,B,C);

[S,V]=dea(bigA,bigB,specent,lambda,nocomp);
S=inv(S*bigB)*S;

[L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

%% compare calculated valued with original
fields = fieldnames(data);
for i=1:length(fields)
    val=eval(fields{i});
    if max(abs(data.(fields{i})(:)-val(:))) > 1e-1
        error([fields{i},' is not equal'])
    end
end

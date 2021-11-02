% Chapter 7.4: EXAMPLE: A Temperature Control Scheme
% Sub-Chapter 7.4.2: Controller Design
% observer based controller design with integral action

clear all;
%close all;

load vertint.mat;

% load aircraft.mat;
% 
[Aa,Ba]=intac(A,B,C);

% 
nocomp=1;
% lambda=[-5.6+-i*4.2 -1 -0.4 -0.7];
lambda=[-5.6 4.2 -1 -0.4 -0.7];
% 4 Eigenvectors (only the first one is complex)
%        1re 1im   2   3   4
specent=[NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN;
           0   0   1   0   1;
           1 NaN   0   1   0;
         NaN   1   0 NaN   0;
         NaN NaN NaN NaN NaN;
         NaN NaN NaN NaN NaN];

%G=-inv(C*inv(A+B*F)*B);
[G,~]=wzobs(A,B,C,[-4 -4.425 -4.5 -5 -4]);

%specent=eigbuil(7,2,nocomp)

[S,V]=dea(Aa,Ba,specent,lambda,nocomp);

S=inv(S*Ba)*S;

Phi=-20*eye(size(Ba,2));

[L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

Gamma=[-0.9 0;
       0 -0.7];

%%
rho=1;
delta=0.001;
[n,m]=size(B);
[p,n]=size(C);

z0=[0 0 0 0 0];

% A=[0,1,0;
%    0,0,1;
%    -0.0001, -0.0082, -0.1029];
% B=[0;0;1];
% C=[0.0001,0.0022,0.0053];
% 
% Gamma=-0.025;
% 
% Phi=-0.1;
% S=[-0.5372,0.0019,0.0822,1.000];
% 
% % design compensator with reduced observer
% 
% % with given Sr
% % Sr=26.1642;
% % [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi,Sr);
% 
% % without given Sr
% [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);
% 
% %%
% open('obsint_mdl.slx');
% sim('obsint_mdl');
% 
% plot(t.Data,er.Data);
% grid on;
% legend({'er'});
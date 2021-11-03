%% clear the workspace
clear all;

%% Load the system data
%load compped.mat;
%load invpen.mat;
%load rpv.mat;

load dcmotor.mat;

Af=A; Bf=B; Cf=C;
[n,m]=size(B);
[p,n]=size(C);

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
% x0(2)=1;
% x0(3)=1;

%% Design of the compensator
rho=0.5; % user defined parameter

% hyperlane design with desired poles
wn=2.0;
zeta=0.95;
sp=roots([1 2*zeta*wn wn^2]); % second order
S=rpp(A,B,sp);

Phi=-5;
[L,P,Lam]=contl(A,B,S,Phi);

% check roots of closed loop
% The condition number of the eigenvectors is an indicator of the
% robustness with respect to unmatched uncertainty (smaller => better)
M=inv(S*B)*S;
[V,D] = eig(A-B*M);
cond(V)

%% Simulate the model
mdl_name='streg_mdl';
if ~bdIsLoaded(mdl_name)
    open('streg_mdl.slx');
end
sim(mdl_name);

%% Plot the results
subplot(3,1,1)
plot(t.Data,u.Data);
grid on;
legend({'u'});
subplot(3,1,2)
plot(t.Data,x.Data);
grid on;
legend({'x'});
subplot(3,1,3)
plot(t.Data,s.Data);
legend({'s'});
grid on;
legend({'s'});

%plot(y.Time,y.Data(:,3),yobs.Time,yobs.Data(:,3));
%% clear the workspace
clear all;

%% Load the system data
%load aircraft.mat;
%load compped.mat;
%load dcmotor.mat;
%load furnobs.mat;
%load hecka.mat;
%load helcopt.mat;
%load huiex.mat;
%load invpen.mat;
%load l1011r.mat;
%load oxymodel.mat;
load vertint.mat;

[n,m]=size(B);
[p,n]=size(C);

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
% x0(2)=1;
% x0(3)=1;

%% Design of the compensator
rho=0.5; % user defined parameter

% hyperlane design with desired poles for system dynamics
wn=2.0; %  [rad/s] cut-off angular frequency
zeta=0.95; % [-] damping ratio less 1 --> complex poles
sp=roots([1 2*zeta*wn wn^2]); % second order for first 2 poles

% additional poles (real) by system order greater than 3
if n-m>2
    % to avoid multiple poles iterate through and increase each value
    for i=1:(n-m-2)
        sp=[sp; -i];
    end
end

% Calculate the switching surface using robust eigenstructure assignment
S=rpp(A,B,sp);

% Calculate L and Lam required for the linear component of the control law
Phi=-5*eye(size(S,1));
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
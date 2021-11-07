% Chapter 5: Sliding Mode Controllers Using Output Information

%% clear the workspace
clear all;

%% Load the system data
% mat_name = 'aircraft';
mat_name = 'compped';  % Sub-Chapter 5.6.2: Design Example 1
mat_name = 'invpen';  % Chap. 5.6.3 Design Example 2: Inverted Pendulum
% mat_name = 'dcmotor';
% mat_name = 'furnobs';
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';
% mat_name = 'l1011r';
%mat_name = 'vertint';
load([mat_name,'.mat']);

[n,m]=size(B);
[p,n]=size(C);

switch mat_name
    case 'compped'
        % define desired dynamics
        ped=-2.5;       % poles for error dynamics
        psm=[-1 -1.5];  % poles for sliding motion
        prsd=-5;        % poles for range space dynamics

        x0=[0 1 0 0]; % [xr, x11, x12]
        xc0=[0 0];

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'invpen'
        % design of sliding mode poles via LQR
        Q=diag([10 1 1 0.1]);
        [~,E]=lqcf(A,B,Q);
        
        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        psm = E; % reduced order sliding motion poles
        prsd = -6; % desired pole(s) for range space dynamics

        x0=[0 0.1 0 0];
        xc0=-0.9855;

        rho=1;
        delta=0.001;

        SimStopTime=5;
end

Lo=ped;

% design compensator with reduced observer
[Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);


%% simulate the model
open('comp_mdl.slx');
set_param('comp_mdl','StopTime',sprintf('%f',SimStopTime));
sim('comp_mdl');

%% plot simulation results
subplot(2,2,1)
%plot(t.Data,x.Data(:,3),t.Data,xm.Data(:,3),t.Data,xobs.Data(:,3));
plot(t.Data,s.Data);
grid on;
legend('switch_fcn','Interpreter','none');
subplot(2,2,3)
plot(t.Data,x.Data);
grid on;
legend({'x1','x2','x3','x4'});

xc=xhat.Data(:,2);
x11=x.Data(:,2);
x12=x.Data(:,3);
ec=xc-x11-Lo*x12;
subplot(2,2,2)
plot(t.Data,ec);
grid on;
legend({'ec'});

xr=x.Data(:,1);
zr=xhat.Data(:,1);
er=zr-xr;
subplot(2,2,4)
plot(t.Data,er);
grid on;
legend({'er'});


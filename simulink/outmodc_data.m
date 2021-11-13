%% clear the workspace
clear all;

%% Load the system data
%mat_name = 'aircraft';
%mat_name = 'compped';
%mat_name = 'dcmotor';
%mat_name = 'furnobs';
%mat_name = 'hecka';
%mat_name = 'helcopt';
%mat_name = 'huiex';
%mat_name = 'invpen';
mat_name = 'l1011r';
%mat_name = 'oxymodel';
%mat_name = 'vertint';

load([mat_name,'.mat']);

%% Get system dimensions
[n,m]=size(B);
[p,n]=size(C);

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
% x0(2)=1;
% x0(3)=1;

%% Design of the compensator
% Chap. 4.5.1 Model-Reference Design

switch mat_name
    case 'l1011r'
        % Design diagonal weighting matrix for the state vector
        Q=diag([5 1 1 5 5]);

        % design of sliding mode poles via LQR
        [S,E]=lqcf(A,B,Q);

        % Design the dynamic compensator
        % ped = -5; %  desired pole(s) for error dynamics (Lo)
        ped=[-4,-4.425,-4.5,-5];    % poles for error dynamics
        %psm=[-0.4470 -2,2372 -2.5061];  % reduced order sliding motion poles
        %psm = E';
        psm=-3;                     % poles for sliding motion
        prsd = -2.0*ones(1,m); % desired pole(s) for range space dynamics

        Am=A+B*Lx;
        Bm=B*-inv(Cm*inv(Am)*B);
    
        [G,F]=wzobs(A,B,C,psm,ped);

        %ped=[-4,-4.425,-4.5,-5];    % poles for error dynamics

        [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped(1),E',prsd);

        %P=0.25*eye(m);

        x0=zeros(1,n);
        x0(1)=1;

        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=10;
end

%% Simulate the model
mdl_name='outmodc_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% Plot the results
subplot(4,1,1)
plot(t.Data,u.Data);
grid on;
legend({'u'});
subplot(4,1,2)
plot(t.Data,x.Data);
grid on;
legend({'x'});
subplot(4,1,3)
plot(t.Data,xm.Data);
grid on;
legend({'xm'});
subplot(4,1,4)
plot(t.Data,s.Data);
legend({'s'});
grid on;
legend({'s'});

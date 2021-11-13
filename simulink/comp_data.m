%% Compensator simulation
% Demonstration of using the function "comrobs"

%% clear the workspace
clear all;

%% Load the system data
mat_name = 'compped';  % Chapter 5.6.2: Design Example 1
mat_name = 'invpen';  % Chap. 5.6.3 Design Example 2: Inverted Pendulum
mat_name = 'l1011r';  % Chap. 5.7.1 Aircraft Example
mat_name = 'dcmotor'; % Chap. 3.6.4 Example: Control of a DC Motor

% mat_name = 'vertint'; % Chap 4.5 Design Study: Pitch-Pointing Flight Controller
% mat_name = 'furnobs'; % Chap. 7.4: EXAMPLE: A Temperature Control Scheme
% mat_name = 'hecka';
% mat_name = 'helcopt';
% mat_name = 'huiex';

load([mat_name,'.mat']);

[n,m]=size(B); % n-state vector length, m-number of inputs
[p,n]=size(C); % n-state vector length, p-number of outputs

switch mat_name
    case 'compped'
        % define desired dynamics
        ped=-2.5;       % poles for error dynamics
        psm=[-1 -1.5];  % poles for sliding motion
        prsd=-5*ones(1,m);        % poles for range space dynamics

        x0=zeros(1,n); % [xr, x11, x12]
        x0(2)=1;
        xc0=zeros(1,n-p);

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
        prsd = -6*ones(1,m); % desired pole(s) for range space dynamics

        x0=zeros(1,n);
        x0(2)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=-0.9855;

        rho=1;
        delta=0.001;

        SimStopTime=5;

    case 'l1011r'
        % Design diagonal weighting matrix for the state vector
        Q=diag([5 1 1 5 5]);

        % design of sliding mode poles via LQR
        [~,E]=lqcf(A,B,Q);

        % Design the dynamic compensator
        ped = -5; %  desired pole(s) for error dynamics (Lo)
        psm = E'; % reduced order sliding motion poles
        prsd = -5*ones(1,m); % desired pole(s) for range space dynamics

        x0=zeros(1,n);
        x0(1)=1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=10;

    case 'dcmotor'
        C=C([1,3],:); % partial output matrix to allow a reduced observer design
        [p,n]=size(C);

        % Design diagonal weighting matrix for the state vector
        Q=diag([5 5 1]);

        % design of sliding mode poles via LQR
        [~,E]=lqcf(A,B,Q);

        % Design the dynamic compensator
        ped = -10; %  desired pole(s) for error dynamics (Lo)
        psm = E'; % reduced order sliding motion poles
        prsd = -20*ones(1,m); % desired pole(s) for range space dynamics

        x0=zeros(1,n);
        x0(1)=0.1;
        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=10;
end

Lo=ped;

% design compensator with reduced observer
[Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);


%% simulate the model
mdl_name='comp_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

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


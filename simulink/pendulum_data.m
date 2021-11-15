%% Chap. 1.5: Double Integrator - State-space Approach
%% Chap. 6.4: Robust Observer Design

clear all;
close all;

%test
A=[0 1; 0 0];
B=[0;1];
C=[1 1];

[n,m]=size(B); % n-state vector length, m-number of inputs
[p,n]=size(C); % n-state vector length, p-number of outputs

% observer design
psm = [-1 -1];
per = -1;

rho=0.1;
[Gl,Gn,P2]=disobs(A,B,C,psm,per);
Tc=[1 0; 1 1];
Ac=Tc*A*inv(Tc);
Bc=Tc*B;
Cc=C*inv(Tc);

% alternative calculation of canonical form
[Ac,Bc,Cc,Tc] = obsfor(A,B,C);


% recompute the Lyapunov matrix with modifid weighting
Q2=2;
a22s=diag(per,0);
P2=lyap(a22s',diag(Q2));



x0=zeros(1,n);
x0(1)=1.5;

rho=1.0;
delta=0.01;

SimStopTime=20;

%% simulate the model
mdl_name='pendulum_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% plot simulation results
subplot(2,2,1)
plot(t.Data,s.Data);
grid on;
legend(get_legend('s'));

xc=xhat.Data(:,2);
x11=x.Data(:,2);
x12=x.Data(:,3);
ec=xc-x11-Lo*x12;

xr=x.Data(:,1);
zr=xhat.Data(:,1);
er=zr-xr;

subplot(2,2,2)
plot(t.Data,ec,t.Data,er);
grid on;
legend([get_legend('ec');get_legend('er')]);

subplot(2,2,3)
plot(t.Data,x.Data);
grid on;
legend(get_legend('x'));

subplot(2,2,4)
plot(t.Data,xhat.Data);
grid on;
legend(get_legend('xhat'));

sgtitle([mdl_name,' - ',mat_name],'Interpreter','None');
% Chapter 5: Sliding Mode Controllers Using Output Information
% Sub-Chapter 5.6.2: Design Example 1

% load system matrices A,B,C
load compped.mat;

% define desired dynamics
ped=-2.5;       % poles for error dynamics
psm=[-1 -1.5];  % poles for sliding motion
prsd=-5;        % poles for range space dynamics

Lo=ped;

% design compensator with reduced observer
[Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);

%%
rho=1;
delta=0.001;
[n,m]=size(B);
[p,n]=size(C);

x0=[0 1 0 0]; % [xr, x11, x12]
xc0=[0 0];

%%
open('comp_mdl.slx');
sim('comp_mdl');

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


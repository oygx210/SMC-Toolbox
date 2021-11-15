%% Chap. 6.4: Robust Observer Design for pendulum

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

[Gl,Gn,P2]=disobs(A,B,C,psm,per);
% Tc=[1 0; 1 1];
% Ac=Tc*A*inv(Tc);
% Bc=Tc*B;
% Cc=C*inv(Tc);

% alternative calculation of canonical form
[Ac,Bc,Cc,Tc] = obsfor(A,B,C);


% recompute the Lyapunov matrix with modified weighting
% Q2=2;
% a22s=diag(per,0);
% P2=lyap(a22s',diag(Q2));

x0=zeros(1,n);
x0(1)=1.5;

z0=zeros(1,n);

% design parameter for the switching function
rho=1.0;
delta=0.01;

SimStopTime=10;

%% simulate the model
mdl_name='pendulum_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% plot simulation results
subplot(2,3,1)
plot(t.Data,u.Data);
grid on;
legend(get_legend('u'));

subplot(2,3,2)
plot(t.Data,ul.Data,t.Data,un.Data);
grid on;
legend([get_legend('ul');get_legend('un')]);

subplot(2,3,3)
plot(t.Data,ey.Data);
grid on;
legend(get_legend('ey'));

subplot(2,3,4)
plot(t.Data,y.Data,t.Data,yo.Data);
grid on;
legend([get_legend('y');get_legend('yo')]);

subplot(2,3,5)
plot(t.Data,x.Data(:,1),t.Data,z.Data(:,1));
grid on;
legend([get_legend('x',1);get_legend('z',1)]);

subplot(2,3,6)
plot(t.Data,x.Data(:,2),t.Data,z.Data(:,2));
grid on;
legend([get_legend('x',2);get_legend('z',2)]);

sgtitle(mdl_name,'Interpreter','None');
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
        prsd = -5*ones(1,m); % desired pole(s) for range space dynamics

        %[Ac,Dc,Cc,Tc]=obsfor(A,B,C,psm);
        %[Acal,Bcal,Ccal,L,B2,TL,r]=canon(A,B,C,psm);

        [L,P,~]=contl(A,B,S,-2.0*eye(m));

        Am=A+B*Lx;
        Bm=B*-inv(Cm*inv(Am)*B);
    
        [G,F]=wzobs(A,B,C,psm,ped);

        %P=0.25*eye(m);

        x0=zeros(1,n);
        x0(1)=1;

        z0=zeros(1,n);
        z0(3)=0.01;

        xc0=zeros(1,n-p);
        xc0(1)=0;

        rho=1;
        delta=0.001;

        SimStopTime=10;
end

%% Simulate the model
mdl_name='obsmod_mdl';
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


%%
% % Kap 5.7.1 and 7.5.1
% clear all;
% %close all;
% 
% load l1011r.mat;
% 
% [V,D]=eigenshuffle(A+B*Lx);
% [n,m]=size(B);
% [p,n]=size(C);
% 
% Q=diag([5 1 1 5 5]);
% [S,E]=lqcf(A,B,Q);


% [Ao,Bo,Co,Tc,r]=outfor(A,B,C);
% Cot=qr(Co);
% 
% %lambda=[-0.5 -1.2 -6];
% lambda=[-0.05 -2+1.5j -2-1.5j -1.5+1.5j -1.5-1.5j];
% Vo=zeros(size(B,1),size(B,2),length(lambda));
% Voc=zeros(size(B,2),size(B,1),length(lambda));
% 
% for i=1:length(lambda)
%     Votemp=-inv(lambda(i)*eye(n)-Ao)*Bo;
%     Vo(:,:,i)=Votemp;
%     Voc(:,:,i)=inv(Votemp'*Votemp)*Votemp';
% end
% 
% cot= Cot(2,:);
% cobt=(Cot(3,:)-1)*-1;
% 
% ro(:,1)=Voc(:,:,1)*cot';
% ro(:,2)=Voc(:,:,2)*cot';
% ro(:,3)=Voc(:,:,3)*cot';
% ro(:,4)=Voc(:,:,4)*cot';
% ro(:,5)=Voc(:,:,5)*cobt';
% 
% not(1,:)=Vo(:,:,1)*ro(:,1);
% not(2,:)=Vo(:,:,2)*ro(:,2);
% not(3,:)=Vo(:,:,3)*ro(:,3);
% not(4,:)=Vo(:,:,4)*ro(:,4);
% not(5,:)=Vo(:,:,5)*ro(:,5);
% 
% Po=not';
% Ro=ro;
% 
% Qo=[Po;Ro];
% 
% Ko=Ro*inv(Po);
% 
% K=Ko*inv(Tc);
% 
% % check roots of closed loop
% %[V,D] = eig(A-B*K);
% [V,D] = eigenshuffle(A-B*K);





% open('obsmod_mdl.slx');
% 
% sim('obsmod_mdl');
% 
% %%
% subplot(3,1,1)
% %plot(t.Data,u.Data);
% plot(t.Data(1:1e4),p_obs_err.Data(1:1e4,:));
% grid on;
% legend({'error'});
% subplot(3,1,2)
% %plot(t.Data,x.Data(:,3),t.Data,xm.Data(:,3),t.Data,xobs.Data(:,3));
% plot(t.Data(1:1e4),switch_fcn.Data(1:1e4,:));
% grid on;
% %legend({'x','xm','xobs'});
% legend('switch_fcn');
% subplot(3,1,3)
% plot(t.Data,yobs.Data);
% grid on;
% legend({'wasched-out yaw rate','roll rate','sideship angle','bank angle'});
% 
% %plot(y.Time,y.Data(:,3),yobs.Time,yobs.Data(:,3));
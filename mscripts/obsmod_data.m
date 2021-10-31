% Kap 5.7.1 and 7.5.1
clear all;
%close all;

load l1011r.mat;

[V,D]=eigenshuffle(A+B*Lx);
[n,m]=size(B);
[p,n]=size(C);

Q=diag([5 1 1 5 5]);
[S,E]=lqcf(A,B,Q);


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



Am=A+B*Lx;
Bm=B*-inv(Cm*inv(Am)*B);
%psm=[-0.4470 -2.2372 -2.5061];

psm=-3;                     % poles for sliding motion
ped=[-4,-4.425,-4.5,-5];    % poles for error dynamics

[Acal,Bcal,Ccal,L,B2,TL,r]=canon(A,B,C,psm);


[G,F]=wzobs(A,B,C,psm,ped);

rho=1;

x0=[0 0 0 0 0];
z0=[0 0 0.1 0 0];
P=0.25*eye(2);

open('obsmod_mdl.slx');

sim('obsmod_mdl');

%%
subplot(3,1,1)
%plot(t.Data,u.Data);
plot(t.Data(1:1e4),p_obs_err.Data(1:1e4,:));
grid on;
legend({'error'});
subplot(3,1,2)
%plot(t.Data,x.Data(:,3),t.Data,xm.Data(:,3),t.Data,xobs.Data(:,3));
plot(t.Data(1:1e4),switch_fcn.Data(1:1e4,:));
grid on;
%legend({'x','xm','xobs'});
legend('switch_fcn');
subplot(3,1,3)
plot(t.Data,yobs.Data);
grid on;
legend({'wasched-out yaw rate','roll rate','sideship angle','bank angle'});

%plot(y.Time,y.Data(:,3),yobs.Time,yobs.Data(:,3));
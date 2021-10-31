%% State Feedback Control Design Using Eigenstructure Decoupling
% Robert Fonod, Pavol Kocsis
% HAL Id: hal-00905233
% https://hal.inria.fr/hal-00905233
% Submitted on 17 Nov 2013
% Additional links:
% http://users.isy.liu.se/en/rt/johans/files/articles/andry_etal_83.pdf
clear all;
A=[ 0  1  0;
    0  0  1;
   -5 -9 -5];
B=[1 3;
   2 1;
   2 5];
C=[1 2 1;
   1 1 0];

[n,m]=size(B);
[p,n]=size(C);

Ao=[-1  10.5  6;
     0 -3.0  -2;
     0  1.0  -1];
 
Bo=[0  0;
    7 10;
    3  4];

Tc=[ 1.0  2.5 -5.5;
    -1.0 -2.5  6.5;
     1.0  3.5 -7.5];
 
Co=C*Tc;
%[Ao,Bo,Co,Tc,r]=outfor(A,B,C);
Cot=qr(Co);

lambda=[-0.5 -1.2 -6];
Vo=zeros(size(B,1),size(B,2),length(lambda));
Voc=zeros(size(B,2),size(B,1),length(lambda));

for i=1:length(lambda)
    Votemp=-inv(lambda(i)*eye(n)-Ao)*Bo;
    Vo(:,:,i)=Votemp;
    Voc(:,:,i)=inv(Votemp'*Votemp)*Votemp';
end

cot= Cot(2,:);
cobt=(Cot(1,:)-1)*-1;

ro(:,1)=Voc(:,:,1)*cot';
ro(:,2)=Voc(:,:,2)*cot';
ro(:,3)=Voc(:,:,3)*cobt';

not(1,:)=Vo(:,:,1)*ro(:,1);
not(2,:)=Vo(:,:,2)*ro(:,2);
not(3,:)=Vo(:,:,3)*ro(:,3);

Po=not';
Ro=ro;

Qo=[Po;Ro];

Ko=Ro*inv(Po);

K=Ko*inv(Tc);

% specent=[  1;
%            0;
%            0];
% 
% [S,V]=dea(Ao,Bo,specent,-0.5);

% check roots of closed loop
%[V,D] = eig(A-B*K);
[V,D] = eigenshuffle(A-B*K);


% load the data from mat-file
clear all;
clc;

data=load('invpen.mat');
A=data.A;
B=data.B;
C=data.C;

[n,m]=size(B);
[p,n]=size(C);

%% Chap. 5.6 Dynamic Compensation
% Chap. 5.6.3 Design Example 2: Inverted Pendulum

[Af,Bf,Cf,Tcan,r]=outfor(A,B,C);

if r==0
    disp('the system (A,B,C) does not have invariant zeros');
end

% design of sliding mode poles via LQR
Q=diag([10 1 1 0.1]);
[~,E]=lqcf(A,B,Q);

% Design the dynamic compensator
ped = -10; %  desired pole(s) for error dynamics (Lo)
%psm = [-4.3241+1i*1.7852,-4.3241-1i*1.7852,-3.1623]; % reduced order sliding motion poles
psm = E'; % reduced order sliding motion poles
prsd = -6; % desired pole(s) for range space dynamics

[Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd);

%% compare calculated valued with original
fields = fieldnames(data);
for i=1:length(fields)
    val=eval(fields{i});
    if max(abs(data.(fields{i})(:)-val(:))) > 1e-1
        error([fields{i},' is not equal'])
    end
end

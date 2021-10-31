% Chapter 7.4: EXAMPLE: A Temperature Control Scheme
% Sub-Chapter 7.4.2: Controller Design
% observer based controller design with integral action

clear all;
%close all;

load furnobs.mat;

A=[0,1,0;
    0,0,1;
    -0.0001, -0.0082, -0.1029];
B=[0;0;1];
C=[0.0001,0.0022,0.0053];

Gamma=-0.025;

Phi=-0.1;
S=[-0.5372,0.0019,0.0822,1.000];

% design compensator with reduced observer

% with given Sr
% Sr=26.1642;
% [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi,Sr);

% without given Sr
[L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

%%
open('obsint_mdl.slx');
sim('obsint_mdl');
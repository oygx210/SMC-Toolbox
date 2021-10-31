%%
clc
close all
clear all
%% Load Data and Retrieve State Space Matricies

%% BWB Example
load BWB
A=BWB.lat.A
B=BWB.lat.B
C=eye(4);


D=zeros(4,8);
x0=[0.1 0 0 0]'
Q=eye(4);




Cc=[1 0 0 0
    0 1 0 0];
Qi=eye(2);

F=-eye(2);

%  [L, P, Lam, S, E, CA]=SMCDesign_init(A, B, Q)
%  [L, Lr, P, Lam, S, E, CA]=SMCDesign_initIA(A, B, Cc, eye(6))
open Example
 


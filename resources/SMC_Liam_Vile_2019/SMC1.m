%[L, R S, E]=SMCDesign(A, B, C, Q) 
A=BWB.A
B=BWB.B
%% Calculate sizes

[nn, mm]= size(B)
%% Create Regular Form 
%
% Regular form:
%
% A=|A_11 A_12|     B=| 0 |
%   |A_21 A_22|       |B_2|
%
[A11,A12,B2,Tr]=regfor(A,B);
%% 


function [S,E]=lqcf(A,B,Q)

% [S,E]=lqcf(A,B,Q)
%
%       Linear quadratic cost function designs a hyperplane S by minimising 
%       the cost function. The problem is of the form described below, where Q
%       is symmetric and positive definite
%       
%               J = Integral x'Qx  dt           (1)
%
%       This program minimises the equation above subject to the equation: 
%               .
%               z1(t) = A11z1(t) + A12z2(t) 
%                
%       The matrix Q is transformed and partitioned compatibly with z so that
%
%       Tr'*Q*Tr=[Q11 Q12] 
%                [Q21 Q22]
%
%       Equation (1) can then be expressed as
%
%       J = Integral z1'*Q11z1 + 2z1'Q12z2 + z2'Q22z2 dt
% 
%       The function returns the hyperplane matrix S and the eigenvalues of
%       the sliding motion in the vector E.

	
%       Chris Edwards, Robert Cortez & Sarah Spurgeon
%       Control Systems Research
%       Leicester University
%       University Road
%       Leicester LE1 7RH
%
%       Email: ce@sun.engg.le.ac.uk
%
%
%       Version 1.0
%       29/11/97
%

% Establish the size of Q
[nnx,nny]=size(Q);
[nn,mm]=size(B);

%----------------------------------------------------------------------------%
% Check whether the plant is controllable, and that the matrices
% conform to various conditions eg size etc...
%----------------------------------------------------------------------------%



%----------------------------------------------------------------------------%
% Check if Q is the right size and positive definite symmetric
%----------------------------------------------------------------------------%
if nn~=nnx | nn~=nny
   error('The matrix Q is not consistent with the state dimension');end

nnx = norm(Q,1);
if any(eig(Q) <= eps*nnx) | (norm(Q'-Q,1)/nnx > eps)
   error('The matrix Q is not symmetric and positive definite');end

%----------------------------------------------------------------------------%
% Obtain the regular form transformation using the function 
%----------------------------------------------------------------------------%

% [A11,A12,B2,Tr]=regfor(A,B);
[A, Bv, Tr, CA]=SMCCanForm(A, B);
[A11, A12]=partiton(A, Bv);
%----------------------------------------------------------------------------%
% Transform weighting matrix to regular form coordinate system
%----------------------------------------------------------------------------%
Qt=Tr*Q*Tr';

% Compatibly partition with regular form description
Q11 = Qt(1:nn-mm,1:nn-mm);
Q12 = Qt(1:nn-mm,nn-mm+1:nn);
Q21 = Qt(nn-mm+1:nn,1:nn-mm);
Q22 = Qt(nn-mm+1:nn,nn-mm+1:nn);

% Form reduced order system description and associated weighting matrix
Qhat=Q11-Q12*inv(Q22)*Q21;
Ahat=A11-A12*inv(Q22)*Q21;

% Solve the LQR problem
[K,P1,E]=lqr(Ahat,A12,Qhat,Q22);

% Obtain the switching function matrix in terms of the original coordinates
M=inv(Q22)*(A12'*P1+Q21);

S=[M eye(mm)]*Tr;

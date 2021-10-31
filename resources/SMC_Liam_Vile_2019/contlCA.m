function [L,P,Lam]=contlCA(A,B,S,Phi)

% [L,P,Lam]=contl(A,B,S,Phi) 
%
%          Returns the variables L and Lam required for the linear component  
%          of the control law along with the Lyapunov matrix which satisfies 
%          P2*Phi+Phi'*P2=-I, where Phi is a stable design matrix required for 
%          the range space dynamics.     

%          Chris Edwards, Robert Cortez & Sarah Spurgeon
%          Control Systems Research
%          Leicester University
%          University Road
%          Leicester LE1 7RH
%
%          Email: ce@sun.engg.le.ac.uk
%
%
%          Version 1.0
%          29/11/97
%
%

[nn,mm]=size(B);
[mx,my]=size(Phi);

%----------------------------------------------------------------------------%
% Check the size of Phi and whether it is a stable design matrix 
%----------------------------------------------------------------------------%
if mx~=mm | my~=mm
   error(' The size of the design matrix Phi is inconsistent');end
 
if any(real(eig(Phi))>=-eps)
      warndlg(' The design matrix is unstable');end

%----------------------------------------------------------------------------%
% Calculate the linear feedback gain, the Lyapunov matrix and Lam
%----------------------------------------------------------------------------%
Lam=S*B;
L=-inv(S*B)*((S*A)-(Phi*S));
P = lyap(Phi',eye(mm));

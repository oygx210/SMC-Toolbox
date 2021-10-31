function [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi,Sr_in)
%
% [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi,Sr)
%
%         Returns the controller parameters in the case where integral action
%         is used for tracking purposes. The matrix S represents the state
%         feedback component of the hyperplane and Sr represents the reference
%         dependant component. Lam is the gain associated with the unit vector
%         P is the Lyapunov matrix which scales the hyperplane in the control
%         law, L represents the linear state feedback gain and Lr and Ldot
%         represent feedforward gains for the reference and derivative of the
%         reference respectively.
%
%         The file can also be used in the from
%         [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi)
%         in which case the reference dependant component of the hyperplane 
%         Sr is computed to ensure, that in the nominal case, the steady 
%         state values of the integral action states tend to zero.

%         Chris Edwards, Robert Cortez & Sarah Spurgeon
%         Control Systems Research
%         Leicester University
%         University Road
%         Leicester LE1 7RH
%
%         Email: ce@sun.engg.le.ac.uk
%
%         Version 2
%         20/8/98
%
%

[pp,nn]=size(C);
[nn,mm]=size(B);

%----------------------------------------------------------------------------%
% Check the size of Phi and whether it is a stable design matrix 
%----------------------------------------------------------------------------%
[mx,my]=size(Phi);
if mx~=mm || my~=mm
   error('The size of the design matrix Phi is inconsistent');end
 
if any(real(eig(Phi))>=-eps)
      warndlg(' The design matrix is unstable');end


%----------------------------------------------------------------------------%
% Check the size of Sr  
%----------------------------------------------------------------------------%
if nargin>5
  [mx,my]=size(Sr_in);
  if (mx~=mm || my~=pp)
     error('The size of the design matrix Sr is inconsistent');
  end
end
 
%----------------------------------------------------------------------------%
% Augment the statespace with integral action states
%----------------------------------------------------------------------------%
AT=[zeros(pp,pp) -C; zeros(nn,pp) A];
[A11,A12,A21,A22,B2,Tr]=regfor(A,B);
BT=[zeros(pp,mm); B];

%----------------------------------------------------------------------------%
% Change coordinates so that regular form is achieved
%----------------------------------------------------------------------------%
Trtilde=[eye(pp) zeros(pp,nn);zeros(nn,pp) Tr];
ATnew=Trtilde*AT*Trtilde';
A11a=ATnew(1:nn,1:nn);
A12a=ATnew(1:nn,nn+1:nn+pp);
Snew=S*Trtilde';

%----------------------------------------------------------------------------%
% Compute the associated feedback gain 
%----------------------------------------------------------------------------%

M=inv(Snew(:,nn+pp-mm+1:nn+pp))*Snew(:,1:nn+pp-mm);
A11s=A11a-A12a*M;

%----------------------------------------------------------------------------%
% Compute the design parameter Sr
%----------------------------------------------------------------------------%
Br=[eye(pp); zeros(nn-mm,pp)];
if nargin==5
    Sr=-Snew(:,nn+pp-mm+1:nn+pp)*inv(Br'*inv(A11s)*A12a)*Br'*inv(A11s)*Br;
else if nargin>5
        Sr=Sr_in;
end
end

%----------------------------------------------------------------------------%
% Compute the controller parameters
%----------------------------------------------------------------------------%
Lam=S*BT;
P=lyap(Phi',eye(mm));
L=-inv(Lam)*(S*AT-Phi*S);
Lr=-inv(Lam)*(Phi*Sr+Snew(:,1:pp));
Lrdot=inv(Lam)*Sr;
   

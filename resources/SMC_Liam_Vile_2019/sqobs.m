function [G,F]=sqobs(A,B,C)

% [G,F]=sqobs(A,B,C)
%
%        This function produces an observer for the special case of square
%        systems. The arguments returned are the linear gain matrix G, along
%        with F. During the execution the user is prompted for a pth order
%        vector to contain stable real poles (where p is the number of outputs).
 

%        Chris Edwards, Robert Cortez & Sarah Spurgeon
%        Control Systems Research
%        Leicester University
%        University Road
%        Leicester LE1 7RH
%
%        Email: ce@sun.engg.le.ac.uk 
%
%
%        Version 2.0
%        7/1/98
%

[nn,mm]=size(B);
[pp,nn]=size(C);                

G=[];
F=[];
%----------------------------------------------------------------------------%
% This function checks the consistency of the matrix sizes etc...
% along with the conditions that B,C and CB are full rank.
%----------------------------------------------------------------------------%
msg=abcchk(A,B,C);
if ~isempty(msg)
   error(msg);
end

% Check if the system is square
if pp~=mm
   error('System must be square');
end

if rank(C*B)~=pp
   error('The Markov parameter CB is rank deficient')
end

%---------------------------------------------------------------------------%
% First put the triple (A,B,C) in regular form
%---------------------------------------------------------------------------%
[A11,A12,B2,Tr]=regfor(A,B);
Ar=Tr*A*Tr';
Br=Tr*B;
Cr=C*Tr';

%--------------------------------------------------------------------------%
% Partition the system matrix to use in the observer gain matrix
%--------------------------------------------------------------------------%
A11=Ar(1:nn-pp,1:nn-pp);
A21=Ar(nn-pp+1:nn,1:nn-pp);
A12=Ar(1:nn-pp,nn-pp+1:nn);
A22=Ar(nn-pp+1:nn,nn-pp+1:nn);
B2=Br(nn-pp+1:nn,:);
C1=Cr(:,1:nn-pp);
C2=Cr(:,nn-pp+1:nn);

%--------------------------------------------------------------------------%
% Check if the invariant zeros are stable
%--------------------------------------------------------------------------%
if any(real(eig(A11-A12*inv(C2)*C1))>=-100*eps)
   error('The invariant zeros are not stable');
end

msg=' ';
pmsg=['Enter ',num2str(pp),' real stable poles for the output error system '];
while msg~=[]
  disp(msg)
  p=input(pmsg);
  p=p(:);
  msg=polechk(p,pp,1);
end
   
%--------------------------------------------------------------------------%
% Setup the linear null space dynamics
%--------------------------------------------------------------------------%
a22s=diag(p,0);
P2=lyap(a22s',eye(pp));

%--------------------------------------------------------------------------%
% Calculate the linear and nonlinear observer gain matrices
%--------------------------------------------------------------------------%
G=Tr'*[A12*inv(C2); (A22*inv(C2)-inv(C2)*a22s)];
F=(P2*C2*B2)';

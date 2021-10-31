function [Gl,Gn,P2]=disobs(A,D,C)

% [Gl,Gn,P2]=disobs(A,D,C) 
%
%         The function returns the linear gain Gl and non-linear gain Gn 
%         along with the Lyapunov matrix P2 which constitutes a discontinuous 
%         observer. For the theory to be valid both C and D should be full 
%         rank  and the Markov parameter C*D must be full rank.
%
%         During the computations a prompt is given for a vector of dimension 
%         n-p-r where n the number of states, p is the number of outputs and 
%         r the number of invariant zeros of (A,D,C). This vector specifies 
%         the desired stable poles of the sliding motion.  
%         In addition the user is prompted for p stable real poles which form
%         part of the state estimation error dynamics.


%         Chris Edwards, Robert Cortez & Sarah Spurgeon
%         Control Systems Research
%         Leicester University
%         University Road
%         Leicester LE1 7RH
%
%         Email: ce@sun.engg.le.ac.uk
%
%         version 1.1 modified from 1.0
%         8/12/97 

%-----------------------------------------------------------------------------%
% This function checks the consistency of the matrices size etc ... and the 
% condition that C*D is full rank.
%-----------------------------------------------------------------------------%
msg=abcchk(A,D,C);
if ~isempty(msg);
  error(msg);
end

[nn,qq]=size(D);
[pp,nn]=size(C);

%-----------------------------------------------------------------------------%
% Bring about a change of coordinates   Acal=[A11 A12] C=[0 I] D=[0 ]
%                                            [A21 A22]           [D2]
%
% where the matrix Acal will now have have a stable top left sub-block
%----------------------------------------------------------------------------%

[Ac,Dc,Cc,Tc]=obsfor(A,D,C);
if isempty(Tc)
   fprintf('No sliding mode observer exists\n')
   P=[];Gl=[];Gn=[];
   return
end


D2=Dc(nn-pp+1:nn,:);

%----------------------------------------------------------------------------%
% Calculate the gains Gl and Gn
%
% The commands below calculate the gain matrices which define the
% observer. The vector p is of dimension pp ( the number of rows of C )
% and specifies the desired stable poles of the error system.

%----------------------------------------------------------------------------%
% Calculate Luenberger gains in the modified coordinates 
%----------------------------------------------------------------------------%
a12=Ac(1:nn-pp,nn-pp+1:nn);
a22=Ac(nn-pp+1:nn,nn-pp+1:nn);

pmsg=['Enter '  num2str(pp) ' output estimation error pole(s) '];
msg=[' '];
while ~isempty(msg);
  p=input(pmsg);
  p=p(:);
  msg=polechk(p,pp,1);
end
a22s=diag(p,0);

%----------------------------------------------------------------------------%
% Compute the Lyapunov matrix to scale the output error
%----------------------------------------------------------------------------%
P2=lyap(a22s',eye(pp));

%----------------------------------------------------------------------------%
% Compute the gain matrices in the original coordinates  
%----------------------------------------------------------------------------%
Gl=inv(Tc)*[a12; (a22-a22s)];
Gn=norm(D2)*inv(Tc)*[zeros(nn-pp,pp); eye(pp)];


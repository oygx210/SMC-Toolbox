function [G,F]=wzobs(A,B,C)

% [G,F]=wzobs(A,B,C) 
%
%         Designs a Walcott and Zak observer for the triple (A,B,C). 
%         The arguments returned are the linear gain G and the matrix F which 
%         in part defines the sliding surface. For the theory to be valid 
%         both C and B should be full rank and the Markov parameter C*B must 
%         be full rank.
%
%         During the computations a prompt is given for a vector of dimension 
%         n-p-r where n the number of states, p is the number of outputs and 
%         r the number of invariant zeros of (A,B,C). This vector specifies 
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
%         version 1.2 modified from 1.1
%         26/5/2000 

%-----------------------------------------------------------------------------%
% This function checks the consistency of the matrices sizes etc ... and the 
% condition that C*B is full rank.
%-----------------------------------------------------------------------------%

msg=abcchk(A,B,C);
if ~isempty(msg);
  error(msg);
end

%-----------------------------------------------------------------------------%
% Bring about a change of coordinates   Acal=[A11 A12] C=[0 I] D=[0 ]
%                                            [A21 A22]           [B2]
%
% where the matrix Acal will now have have a stable top left sub-block
%----------------------------------------------------------------------------%

[nn,mm]=size(B);
[pp,nn]=size(C);
[Ac,Bc,Cc,Tc]=obsfor(A,B,C);
if isempty(Tc)
   fprintf('No sliding mode observer exists\n')
   G=[];F=[];
   return
end

B2=Bc(nn-pp+1:nn,:);
%----------------------------------------------------------------------------%
% Calculate Luenberger gains in the modified coordinates 
%----------------------------------------------------------------------------%
a12=Ac(1:nn-pp,nn-pp+1:nn);
a22=Ac(nn-pp+1:nn,nn-pp+1:nn);

pmsg=['Enter ',num2str(pp),' desired output estimation error pole(s) '];
msg=[' '];
while ~isempty(msg)
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
G=inv(Tc)*[a12; (a22-a22s)];
F=B2'*P2;


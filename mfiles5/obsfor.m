function [Ac,Dc,Cc,Tc]=obsfor(A,D,C,psm)

% [Ac,Dc,Cc,Tc]=obsfor(A,D,C) 
%
%         The function returns the nominal triple (A,D,C) in the canonical
%         form used in the construction of the observer gains (Ac,Dc,Cc). 
%         The transformation matrix to induce this form is given in Tc.

%         Chris Edwards, Robert Cortez & Sarah Spurgeon
%         Control Systems Research
%         Leicester University
%         University Road
%         Leicester LE1 7RH
%
%         Email: ce@sun.engg.le.ac.uk
%
%         version 2.1
%
%         Bug fixed associated with the calculation of L 18/1/2000
%         25/3/98 


%-----------------------------------------------------------------------------%
% This function checks the consistency of the matrices size etc... and the 
% condition that C*D is full rank.
%-----------------------------------------------------------------------------%
msg=abcchk(A,D,C);
if ~isempty(msg)
  error(msg);
end

%-----------------------------------------------------------------------------%
% Bring about a change of coordinates   Acal=[A11 A12] C=[0 I] D=[0 ]
%                                            [A21 A22]           [D2]
%
% where the matrix Acal will now have have a stable top left sub-block
%----------------------------------------------------------------------------%

[nn,qq]=size(D);
[pp,nn]=size(C);
[Af,Df,Cf,Tout,r]=outfor(A,D,C);
if isempty(Tout)
   fprintf('The observer canonical form is not attainable\n')
   Ac=[];Dc=[];C=[];Tc=[];
   return
end

if nn-pp-r>0
    if nargin==3
        pmsg=['Enter '  num2str(nn-pp-r) ' desired stable sliding mode pole(s) '];
        msg=' ';
        while ~isempty(msg)
            p1=input(pmsg);
            p1=p1(:);
            msg=polechk(p1,nn-pp-r);
        end
    else if nargin==4
            msg=polechk(psm,nn-pp-r);
            if ~isempty(msg)
                error(msg);
                return;
            else
                p1=psm(:);
            end
            
        end
    end
    cplxpair(p1);
    A22o=Af(r+1:nn-pp,r+1:nn-pp);
    A21o=Af(nn-pp+1:nn-qq,r+1:nn-pp);
    Ltemp=vplace(A22o',A21o',p1)';
    if r>0
        L=-[zeros(r,pp-qq); Ltemp];
    else
        L=-Ltemp;
    end
else
    L=zeros(nn-pp,pp-qq);
end
Lbar=[L zeros(nn-pp,qq)];

%---------------------------------------------------------------------------%
% Build the final transformation needed to obtain the canonical form
%---------------------------------------------------------------------------%
T=Cf(:,nn-pp+1:nn);
TL=[eye(nn-pp) Lbar; zeros(pp,nn-pp) T];
Ac=TL*Af*inv(TL);                 % In this state-space representation 
Dc=TL*Df;                         % C=[0 I]; D=[0 D2']' and the A matrix  
Cc=Cf*inv(TL);                    % now has a stable top left sub-block 
Tc=TL*Tout;
end


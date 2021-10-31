function [Acal,Dcal,Ccal,Lbar,D2,TL,r]=canon(A,D,C,p1)

%  [L,D2,TL,r]=canon(A,B,C)
%  This function computes where possible the canonical form for
%  the observer in section 6.3.1. If the dimension of the
%  observable subspace is non zero, p1 is a vector of dimension
%  nn-pp-r which specifies the desired stable poles of the subsystem
%  where r is the number of invariant zeros af (A,D,C)
%  Chris Edwards, Robert Cortez & Sarah Spurgeon
%  Control Systems Research
%  Leicester University
%  University Road
%  Leicester LE1 7RH
%
%  Email: ce@sun.engg.le.ac.uk
%
% Version 1.0
% 29/11/97
%

%----------------------------------------------------------------------------%
% Establish the number of inputs and outputs
%----------------------------------------------------------------------------%
[nn,qq]=size(D);
[pp,nn]=size(C);

msg=abcchk(A,D,C);
if ~isempty(msg)
   disp(msg);
   Acal=[];Dcal=[];Ccal=[];L=[];D2=[];TL=[];r=0; 
   return
end

if rank(C*D)~=qq      
   fprintf('The Markov Parameter CB is not full rank\n');
   Acal=[];Dcal=[];Ccal=[];L=[];D2=[];TL=[];r=0; 
   return
end

%----------------------------------------------------------------------------%
% Change coordinates so the output distribution matrix is [0 I]
%----------------------------------------------------------------------------%
nc = null(C);
Tc=[nc'; C];
Ac=Tc*A*inv(Tc);
Dc=Tc*D;
Cc=C*inv(Tc);

%----------------------------------------------------------------------------%
% Partition the input distribution matrix conformably 
%----------------------------------------------------------------------------%
Dc1=Dc(1:nn-pp,:);
Dc2=Dc(nn-pp+1:nn,:);


%----------------------------------------------------------------------------%
% Finds a transformation to bring about a special structure in
% the input and output distribution matrices
%----------------------------------------------------------------------------%
[T,temp]=qr(Dc2);
T=(flipud(T'))';
clear temp;
Tb=[eye(nn-pp) -Dc1*inv(Dc2'*Dc2)*Dc2'; zeros(pp,nn-pp) T'];


Aa=Tb*Ac*inv(Tb);                  % In this new coordinate system
Da=Tb*Dc;                          % we have C=[0 T] and B=[0 B2']'
Ca=Cc*inv(Tb);                     % A has no special structure yet
                    
A11=Aa(1:nn-pp,1:nn-pp);           
A211=Aa(nn-pp+1:nn-qq,1:nn-pp);


%----------------------------------------------------------------------------%
% The aim is to put (A1111,A1121) into the observability canonical form
%----------------------------------------------------------------------------%
[Ab,Bb,Cb,Tobs,k]=obsvf(A11,zeros(nn-pp,1),A211,1000*eps);


%----------------------------------------------------------------------------%
% r is the dimension of the unobservable subspace and
% the number of invariant zeros of the system (A,B,C)
%----------------------------------------------------------------------------%
r=nn-pp-sum(k);

fprintf('Dimension of the unobservable subspace is %.0f \n',r);

Ta=[Tobs zeros(nn-pp,pp);zeros(pp,nn-pp) eye(pp)];

Af=Ta*Aa*inv(Ta);
Df=Ta*Da;
Cf=Ca*inv(Ta);

%----------------------------------------------------------------------------%
% Calculates a gain matrix L so that A11+L*A211 is stable
%----------------------------------------------------------------------------%  
 
if nn-pp-r>0
   A22o=Af(r+1:nn-pp,r+1:nn-pp);
   A21o=Af(nn-pp+1:nn-qq,r+1:nn-pp);
   Ltemp=place(A22o',A21o',p1)';
   L=-inv(Tobs)*[zeros(pp-qq,1:r),Ltemp];
else
    L=zeros(nn-pp,pp-qq);
end
Lbar=[L zeros(nn-pp,qq)];

% Build the final transformation needed to obtain the canonical form
TL=[eye(nn-pp) Lbar; zeros(pp,nn-pp) T];

Acal=TL*Af*inv(TL);     % In this state space representation
Dcal=TL*Df;             % C=[0 I]; D=[0 D2']' and the A matrix 
Ccal=Cf*inv(TL);        % now has a stable top left sub-block

D2=Dcal(nn-pp+1:nn,:);
end
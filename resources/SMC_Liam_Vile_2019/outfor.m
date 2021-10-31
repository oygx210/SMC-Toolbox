function [Af,Bf,Cf,Tcan,r]=outfor(A,B,C)

%  [Af,Bf,Cf,Tcan,r]=outfor(A,B,C)
%          This function computes where possible the canonical form for
%          output feedback design. This requires C*B to be full rank and 
%          the system to be minimum phase. If this is not possible the
%          function returns empty matrices. Otherwise the triple (Af,Bf,Cf)
%          is returned in canonical form. Tcan represents the coordinate
%          transformation necessary to bring this about. The variable r
%          represents the number of invariant zeros


%          Chris Edwards, Robert Cortez & Sarah Spurgeon
%          Control Systems Research
%          Leicester University
%          University Road
%          Leicester LE1 7RH
%
%          Email: ce@sun.engg.le.ac.uk
%
%          Version 1.0
%          29/11/97
%

%----------------------------------------------------------------------------%
% Establish the number of inputs and outputs
%----------------------------------------------------------------------------%
[nn,mm]=size(B);
[pp,nn]=size(C);

msg=abcchk(A,B,C);
if ~isempty(msg)
   disp(msg);
   Af=[];Bf=[];Cf=[];Tcan=[];r=0; 
   return
end

if rank(C*B)~=mm      
   fprintf('The Markov Parameter CB is not full rank\n');
   Af=[];Bf=[];Cf=[];Tcan=[];r=0; 
   return
end

%----------------------------------------------------------------------------%
% Change coordinates so the output distribution matrix is [0 I]
%----------------------------------------------------------------------------%
nc = null(C);
Tc=[nc'; C];
Ac=Tc*A*inv(Tc);
Bc=Tc*B;
Cc=C*inv(Tc);

%----------------------------------------------------------------------------%
% Partition the input distribution matrix conformably 
%----------------------------------------------------------------------------%
Bc1=Bc(1:nn-pp,:);
Bc2=Bc(nn-pp+1:nn,:);


%----------------------------------------------------------------------------%
% Finds a transformation to bring about a special structure in
% the input and output distribution matrices
%----------------------------------------------------------------------------%
[T,temp]=qr(Bc2);
T=(flipud(T'))';
clear temp;
Tb=[eye(nn-pp) -Bc1*inv(Bc2'*Bc2)*Bc2'; zeros(pp,nn-pp) T'];


Aa=Tb*Ac*inv(Tb);                  % In this new coordinate system
Ba=Tb*Bc;                          % we have C=[0 T] and B=[0 B2']'
Ca=Cc*inv(Tb);                     % A has no special structure yet
                    
A1111=Aa(1:nn-pp,1:nn-pp);           
A1121=Aa(nn-pp+1:nn-mm,1:nn-pp);


%----------------------------------------------------------------------------%
% The aim is to put (A1111,A1121) into the observability canonical form
%----------------------------------------------------------------------------%
[Ab,Bb,Cb,Tobs,k]=obsvf(A1111,zeros(nn-pp,1),A1121,1000*eps);


%----------------------------------------------------------------------------%
% r is the dimension of the unobservable subspace and
% the number of invariant zeros of the system (A,B,C)
%----------------------------------------------------------------------------%
r=nn-pp-sum(k);   

Tf=[Tobs zeros(nn-pp,pp);zeros(pp,nn-pp) eye(pp)];

Af=Tf*Aa*inv(Tf);
Bf=Tf*Ba;
Cf=Ca*inv(Tf);

%----------------------------------------------------------------------------%
% If there are any invariant zeros, i.e. r>0, then they must lie in the LHP
%----------------------------------------------------------------------------%  
 
if r>0
   A11o=Af(1:r,1:r);
   if any(real(eig(A11o))>-100*eps)
      fprintf('Unstable invariant zeros are present\n');
      Af=[];Bf=[];Cf=[];Tcan=[];r=0; 
      return
   else
      fprintf('\n')
      fprintf('The system has %.0f stable invariant zero(s)\n',r); 
   end
end
Tcan=Tf*Tb*Tc;

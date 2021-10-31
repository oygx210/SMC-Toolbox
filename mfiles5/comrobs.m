function  [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C,ped,psm,prsd)

% [Hhat,Dhat,S,L,P,Lam,T]=comrobs(A,B,C) 
%
%      Generates a compensator and associated control law. During the 
%      execution the user is prompted to supply 3 vectors:
%
%      an n-p-r dimensional vector (the poles of the estimation error system) 
%      an n-m-r dimensional vector (the poles of effectively the sliding motion)
%      an m dimensional vector (the poles of the range space dynamics)
%
%      where n,m,p and r are the number of states, input, outputs and 
%      invariant zeros respectively.
%
%      The matrices Hhat, Dhat which are returned are the compensators state 
%      and input distribution matrices respectively; S is the sliding surface, 
%      L is the linear feedback component in the control law; P is a Lyapunov
%      matrix which together with Lam constitutes the nonlinear control 
%      component. The orthogonal matrix T is obtained from the canonical form
%      for output feedback design and is used to scale the outputs.

%      ped - desired pole(s) for error dynamics (Lo)
%      psm - desired pole(s) for sliding motion
%      prsd - desired pole(s) for range space dynamics

%      Chris Edwards, Robert Cortez & Sarah Spurgeon
%      Control Systems Research
%      Leicester University
%      University Road
%      Leicester LE1 7RH
%
%      Email: ce@sun.engg.le.ac.uk
%
%      Version 1.1 modified from 1.0
%      9/12/97
%
%
 
%----------------------------------------------------------------------------%
% Check whether the plant is controllable and that the matrices
% conform to various conditions eg size etc...
%----------------------------------------------------------------------------%
msg=abcchk(A,B,C);
if ~isempty(msg)
  error(msg);
end	
				
[nn,mm]=size(B);
[pp,nn]=size(C);
Co=ctrb(A,B);
if rank(Co)~=nn
   wrndlg('The pair (A,B) may not be controllable');
end 	

%----------------------------------------------------------------------------%
% Transform the system triple into the appropriate canonical form
% Bf'=[0 B2']	Cf=[0 T]
%----------------------------------------------------------------------------%

[Af,Bf,Cf,Tcan,r]=outfor(A,B,C);
if isempty(Af)
   error('Unable to attain the canonical form: this approach is not suitable')
end

%----------------------------------------------------------------------------% 
% Prompt the user for vectors p1, p2, and p3 and then check the sizes 
%----------------------------------------------------------------------------%

if nargin==3
    t1=['Enter ' num2str(nn-pp-r) ' pole(s) for the error dynamics '];
    msg=' ';
    while ~isempty(msg)
        p1=input(t1); p1=p1(:);
        msg=polechk(p1,nn-pp-r);
        disp(msg)
    end

    t2=['Enter ' num2str(nn-mm-r) ' pole(s) for sliding motion '];
    msg=' ';
    while ~isempty(msg)
        p2=input(t2); p2=p2(:);
        msg=polechk(p2,nn-mm-r);
        disp(msg)
    end

else if nargin==6
        % error dynamics
        msg=polechk(ped,nn-pp-r);
        if ~isempty(msg)
            error(msg);
            return;
        else
            p1=ped(:);
        end

        % sliding motion
        msg=polechk(psm,nn-mm-r);
        if ~isempty(msg)
            error(msg);
            return;
        else
            p2=psm(:);
        end
end
end

%----------------------------------------------------------------------------%
% If the poles chosen are complex, check if they are in complex conjugate pairs
%----------------------------------------------------------------------------%
p1=cplxpair(p1); 
p2=cplxpair(p2);

if nargin==3
    t3=['Enter ' num2str(mm) ' pole(s) for the range space dynamics '];
    msg=' ';
    while ~isempty(msg)
        p3=input(t3); p3=p3(:);
        msg=polechk(p3,mm,1);
        disp(msg)
    end
else if nargin==6
        % range space dynamics
        msg=polechk(prsd,mm,1);
        if ~isempty(msg)
            error(msg);
            return;
        else
            p3=prsd(:);
        end
end

end

%----------------------------------------------------------------------------%
% Generates the compensator based on a reduced order observer
%----------------------------------------------------------------------------%
T=Cf(:,nn-pp+1:nn);
A11=Af(1:nn-mm,1:nn-mm);
A11o=Af(1:r,1:r);
A12o=Af(1:r,r+1:nn-pp);
A22o=Af(r+1:nn-pp,r+1:nn-pp);
A21o=Af(nn-pp+1:nn-mm,r+1:nn-pp);
A121m=Af(1:r,nn-pp+1:nn-mm);              % Extracts the various sub-blocks
A122m=Af(r+1:nn-pp,nn-pp+1:nn-mm);        % from the system matrix A under
A22m=Af(nn-pp+1:nn-mm,nn-pp+1:nn-mm);     % the assumption the a coordinate
A121=Af(1:r,nn-mm+1:nn);                  % has bee used to obtain the 
A122=Af(r+1:nn-mm,nn-mm+1:nn);            % canonical form of Lemma 5.3
A211=Af(nn-mm+1:nn,1:r);
A212=Af(nn-mm+1:nn,r+1:nn-pp);
A213=Af(nn-mm+1:nn,nn-pp+1:nn-mm);
A22=Af(nn-mm+1:nn,nn-mm+1:nn);

A11tilde=[A22o A122m; A21o A22m];
A1221=A122(1:nn-pp-r,:);
A1222=A122(nn-pp-r+1:nn-mm-r,:);

Lo=vplace(A22o',-A21o',p1);
calK=vplace(A11tilde,A122,p2); 
K1=calK(:,1:nn-r-pp);
K2=calK(:,nn-pp-r+1:nn-mm-r);
Lo=Lo';
H=A22o+Lo*A21o;
D1=A122m+Lo*A22m-H*Lo;
D2=A1221+Lo*A1222;
K=K2-K1*Lo;
Kc=K1;

%----------------------------------------------------------------------------%
% Form the new compensators
%----------------------------------------------------------------------------%
Hhat=[A11o A12o; zeros(nn-pp-r,r) H];
Dhat=[A121m-A12o*Lo A121; D1 D2]*T';

%----------------------------------------------------------------------------%
% Calculate the hyperplane
%----------------------------------------------------------------------------%

S2=eye(mm);
S=S2*[zeros(mm,r) Kc K eye(mm)];

Ahat=[A11o A12o (A121m-A12o*Lo) A121;  zeros(nn-pp-r,r) H D1 D2;
     zeros(pp-mm,r) A21o (A22m-A21o*Lo) A1222; A211 A212 (A213-A212*Lo) A22];

Lam=S*Bf;
Phi=diag(p3);

%---------------------------------------------------------------------------%
% Calculate the linear feedback component
%---------------------------------------------------------------------------%
L=-inv(Lam)*S*Ahat + inv(Lam)*Phi*S; 
P=lyap(Phi',eye(mm));

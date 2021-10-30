function S=rpp(A,B,p)

% S = rpp(A,B,p) 
%
%       Calculates the switching surface using robust eigenstructure assignment
%       methods (under the assumption of full state feedback). The matrix pair
%       (A,B)'s original coordinates are first re-ordered to bring about regular
%       form. The function M = place(A11,A12,p) then computes the state
%       feedback matrix M such that the eigenvalues of  A11-A12*M  are those 
%       specified in vector p. 
%
%       The complex eigenvalues in the vector p must appear in consecutive
%       complex conjugate pairs. These eigenvalues cannot be placed with
%       multiplicity greater than the number of inputs. 
%

%       Chris Edwards, Robert Cortez & Sarah Spurgeon
%       Control Systems Research
%       Leicester University
%       University Road
%       Leicester LE1 7RH
%
%       Email: ce@sun.engg.le.ac.uk 
%
%       Version 1.0
%       29/11/97
%

[nn,mm]=size(B);

%----------------------------------------------------------------------------%
% Check whether the plant is controllable, and that the matrices
% conform to various conditions eg size etc...
%----------------------------------------------------------------------------%
msg=abcchk(A,B);
if ~isempty(msg)
   error(msg);
   return
end

%----------------------------------------------------------------------------%
% Obtain regular form transformation. 
%----------------------------------------------------------------------------%
[A11,A12,B2,Tr]=regfor(A,B);

%----------------------------------------------------------------------------%
% Check if the correct number of poles has been specified,
%----------------------------------------------------------------------------%
p=p(:);
if length(p)~=(nn-mm)
   error(' The number of sliding mode poles chosen is incorrect');
end

% Ensures complex eigenvalues appear in complex conjugate pairs
cplxpair(p);

%----------------------------------------------------------------------------%
% Direct use of the matlab place command gives a robust reduced order
% dynamic via robust eigenstructure assignment
%----------------------------------------------------------------------------%
rnk=rank(A12);
if rnk==mm,
   M=place(A11,A12,p); 
else
  [Q,A12tilde]=qr(A12');
  A12tilde=A12tilde';
  A12tilde=A12tilde(:,1:rnk);
  Mtilde=place(A11,A12tilde,p);
  M=Q*[Mtilde;zeros(mm-rnk,nn)];
end 

% --------------------------------------------------------------------------%
% Recover the switching function matrix in the original coordinates
% --------------------------------------------------------------------------%
S=[M eye(mm)]*Tr;

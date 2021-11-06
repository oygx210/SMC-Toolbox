function [S,V]=dea(A,B,specent,lambda,nocomp)

%    [S,V]=dea(A,B,specent,lambda,nocomp)
%
%          The function DEA ( Direct Eigenstructure Assignment) calculates a
%          hyperplane for the system pair (A,B) and the eigenvectors of the
%          sliding motion. The matrix specent (specified entries) contains the
%          desired eigenvector entries along with any arbitrary values denoted
%          by the variable NaN). If complex entries are desired, the number 
%          nocomp defines the number of pairs of complex conjugate poles and
%          corresponding eigenvector entries. The vector lambda defines the      
%          desired sliding mode poles; the first 2i-1 (i=1, nocomp) entries 
%          contain the real parts of the complex conjugate poles and the 
%          first 2i (i=1,nocomp) contain the imaginary parts. 
%
%          Non-complex eigenvalue calculations are obtained from the function
%          call [S,V]=dea(A,B,specent,p)


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

[nn,mm]=size(B);
[spn,spm]=size(specent);

if nargin<5
   nocomp=0;end
if isempty(nocomp)
   nocomp=0;end

%----------------------------------------------------------------------------%
% Check whether the plant is controllable, and that the matrices
% conform to various conditions eg size etc... 
%----------------------------------------------------------------------------%
msg=abcchk(A,B);                            
if ~isempty(msg)
  error(msg);
end
 
%----------------------------------------------------------------------------%
% Check if the number of sliding mode poles are correct and that no complex
% imaginary numbers are present. 
%----------------------------------------------------------------------------% 

lambda=lambda(:);
if length(lambda)~=(nn-mm)
   error(' The number of sliding mode poles chosen is incorrect');end

for i=1:nn-mm
    if isreal(lambda(i))~=1
       error(' The entries should not have imaginary parts.');
    end
end

%----------------------------------------------------------------------------%
% Create the matrix that indicates whether values are arbitrary or specified
%----------------------------------------------------------------------------%
specpos=~isnan(specent);

% check if the specified matrix is the correct size
if spn~=nn || spm~=(nn-mm)
   error(' The matrix for eigenvector entries is not consistent');end 


% Check if each column of the specified matrix as at least one
% specified entry, since by definition eigenvectors are non-zero                                                        
for i=1:nn-mm
    if norm(double(specpos(:,i)))<100*nn*eps
       error('Specified entries need to be made in each column of the an eigenvector');end
end  

%---------------------------------------------------------------------------%
% If the desired entries are complex then determine the space in which the
% eigenvector corresponding to a specific desired complex conjugate pair of 
% poles must lie.
%---------------------------------------------------------------------------%
for i=1:nocomp

   glambda1=[lambda(2*i-1)*eye(nn)-A lambda(2*i)*eye(nn) B];
   [u,v,w]=svd(glambda1);                                       
   nlambda1=w(1:nn,nn+1:2*nn+mm);                       
   plambda1=w(nn+1:2*nn,nn+1:2*nn+mm);                          
   alpha=[nlambda1; -plambda1]';                
   beta=[plambda1; nlambda1]';
   [u1,v1,w1]=svd(alpha);
   [u2,v2,w2]=svd(beta);
   gamma=[w1(:,nn+mm+1:2*nn)'; w2(:,nn+mm+1:2*nn)'];
   [u3,v3,w3]=svd(gamma);
   kgamma=w3(:,2*nn-2*mm+1:2*nn);

   % Find subvector of specified entries 
   
   despos=find(specpos(:,2*i-1));
   numspecr=length(despos);
   n1=[]; desent=[];
   for j=1:numspecr
      n1(j,:)=kgamma(despos(j),:);
      desent(j)=specent(despos(j),2*i-1);
   end
   despos=find(specpos(:,2*i));
   numspecc=length(despos);
   for j=1:numspecc
      n1(j+numspecr,:)=kgamma(despos(j)+nn,:);
      desent(j+numspecr)=specent(despos(j),2*i);
   end

   % Find solution for delta

   delta=n1\desent';
   if norm(delta) <=100*mm*eps
      %warndlg(' The specified entries produce a zero eigenvector');
   end

   vector=kgamma*delta;
   
   V(1:nn,2*i-1)=vector(1:nn);
   V(1:nn,2*i)=vector(nn+1:2*nn);
   
end

%---------------------------------------------------------------------------%
% Determine the space in which the eigenvector corresponding to a specific real
% desired pole must lie.
%---------------------------------------------------------------------------%
for i=2*nocomp+1:(nn-mm)
 
   glambda=[lambda(i)*eye(nn)-A B];
   [u,v,w]=svd(glambda);
   nlambda=w(1:nn,nn+1:nn+mm);

   % Find subvector of specified entries
   despos=find(specpos(:,i));                 
   numspec=length(despos);                    
   n2=[]; desent2=[];                           
   for j=1:numspec
    n2(j,:)=nlambda(despos(j),:);
    desent2(j)=specent(despos(j),i);      
   end
      
  % Find solution for delta
   delta=n2\desent2';
   if norm(delta) < 100*mm*eps
      %warndlg(' The specified entries have produced a zero eigenvector');
   end
   V(:,i)=nlambda*delta;

end

% Use a singular value decomposition to determine the switching 
% function matrix from the selected eigenvectors

[u,v,w]=svd(V);
S=u(:,nn-mm+1:nn)';

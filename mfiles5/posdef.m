function   t=posdef(X,msg)

%    posdef(X,msg)
%
%         Tests whether the matrix X is symmetric positive definite. 
%         If it is the value 1 is returned, otherwise 0 is returned.
%         msg is a string variable which specifies the variable name  
%         of the matrix tested in a dialogue box.


%         Chris Edwards, Robert Cortez & Sarah Spurgeon
%         Control Systems Research
%         Leicester University
%         University Road
%         Leicester LE1 7RH
%
%         Email: ce@sun.engg.le.ac.uk 
%
%         Version 1.0
%         29/11/97
%

if nargin>1
   msg=[msg,' '];
else
   msg=[];
end

% Test for symmetry and positive definiteness
nnx = norm(X,1);
if any(eig(X) <= eps*nnx) | (norm(X'-X,1)/nnx > eps)
  warndlg(['The matrix ',msg,'is not symmetric and positive definite']);
  t=1;
else
  t=0;
end

function   t=stable(X,msg)

%    stable(X,msg)
%
%         Tests whether the matrix X is stable. If it is stable
%         the value 1 is returned, otherwise 0 is returned.
%         msg is a string variable which specifies the variable name of 
%         the unstable matrix in a dialogue box.

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

% Test for stability
if any(real(eig(X))>-eps)
  warndlg(['The matrix ',msg,'is unstable']);
  t=1;
else
  t=0;
end

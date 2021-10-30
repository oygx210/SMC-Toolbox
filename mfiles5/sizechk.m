function   t=sizechk(X,n,m,msg)

%    sizechk(X,n,m,msg)
%
%         Tests whether the matrix X has n rows and m columns if it has 
%         then the value 1 is returned, otherwise 0 is returned.
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

if nargin>3
   msg=[msg,' '];
else
   msg=[];
end

% Test the size of the matrix
[nx,mx]=size(X);
if (nx~=n | mx~=m)
  warndlg(['The matrix ',msg,'is of inappropriate size']);
  t=1;
else
  t=0;
end

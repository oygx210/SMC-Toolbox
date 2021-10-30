function   msg=polechk(p,n,flagc)

% polechk(p,n,flagc)
%
%       Checks whether a given vector is suitable as a set of closed loop poles.
%       Specifically it checks that the correct number of poles have been given
%       i.e. length(p)=n, that they are all stable and, if complex poles have
%       been specified, that they occur in complex conjugate pairs.
%
%       The flagc is optional and if it is included it ensures that only real
%       poles have been entered.
%
%       If the vector p passes all these restriction then an empty string is
%       returned, otherwise a string indicating the reason for failure is  
%       returned.
%

%       Chris Edwards, Robert Cortez & Sarah Spurgeon
%       Control Systems Research
%       Leicester University
%       University Road
%       Leicester LE1 7RH
%
%       Email: ce@sun.engg.le.ac.uk
%
%        Version 1.0
%        29/11/97
%

if nargin<3
  flagc=0;
else 
  flagc=1;
end

p=p(:);
msg=[];

if flagc==1 & isreal(p)~=1 
   msg='Only real poles are allowed';
else
   if length(p)~=n   
      msg=['Incorrect number of poles specified: ',num2str(n), ' were expected'];
   else
      if any(real(p)>-100*eps)
         msg='Only stable poles can be entered';
      else 
         if isempty(cplxchk(p))
            msg='Complex poles exist which cannot be paired';
          end
      end
   end
end

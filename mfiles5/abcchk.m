function msg=abcchk(A,B,C)

% msg=abcchk(A,B,C) 
%
%          When only two arguments are supplied the mfile checks whether 
%          the matrices A,B are of consistent dimension, B has full rank
%          and the pair (A,B) is controllable.     
%          When three arguments are supplied in addition to the checks above
%          the rank and dimension of C is checked.    

%          Chris Edwards, Robert Cortez & Sarah Spurgeon
%          Control Systems Research
%          Leicester University
%          University Road
%          Leicester LE1 7RH
%
%          Email: ce@sun.engg.le.ac.uk
%
%          version 1.2
%         
%          Modified 14/2/01 with regard to controlability checking
%          The controlability matrix is now no longer formed. The
%          controlability is now done through the `staircase form' 

[nx,ny]=size(A);
[nn,mm]=size(B);

msg=[];
if isempty(A) || isempty(B)
   msg='A and B cannot be empty';
   return
end

% Check size of A and B with number of states nx, and inputs mm

if  nx~=ny
  msg='A should be square';
  return
end
if (nx~=nn) || (nx<mm)
  msg='The dimensions of A and B are not consistent';
  return
end



% check whether B is full rank
rb=rank(B);

if rb~=mm
   msg='The input distribution matrix is rank deficient';
   return
end

  
if nargin==3
   [pp,np]=size(C);       
   if np~=nn || pp>nn
      msg='The output distribution matrix is not consistent';
      return;
   end
   if pp<mm 
      msg='The output distribution matrix has too few arguments';
      return
   end 
   rc=rank(C);
   if rc~=pp
      msg='The ouput distribution matrix is rank deficient\n';
      return
   end
else
  [Ac,Bc,Cc,Tc,k]=ctrbf(A,B,zeros(1,nn),1000*eps);
  u=nn-sum(k);
  if u>0
     msg='The matrix pair (A,B) may not be controllable';
     if u>0
        Ac11=Ac(1:u,1:u);
        if any(real(eig(Ac11))>-100*eps)
           msg=[msg ': worse its unstabilisable'];      
        end
     end
     warndlg(msg);
     msg=[];
  end 
  return
end   

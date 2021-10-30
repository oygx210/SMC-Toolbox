function specent=eigbuil(n,m,nocomp)
%
%  specent=eigbuil(n,m,nocomp)
%
%       Builds up a matrix containing information about the eigenvector 
%       elements in a form suitable for DEA. The scalar n is the size of the
%       state space vector, m the number of inputs and nocomp the complex pairs
%       of eigenvalues. If only two right arguments are presented, nocomp is
%       taken to be zero. The size of the returned vector is n by n-m. 

%       Chris Edwards, Robert Cortez & Sarah Spurgeon
%       Control Systems Research
%       Leicester University
%       University Road
%       Leicester LE1 7RH
%
%       Email: ce@sun.engg.le.ac.uk
%
if nargin==2
  nocomp=0
end

if n<=m
  sprintf('Too many inputs ((%.f)) for the number of states ((%.f))',m,n)
end

for jj=1:n-m
  for ii=1:n
    specent(ii,jj)=NaN;
  end
end

for jj=1:n-m
  if nocomp==0
    msg=['Eigenvector ',num2str(jj)];
  else
    jjc=fix((jj+1)/2);
    if jjc<=nocomp
      if 2*jjc==jj
        msg=['Eigenvector ',num2str(jjc),' (Imaginary part)'];
      else
        msg=['Eigenvector ',num2str(jjc),' (Real Part)'];
      end
    else
      msg=['Eigenvector ',num2str(jj-nocomp)];
    end
  end
  disp(msg)
  disp('Enter position of specified entry [0 to terminate]');
  disp('Followed by the value of specified entry');
  ii=1;
  while ii>0
    ii=input('Enter position ');
    ii=fix(ii);
    if ii>0 & ii<=n
      specent(ii,jj)=input('Enter associated value ');
    end
  end
  disp('---------------------------------------------------')
  if sum(isnan(specent(:,jj)))>n-2
    warndlg('Less than 2 elements have been specified in this eigenvector');
  end
  fprintf('\n')
end

return

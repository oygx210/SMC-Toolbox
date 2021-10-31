function K=vplace(A,B,p)

%    k=vplace(A,B,p)
%
%        This command is a modification of the pole placement routine PLACE
%        An initial transformation is made to remove any redundant columns of
%        the input distribution matrix.


%        Chris Edwards, Robert Cortez & Sarah Spurgeon
%        Control Systems Research
%        Leicester University
%        University Road
%        Leicester LE1 7RH
%
%        Email: ce@sun.engg.le.ac.uk
%
%        Version 1.0
%        29/11/97
%

[nn,mm]=size(B);

rnk=rank(B);
if rnk==mm,
   K=place(A,B,p); 
else
  [Q,Btilde]=qr(B');
  Btilde=Btilde';
  Btilde=Btilde(:,1:rnk);
  Ktilde=place(A,Btilde,p);
  K=Q*[Ktilde;zeros(mm-rnk,nn)];
end 

return

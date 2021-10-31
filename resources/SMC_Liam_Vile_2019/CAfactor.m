function [Bv, B2, Tv]=CAfactor(B)


 [nn, mm]=size(B);
 
 [Tr temp]=qr(B); 
 
 Tv=fliplr(Tr');
 Breg=Tr*B;
     
Bsize=[];
    for n=1:nn
        
        Bsize=[Bsize; norm(B(n,:))];
        
    end
 
ll=nn-size(find((Bsize/max(Bsize))<=0.001),1); %Calculate Negligable rows

B2 = Breg(nn-ll+1:nn,1:mm);                   %Calculate B2

Bv=[zeros(nn-ll, ll); eye(ll)];               %New Virtual Control Matrix 

end
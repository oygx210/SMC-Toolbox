function [A, Bv, Tr, CA]=SMCCanForm(A, B)

[nn mm]=size(B);


%--------------------------------------------------------------------------%
% Perform QR decomposition on the input distribution matrix
%--------------------------------------------------------------------------%
[Tr temp]=qr(B); 
Tr=flip(Tr');


%--------------------------------------------------------------------------%
% Obtain (Areg,Breg); regular form description
%--------------------------------------------------------------------------%
A=Tr*A*Tr';
B=Tr*B;

BB=round(B/(max(max(abs(B)))), 2);

ll=rank(BB);
%--------------------------------------------------------------------------%
% Obtain matrix sub-blocks for sliding mode controller design
%--------------------------------------------------------------------------%
A11 = A(1:nn-ll,1:nn-ll);
A12 = A(1:nn-ll,nn-ll+1:nn);
B2 = B(nn-ll+1:nn,1:end);

Bv=[zeros(nn-ll ,ll)
    eye(ll)];

CA=pinv(B2);

end

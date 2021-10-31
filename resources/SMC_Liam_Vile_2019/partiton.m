function [A11, A12]=partiton(A, Bv)


[nn ll]=size(Bv);

A11=A(1:nn-ll, 1:nn-ll);
A12=A(1:nn-ll, nn-ll+1:nn); 

end
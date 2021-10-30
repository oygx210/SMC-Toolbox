function [A11,A12,B2,Tr]=regfor(A,B)

% [A11,A12,B2,Tr]=regfor(A,B) 
%        Returns the sub-matrices and the orthogonal transformation matrix that 
%        brings about the canonical form. The pair A,B are partitioned in the 
%        form
%          
%        Tr*A*Tr'=[A11 A12]	Tr*B=[0 ] 
%                 [A21 A22]	     [B2]
%
%

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

%--------------------------------------------------------------------------%
% Establish the size of the input distribution matrix
%--------------------------------------------------------------------------%
[nn,mm]=size(B);

%--------------------------------------------------------------------------%
% Perform QR decomposition on the input distribution matrix
%--------------------------------------------------------------------------%
[Tr temp]=qr(B); 
Tr=Tr';
Tr=[Tr(mm+1:nn,:);Tr(1:mm,:)];

%--------------------------------------------------------------------------%
% Obtain (Areg,Breg); regular form description
%--------------------------------------------------------------------------%
Areg=Tr*A*Tr';
Breg=Tr*B;

%--------------------------------------------------------------------------%
% Obtain matrix sub-blocks for sliding mode controller design
%--------------------------------------------------------------------------%
A11 = Areg(1:nn-mm,1:nn-mm);
A12 = Areg(1:nn-mm,nn-mm+1:nn);
B2 = Breg(nn-mm+1:nn,1:mm);

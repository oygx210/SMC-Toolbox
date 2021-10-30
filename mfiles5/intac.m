function [AT,BT]=intac(A,B,C)

% [AT,BT]=intac(A,B,C)
%
%          For an integral action based approach in sliding mode theory the 
%          the augmented system 
%                       AT=[0 -C]  BT=[0] 
%                          [0  A]     [B]
%          is first formed around which the hyperplane is designed. 
%          If this approach is not applicable because of the presence 
%          of zeros at the origin empty matrices are returned and an
%          appropriate error message is displayed


%          Chris Edwards, Robert Cortez & Sarah Spurgeon
%          Control Systems Research
%          Leicester University
%          University Road
%          Leicester LE1 7RH
%
%          Email: ce@sun.engg.le.ac.uk
%
%
%          Version 1.0
%          29/11/97
%

[pp,nn]=size(C);
[nn,mm]=size(B);

%-----------------------------------------------------------------------------%
% Check whether the plant is controllable and that the matrices
% conform to various conditions eg size etc...
%-----------------------------------------------------------------------------%

AT=[];
BT=[];

msg=abcchk(A,B);                             
if ~isempty(msg)
  error(msg)
end

AT=[zeros(pp,pp) -C;zeros(nn,pp) A];
BT=[zeros(pp,mm);B];

%----------------------------------------------------------------------------
% Check that (AT,BT) is controllable. If it is not then the triple (A,B,C)
% has invariant zeros at the origin
%----------------------------------------------------------------------------

if rank(ctrb(AT,BT))~=nn+mm
  error('The triple (A,B,C) may have zeros at the origin');
end 


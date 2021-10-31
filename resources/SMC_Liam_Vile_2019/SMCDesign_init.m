function [L, P, Lam, S, E, CA]=SMCDesign_init(A, B, Q)

[At, Bv, Tr, CA]=SMCCanForm(A, B);

[S,E]=lqcfCA(At,Bv,Q);

Phi=-eye(size(S,1));

[L,P,Lam]=contl(At,Bv,S,Phi);

L=L*Tr;
S=S*Tr;

end

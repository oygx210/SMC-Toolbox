load compped.mat;
[Af,Bf,Cf,Tcan,r]=outfor(A,B,C);

[n,m]=size(B);
[p,n]=size(C);

x0=zeros(1,n);
x0=[0 0 0.1 0];

rho=0.1;

S=[1 1 1 1];
Phi=-5;
[L,P,Lam]=contl(A,B,S,Phi);

%%
mdl_name='streg_mdl';
if ~bdIsLoaded(mdl_name)
    open('streg_mdl.slx');
end
sim(mdl_name);

%%
subplot(3,1,1)
plot(t.Data,u.Data);
grid on;
legend({'u'});
subplot(3,1,2)
plot(t.Data,x.Data);
grid on;
legend({'x'});
subplot(3,1,3)
plot(t.Data,s.Data);
legend({'s'});
grid on;
legend({'yobs'});

%plot(y.Time,y.Data(:,3),yobs.Time,yobs.Data(:,3));
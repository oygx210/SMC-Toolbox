% load the data from mat-file
clear all;
clc;

data=load('l1011r.mat');
A=data.A;
B=data.B;
C=data.C;

[n,m]=size(B);
[p,n]=size(C);

%% Chap. 5.7.1 Model-Reference System Using only Outputs

% design the swithcing surface
lambda = [-0.05 -2 1.5 -1.5 1.5];
nocomp = 1;
specent=[  1 NaN NaN NaN NaN;
         NaN   1   0 NaN NaN;
         NaN   0   1 NaN NaN;
         NaN NaN NaN   1 NaN;
         0 NaN NaN NaN   1];

%[S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
[S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
S=inv(S*B)*S;

%S=rpp(A,B,lambda(1:n-m));

% check if the poles are complex
p = lambda(n-m+1:end);
if any(p>0)
    cp_idx=find(p>0);
    D=eye(2*length(cp_idx));
    for i=1:length(cp_idx)
        D(i,i)=p(cp_idx(i)-1)+1i*p(cp_idx(i));
        D(i+1,i+1)=p(cp_idx(i)-1)-1i*p(cp_idx(i));
    end
    [V,D]=cdf2rdf(eye(m),D);
    Phi = D;
else
    Phi=diag(p);
end

Phi=diag([-5 -5]);

[Lx,P,Lam]=contl(A,B,S,Phi);

% Design diagonal weighting matrix for the state vector
Q=diag([5 1 1 5 5]); % compare with Eq. 4.139

% Designs a hyperplane S by minimising the linear quadratic cost function
% [S,E]=lqcf(A,B,Q);
% L=-inv(S*B)*(S*Am-Phi*S);



dLx = 1;

iter=0;

while dLx > 1e-1
    
    nanfound=1;
    while nanfound
        specent=randi([0,2],[n,n-m]);
        specent(specent==2)=NaN;
        colallnan=0;
        for i=1:(n-m)
            if all(isnan(specent(:,i)))
                colallnan=1;
            elseif all(specent(:,i) == specent(1,i))
                    colallnan=1;
            elseif ~any(isnan(specent(:,i)))
                    colallnan=1;
            end
        end

        if ~colallnan
            nanfound=0;
        end
    end
    
    try
        [S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
        S=inv(S*B)*S;
        [Lx,~,~]=contl(A,B,S,Phi);
        if ~isnan(Lx(:))
            dLx=max(abs(data.Lx(:)-Lx(:)));
        end
    catch ME
        rethrow(ME);
    end

    iter=iter+1;
    disp(iter);
    %disp(specent);
end
disp(specent);
error('specent found');

% try random matrix
numbers=[0,1,NaN];
for i=1:size(specent,1)
    for j=1:size(specent,2)
        for k=1:length(numbers)
            specent(i,j) = numbers(k);
            [Lx,P,Lam]=contl(A,B,S,Phi); % compare Lx (page 179) with F (Eq. 4.137)
            if abs(data.Lx-Lx) < 1e-2
                disp(specent);
                error('specent found');
            end
        end
    end
end

% H=[1 0 1 0 0;
%    1 0 0 0 0];
% G=-inv(H*inv(A+B*F)*B);

% Model-Reference matrices
Am=A+B*Lx;
eig(Am)
%Bm=B*G;

% Design diagonal weighting matrix for the state vector
Q=diag([5 1 1 5 5]);

% Designs a hyperplane S by minimising the linear quadratic cost function
[S,E]=lqcf(A,B,Q);
L=-inv(S*B)*(S*Am-Phi*S);

% Design parameter for switching function 
rho=1;
delta=0.001;

% Initial state vector
x0=zeros(n);

%% compare calculated valued with original
fields = fieldnames(data);
for i=1:length(fields)
    if abs(data.(fields{i})-eval(fields{i})) > 1e-2
        error([fields{i},' is not equal'])
    end
end

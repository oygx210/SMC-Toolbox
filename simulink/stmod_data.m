%% clear the workspace
clear all;

%% Load the system data
%mat_name = 'aircraft';
%mat_name = 'compped';
%mat_name = 'dcmotor';
mat_name = 'furnobs';
%mat_name = 'hecka';
%mat_name = 'helcopt';
%mat_name = 'huiex';
%mat_name = 'invpen';
%mat_name = 'l1011r';
%mat_name = 'oxymodel';
%mat_name = 'vertint';

load([mat_name,'.mat']);

%% Get system dimensions
[n,m]=size(B);
[p,n]=size(C);

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
% x0(2)=1;
% x0(3)=1;

%% Design of the compensator
% Chap. 4.5.1 Model-Reference Design

switch mat_name
    case 'furnobs'
        % desired sliding motion poles
        psm=[-0.025,-0.03+1i*0.025,-0.03-1i*0.025];

        [~,~,~,~,r]=outfor(A,B,C);
        if n-p-r>0
            [G,F]=wzobs(A,B,C,psm,-0.2);
        else
            [G,F]=wzobs(A,B,C,[],-0.2);
        end

        Gamma=diag([-0.025]);

        Phi=-0.1;
        [L,P,Lam]=contl(A,B,S,Phi);
        SimStopTime=300;

    case 'vertint'
        % design the swithcing surface
        lambda = [-5.6 4.2 -1 -20 -20];
        nocomp = 1;
        specent=[  0   0   1 NaN NaN;
            1 NaN   0 NaN NaN;
            NaN   1   0 NaN NaN;
            NaN NaN NaN   1 NaN;
            NaN NaN NaN NaN   1];

        [S,V]=dea(A,B,specent(:,1:n-m),lambda(1:n-m),nocomp);
        S=inv(S*B)*S;

        Phi=diag(lambda(n-m+1:end));
        [F,P,Lam]=contl(A,B,S,Phi);

        H=[1 0 1 0 0;
            1 0 0 0 0]; % compare with Eq. 4.131
        G=-inv(H*inv(A+B*F)*B); % compare with Eq. 4.133

        % Model-Reference matrices
        Am=A+B*F;
        Bm=B*G;

        % Design diagonal weighting matrix for the state vector
        Q=diag([10 5 5 20 20]); % compare with Eq. 4.139

        % Designs a hyperplane S by minimising the linear quadratic cost function
        [S,E]=lqcf(A,B,Q); % compare with Eq. 4.140
        L=-inv(S*B)*(S*Am-Phi*S);  % compare with Eq. 4.141

        % Design parameter for switching function
        rho=1;
        delta=0.001;

        % Initial state vector
        x0=zeros(1,n);
        xm0=zeros(1,n);

        % check roots of closed loop
        % The condition number of the eigenvectors is an indicator of the
        % robustness with respect to unmatched uncertainty (smaller => better)
        M=inv(S*B)*S;
        [V,D] = eig(A-B*M);
        cond(V)

        SimStopTime=10;
end

%% Simulate the model
mdl_name='stmod_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% Plot the results
subplot(4,1,1)
plot(t.Data,u.Data);
grid on;
legend({'u'});
subplot(4,1,2)
plot(t.Data,x.Data);
grid on;
legend({'x'});
subplot(4,1,3)
plot(t.Data,xm.Data);
grid on;
legend({'xm'});
subplot(4,1,4)
plot(t.Data,s.Data);
legend({'s'});
grid on;
legend({'s'});
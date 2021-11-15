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


[n,m]=size(B);
[p,n]=size(C);

%% Set initial conditions
x0=zeros(1,n);
x0(1)=1; % set first state initial value
% x0(2)=1;
% x0(3)=1;

switch mat_name
    case 'vertint'

%         C=C(2,:); % partial output matrix to allow a reduced observer design
%         C=[C; [0,1,0,0,0]];
%         [p,n]=size(C);
% 
%         sys_info(ss(A,B,C,zeros(p,m)),mat_name);

        [Aa,Ba]=intac(A,B,C);

        %
        nocomp=1;
        lambda=[-5.6 4.2 -1 -0.4 -0.7];
        %specent=eigbuil(7,2,nocomp)
        % 4 Eigenvectors (only the first one is complex)
        %        1re 1im   2   3   4
        specent=[NaN NaN NaN NaN NaN;
                 NaN NaN NaN NaN NaN;
                   0   0   1   0   1;
                   1 NaN   0   1   0;
                 NaN   1   0 NaN   0;
                 NaN NaN NaN NaN NaN;
                 NaN NaN NaN NaN NaN];

        [S,V]=dea(Aa,Ba,specent,lambda,nocomp);

        Phi=-20*eye(size(Ba,2));

        S=inv(S*Ba)*S;

        % Calculate the switching surface using robust eigenstructure assignment
        Sob=rpp(A,B,[-5.6+1i*4.2,-5.6-1i*4.2,-1]);
        [F,~,~]=contl(A,B,Sob,Phi);

        G=-inv(C*inv(A+B*F)*B);
        %C=[C; [0,1,0,0,0];[0,0,1,0,0];[0,0,0,1,0]];
        %[G,F]=wzobs(A,B,C,[-5.6+1i*4.2,-5.6-1i*4.2,-1,-0.4,-0.7],-10);

        x0=zeros(1,n); % [xr, x11, x12]
        x0(2)=1;
        z0=zeros(1,n);

        Gamma=diag([-0.9 -0.7]);

        rho=1;
        delta=0.001;

        SimStopTime=10;

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

        nn=n;pp=p;mm=m;
        %----------------------------------------------------------------------------%
        % Augment the statespace with integral action states
        %----------------------------------------------------------------------------%
        AT=[zeros(pp,pp) -C; zeros(nn,pp) A];
        A11a=AT(1:nn,1:nn);
        A12a=AT(1:nn,nn+1:nn+pp);

        Co=ctrb(A11a,A12a);
        unco = length(A) - rank(Co);

        if unco==0
            Tr=eye(nn);
        else
            [~,~,~,~,~,Tr]=regfor(A,B);
            %BT=[zeros(pp,mm); B];
        end

        %----------------------------------------------------------------------------%
        % Change coordinates so that regular form is achieved
        %----------------------------------------------------------------------------%
        Trtilde=[eye(pp) zeros(pp,nn);zeros(nn,pp) Tr];
        ATnew=Trtilde*AT*Trtilde';
        A11a=ATnew(1:nn,1:nn);
        A12a=ATnew(1:nn,nn+1:nn+pp);

        M=place(A11a,A12a,psm);
        S=[M eye(mm)]*Trtilde';

        x0=zeros(1,n); % [xr, x11, x12]
        x0(1)=0;
        z0=zeros(1,n);

        rho=1;
        delta=1.0;%0.001;

        SimStopTime=500;
end

% Design of observer based controller with integral action
[L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);

% %% Design of the compensator
% rho=0.5; % user defined parameter
% 
% % hyperlane design with desired poles for system dynamics
% wn=2.0; %  [rad/s] cut-off angular frequency
% zeta=0.95; % [-] damping ratio less 1 --> complex poles
% sp=roots([1 2*zeta*wn wn^2]); % second order for first 2 poles
% 
% % additional poles (real) by system order greater than 3
% if n-m>2
%     % to avoid multiple poles iterate through and increase each value
%     for i=1:(n-m-2)
%         sp=[sp; -i];  % [-1+i1; -1-i1; -1; -2; -3 ...]
%     end
% end
% 
% % Calculate the switching surface using robust eigenstructure assignment
% S=rpp(A,B,sp);
% 
% % Calculate L and Lam required for the linear component of the control law
% Phi=-5*eye(size(S,1));
% 
% % Design of observer based controller with integral action
% [L,Lr,Lrdot,Sr,Lam,P]=contlia(A,B,C,S,Phi);
% 
% % check roots of closed loop
% % The condition number of the eigenvectors is an indicator of the
% % robustness with respect to unmatched uncertainty (smaller => better)
% M=inv(S*B)*S;
% [V,D] = eig(A-B*M);
% cond(V)

%% Simulate the model
mdl_name='stint_mdl';
if ~bdIsLoaded(mdl_name)
    open([mdl_name,'.slx']);
end
set_param(mdl_name,'StopTime',sprintf('%d',SimStopTime));
sim(mdl_name);

%% Plot the results
subplot(2,3,1)
plot(t.Data,u.Data);
grid on;
legend(get_legend('u'));

subplot(2,3,2)
plot(t.Data,s.Data);
grid on;
legend(get_legend('s'));

subplot(2,3,3)
plot(t.Data,r.Data);
grid on;
legend(get_legend('r'));

subplot(2,3,4)
plot(t.Data,y.Data);
grid on;
legend(get_legend('y'));

subplot(2,3,5)
plot(t.Data,ym.Data);
grid on;
legend(get_legend('ym'));

subplot(2,3,6)
plot(t.Data,ey.Data);
grid on;
legend(get_legend('ey'));

sgtitle([mdl_name,' - ',mat_name],'Interpreter','None');

%% MPC on HE
clc; clear; yalmip('clear'); close all;

% When generates flat figures
set(0, 'DefaultFigureRenderer', 'painters');

%% Load polytope and observer matrices
% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
RUIO = struct;
UIOO = struct;
load linObs

% %% LTI model and observers
% % This section is commented to reduce simulation time (using pre-calculated observer matrices)
% Ts = 0.05;               % Sample time [min]
% Theta_1s = 498;     % Temperatura de salida de fluido 1 (K) (set-point) (498)
% Theta_2s = 690;     % Temperatura de salida de fluido 2 (K)
% 
% run HE;
%
% %% Reduced-order unknown input observer
% N = 2;
% run HE_RUIO;
% 
% %% Unknown input output observer
% run HE_UIOO;
% 
% % Save observers' data
% save linObs.mat

%% Simulation parameters
Time = 50;                              % Simulation end time 
Nsim = Time/Ts;                     % Simulation steps
t = 0:Ts:Time-Ts;                    % Simulation time

Fail_Q1 = 5; Fail_Q2 = 0.43;    % Actuator fault magnitude [5%, 5%]
Fail_S1 = 1.5; Fail_S2 = -1.5;	% Sensor fault magnitude [0.5% 0.5%]
x0 = [495; 680; 570];             % Start-point
xsp = [Theta_1s; Theta_2s; Theta_p];% Set-point
Np = 5;                                   % Prediction horizon

% Constraints
xmin = [495; 650; 530];
xmax = [500; 750; 590];
umin = [90; 6];
umax = [110; 10];

% Weight matrix
Qx = eye(nx);
% Qx = diag([1 1 1]);
% Rx = eye(nu);
Rx = diag([1 0.1]);
gamma = 1e2*diag([1e6 1 1]); 

%% Noise
sig = 3e-3*([1 1 1])';         % Ouput noise sigma

rng default;                        % Random seed start
v = sig*randn(1, Nsim);    % Measurement noise v~N(0, sig)

%% Max reachable set
Xpoly = Polyhedron('lb', xmin-X_lin, 'ub', xmax-X_lin);     % State polyhedron
Upoly = Polyhedron('lb', umin-U_lin, 'ub', umax-U_lin);    % Input polyhedron
Acl = [Ad, Bd; zeros(nu, nx) eye(nu)];
Maxiter = 50;
Xs = Xpoly;
for k = 1:Maxiter
    Xo = Xs;
    Z = Xo*Upoly;   % Extended Set
    
    % Next set
    S = Z.invAffineMap(Acl);
    S = S.intersect(Z);
    Xs = S.projection(1:nx).minHRep();
    Xs = Xs.intersect(Xo).minHRep();
    
    % Check Invariant
	if Xs == Xo
        break;
	end
end

%% Error detection threshold
Tau = 2;               % Convergence period
mag_1 = 1e-1;     % Value Q1
mag_2 = 5e-2;     % Value Q2
mag_3 = 5e-4;     % Value O1
mag_4 = 2e-3;     % Value O2

threshold = zeros(4, Nsim);

for k = 1:Nsim
    threshold(1, k) = mag_1 + 1000*exp(-(k-1)/Tau);  % Q1
    threshold(2, k) = mag_2 + 900*exp(-(k-1)/Tau);    % Q2
    threshold(3, k) = mag_3 + 100*exp(-(k-1)/Tau);    % O1
    threshold(4, k) = mag_4 + 1000*exp(-(k-1)/Tau);  % O2
end

%% MPC controller
run MPC;

%% Simulation Setup
% Fault Tolerant Control System (1 = enable, 2 = disable).Matrix 
FTCS = struct;
for j = 1:2
    FTCS(j).U = zeros(nu, Nsim);                   % Control Input
    FTCS(j).Ufail = zeros(nu, Nsim);              % Faulty control Input
    FTCS(j).Uff = zeros(nu, Nsim+1);            % Feedforward control Input
    FTCS(j).Ufails = zeros(nu, Nsim);             % Fails of control inputs
    FTCS(j).X = zeros(nx, Nsim+1);               % States
    FTCS(j).Y = zeros(ny, Nsim);                    % Measure outputs
    FTCS(j).Yfail = zeros(ny, Nsim);                % Faulty measure outputs
    FTCS(j).Obj = zeros(1, Nsim);                   % Objective cost
end

% RUIOs
for j = 1:N
    RUIO(j).Phi = zeros(N, Nsim+1);     % Observer states
    RUIO(j).X = zeros(nx, Nsim);           % Estimated states
    RUIO(j).error = zeros(1, Nsim);       % Error
    RUIO(j).Fact = zeros(1, Nsim);        % Estimated control Input
    RUIO(j).FQ = zeros(1, Nsim);           % Fault detect Q
    RUIO(j).delay = 0;                            % Detection delay
end

% UIOOs
for j = 1:N
    UIOO(j).Z = zeros(nx, Nsim+1);      % Observer states
    UIOO(j).Ymon = zeros(N, Nsim);     % Monitorated outputs
    UIOO(j).X = zeros(nx, Nsim);           % Estimated states
    UIOO(j).Y = zeros(nx, Nsim);           % Estimated outputs
    UIOO(j).res = zeros(nx, Nsim);        % Residue
    UIOO(j).error = zeros(1, Nsim);       % Error
    UIOO(j).Fsen = zeros(1, Nsim);       % Estimated sensor fault
    UIOO(j).FO = zeros(1, Nsim);          % Fault detect S
end

% Initial states and inputs
for j = 1:2
    FTCS(j).X(:, 1) = x0;
    FTCS(j).Y(:, 1) = C*x0;
    FTCS(j).U(:, 1) = U_lin;
    RUIO(j).X(:, 1) = x0;
    UIOO(j).X(:, 1) = x0;
end
              
%% Simulation
for FT = 2:-1:1 % 1 - FT is on; 2 - FT is off
    
    % Vector initialization for plots
    FTCS(FT).Y_sim = C*x0; FTCS(FT).Umin = umin; FTCS(FT).Umax = umax;
    
    for k = 1:Nsim
        tk = k*Ts; % Simulation time
        
        if FT == 1
            if k == 1
                YY = FTCS(FT).X(:, k);
                UF = FTCS(FT).Uff(:, k);
            else
                YY = [FTCS(FT).Yfail(1, k-1)-UIOO(1).Fsen(k); FTCS(FT).Yfail(2, k-1)-UIOO(2).Fsen(k); FTCS(FT).Yfail(3, k-1)];
                UF = FTCS(FT).Uff(:, k);
            end
        else
            if k == 1
                YY = FTCS(FT).X(:, k);
                UF = FTCS(FT).Uff(:, k);
            else
                YY = FTCS(FT).Yfail(:, k-1);%X(:, k);%
                UF = FTCS(FT).Uff(:, k);
            end       
        end
        
        % Run MPC controller
        [sol, diag] = mpc{YY, xsp, UF};
        if diag
            msg = ['Infeasible problem at t = ', num2str(k*Ts)];
            disp(msg)
            return;
        end
        FTCS(FT).U(:, k) = sol{1}; FTCS(FT).Obj(:, k) = sol{2};

        %% Actuator fault income
%         if k*Ts > 5 && k*Ts < 15
%             Ufail(:, k) = U(:, k) + f1;
%             FTCS(FTC).Umax(:, k) = FTCS(FTC).Umax + f1;
%             FTCS(FTC).Umin(:, k) = FTCS(FTC).Umin + f1;
%         elseif k*Ts > 20 && k*Ts < 25
%             Ufail(:, k) = U(:, k) + f2;
%             FTCS(FTC).Umax(:, k) = FTCS(FTC).Umax + f2;
%             FTCS(FTC).Umin(:, k) = FTCS(FTC).Umin + f2;
%         else
%             Ufail(:, k) = U(:, k);
%         end

        FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k);      
        FTCS(FT).Ufails(:, k) = [0; 0];
        FTCS(FT).Umax(:, k) = umax;
        FTCS(FT).Umin(:, k) = umin;
 
        if tk > 5 && tk < 15
            FTCS(FT).Ufails(:, k) = [Fail_Q1; 0];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
        end

        if tk > 20 && tk < 25
            FTCS(FT).Ufails(:, k) = [0; -Fail_Q2+Fail_Q2*(exp(-2*(tk-20)/1))];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
        end

        % Natural system saturation
        for j = 1:nu
            if FTCS(FT).Ufail(j, k) >= FTCS(FT).Umax(j, k)
                FTCS(FT).Ufail(j, k) = FTCS(FT).Umax(j, k);
            elseif FTCS(FT).Ufail(j, k) <= FTCS(FT).Umin(j, k)
                FTCS(FT).Ufail(j, k) = FTCS(FT).Umin(j, k);
            end
        end
        
        % Continuous-time simulation (reality)
        %TODO: Evaluate use ODE tool
        Dt  = linspace(0, Ts, 10)';
        u = (FTCS(FT).Ufail(:, k) - U_lin)*ones(1, numel(Dt)); % Zero-order hold input (with Input offset)
        X_sim = lsim(sys, u, Dt, FTCS(FT).X(:, k)-X_lin);   % Nominal continuous system
        X_sim = X_sim' + X_lin;                        % (with state offset)
        FTCS(FT).Y_sim = [FTCS(FT).Y_sim C*X_sim];                    % System output with offset (Continuous time)
        FTCS(FT).X(:, k+1) = X_sim(:, end);                     % Discrete-time state
        
        % Natural state limits (maybe should be erased)
        for j = 1:nu
            if FTCS(FT).X(j, k) >= xmax(j)
                FTCS(FT).X(j, k) = xmax(j);
            elseif FTCS(FT).X(j, k) <= xmin(j)
                FTCS(FT).X(j, k) = xmin(j);
            end
        end
        
        FTCS(FT).Y(:, k) = C*FTCS(FT).X(:, k);                                  % Discrete-time output

        %% Sensor fault income
        %TODO: Add output fails
        FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k);
        
%         if k*Ts > 28 && k*Ts < 38
%             Yfail(:, k) = Y(:, k) + [0; -3.5; 0];
%         end
% 
%         if k*Ts >40 && k*Ts < 50
%             Yfail(:, k) = Y(:, k) + [3; 0; 0];
%         end

        if tk > 28 && tk < 38
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; Fail_S2-Fail_S2*(exp(-5*(tk-28)/1)); 0];
        end

        if tk >40 && tk < 51
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [Fail_S1-Fail_S1*(exp(-4*(tk-40)/1)); 0; 0];
        end
        
        %% RUIO 1
        RUIO(1).Phi(:, k+1) = RUIO(1).K*RUIO(1).Phi(:, k) + RUIO(1).L_ast*FTCS(FT).Yfail(:, k) + RUIO(1).B_bar_1*FTCS(FT).U(:, k) + RUIO(1).delta_bar_1;
        RUIO(1).X(:, k) = RUIO(1).T*[RUIO(1).Phi(:, k); RUIO(1).U_1*FTCS(FT).Yfail(:, k)-RUIO(1).U_1*RUIO(1).C_tilde_1*RUIO(1).Phi(:, k)];

        RUIO(1).Fact(k) = RUIO(1).U_1*(FTCS(FT).X(:, k+1) - RUIO(1).C_tilde_1*RUIO(1).Phi(:, k+1)) + RUIO(1).A_bar_22*RUIO(1).U_1*(RUIO(1).C_tilde_1*RUIO(1).Phi(:, k) - FTCS(FT).Yfail(:, k)) - RUIO(1).A_bar_21*RUIO(1).Phi(:, k) - RUIO(1).B_bar_2*FTCS(FT).U(:, k) - RUIO(1).delta_bar_2;

        % Error norm 1
        RUIO(1).error(k) = sqrt((RUIO(1).X(1, k)-FTCS(FT).Yfail(1, k))^2 + (RUIO(1).X(2, k)-FTCS(FT).Yfail(2, k))^2 + (RUIO(1).X(3, k)-FTCS(FT).Yfail(3, k))^2);

        if RUIO(1).error(k) > threshold(1, k)
            RUIO(1).FQ(k) = true;
        else
            RUIO(1).FQ(k) = false;
        end        

        %% RUIO 2
        RUIO(2).Phi(:, k+1) = RUIO(2).K*RUIO(2).Phi(:, k) + RUIO(2).L_ast*FTCS(FT).Yfail(:, k) + RUIO(2).B_bar_1*FTCS(FT).U(:, k) + RUIO(2).delta_bar_1;
        RUIO(2).X(:, k) = RUIO(2).T*[RUIO(2).Phi(:, k); RUIO(2).U_1*FTCS(FT).Yfail(:, k)-RUIO(2).U_1*RUIO(2).C_tilde_1*RUIO(2).Phi(:, k)];

        RUIO(2).Fact(k) = RUIO(2).U_1*(FTCS(FT).X(:, k+1) - RUIO(2).C_tilde_1*RUIO(2).Phi(:, k+1)) + RUIO(2).A_bar_22*RUIO(2).U_1*(RUIO(2).C_tilde_1*RUIO(2).Phi(:, k) - FTCS(FT).Yfail(:, k)) - RUIO(2).A_bar_21*RUIO(2).Phi(:, k) - RUIO(2).B_bar_2*FTCS(FT).U(:, k) - RUIO(2).delta_bar_2;

        % Error norm 2
        RUIO(2).error(k) = sqrt((RUIO(2).X(1, k)-FTCS(FT).Yfail(1, k))^2 + (RUIO(2).X(2, k)-FTCS(FT).Yfail(2, k))^2 + (RUIO(2).X(3, k)-FTCS(FT).Yfail(3, k))^2);
        
        if RUIO(2).error(k) > threshold(2, k)
            RUIO(2).FQ(k) = true;
        else
            RUIO(2).FQ(k) = false;
        end
        
        %% UIOO 1
        UIOO(1).Ymon(:, k) = UIOO(1).T2*FTCS(FT).Yfail(:, k);
        UIOO(1).Z(:, k+1) = zeros(nx, 1);      

        UIOO(1).Z(:, k+1) = UIOO(1).N*UIOO(1).Z(:, k) + UIOO(1).L*UIOO(1).Ymon(:, k) + UIOO(1).G*FTCS(FT).U(:, k) + UIOO(1).Tg;

        UIOO(1).X(:, k) = UIOO(1).Z(:, k) - UIOO(1).E*UIOO(1).Ymon(:, k);
        UIOO(1).Y(:, k) = Cd*UIOO(1).X(:, k);

        % Residue 1
        UIOO(1).res(:, k) = FTCS(FT).Yfail(:, k) - UIOO(1).Y(:, k);

        % Error norm 1
        UIOO(1).error(k) = sqrt(UIOO(1).res(1, k)^2);

        if UIOO(1).error(k) > threshold(3, k)
            UIOO(1).FO(k) = true;
        else
            UIOO(1).FO(k) = false;
        end        

        %% UIOO 2
        UIOO(2).Ymon(:, k) = UIOO(2).T2*FTCS(FT).Yfail(:, k);
        UIOO(2).Z(:, k+1) = zeros(nx, 1);      

        UIOO(2).Z(:, k+1) = UIOO(2).N*UIOO(2).Z(:, k) + UIOO(2).L*UIOO(2).Ymon(:, k) + UIOO(2).G*FTCS(FT).U(:, k) + UIOO(2).Tg;

        UIOO(2).X(:, k) = UIOO(2).Z(:, k) - UIOO(2).E*UIOO(2).Ymon(:, k);
        UIOO(2).Y(:, k) = Cd*UIOO(2).X(:, k);

        % Residue 2
        UIOO(2).res(:, k) = FTCS(FT).Yfail(:, k) - UIOO(2).Y(:, k);

        % Error norm 2
        UIOO(2).error(k) = sqrt(UIOO(2).res(2, k)^2);

        if UIOO(2).error(k) > threshold(4, k)
            UIOO(2).FO(k) = true;
        else
            UIOO(2).FO(k) = false;
        end
        
        %% Actuator fault estimation
        if RUIO(1).FQ(k) && ~RUIO(2).FQ(k) && ~UIOO(1).FO(k) && UIOO(2).FO(k) % Actuator fault 1
            if RUIO(2).delay
                RUIO(2).Fact(k) = RUIO(2).Fact(k);
            else
                RUIO(2).delay = 1;
                RUIO(2).Fact(k) = 0;
            end
        elseif ~RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k) % Actuator fault 2
            if RUIO(1).delay
                RUIO(1).Fact(k) = RUIO(1).Fact(k);
            else
                RUIO(1).delay = 1;
                RUIO(1).Fact(k) = 0;
            end
        else
            RUIO(1).delay = 0;
            RUIO(2).delay = 0;
            RUIO(1).Fact(k) = 0;
            RUIO(2).Fact(k) = 0;
        end

        %% Sensor fault estimation
        % Sensor fault 1
        if RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k)
            UIOO(1).Fsen(k) = UIOO(2).res(1, k);
        else
            UIOO(1).Fsen(k) = zeros(size(UIOO(2).res(1, k)));
        end

        % Sensor fault 2
        if RUIO(1).FQ(k) && RUIO(2).FQ(k) && ~UIOO(1).FO(k) && UIOO(2).FO(k)
            UIOO(2).Fsen(k) = UIOO(1).res(2, k);
        else
            UIOO(2).Fsen(k) = zeros(size(UIOO(1).res(2, k)));
        end
        
        % If FT-MPC is enabled
        if FT == 1
            FTCS(FT).Uff(:, k+1) = [RUIO(1).Fact(k); RUIO(2).Fact(k)];
        else
            FTCS(FT).Uff(:, k+1) = [0; 0];
        end
        
        FTCS(FT).RUIO(1).error(k) = RUIO(1).error(k);
        FTCS(FT).RUIO(1).Fact(k) = RUIO(1).Fact(k);
        FTCS(FT).RUIO(2).error(k) = RUIO(2).error(k);
        FTCS(FT).RUIO(2).Fact(k) = RUIO(2).Fact(k);
        
        FTCS(FT).UIOO(1).error(k) = UIOO(1).error(k);
        FTCS(FT).UIOO(1).Fsen(k) = UIOO(1).Fsen(k);
        FTCS(FT).UIOO(2).error(k) = UIOO(2).error(k);
        FTCS(FT).UIOO(2).Fsen(k) = UIOO(2).Fsen(k);
        
    end
end

save FTCS.mat

run enPlotHE
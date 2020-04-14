%% MPC on CSTR
clc; clear; yalmip('clear'); close all;

% When generates flat figures
set(0, 'DefaultFigureRenderer', 'painters');

%% Load polytope and observer matrices
% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
RUIO = struct;
UIOO = struct;
load linObs

%% LTI model and observers
% This section is commented to reduce simulation time (using pre-calculated observer matrices)
% Ts = 0.05;                       % Sample time [min]
% Vr = 100;                        % Reactor volume [l]
% Tr = 440;                         % Temperature [K]
% 
% run CSTR;
% 
% %% Reduced-order unknown input observer
% N = 2;
% run CSTR_RUIO;
% 
% %% Unknown input output observer
% run CSTR_UIOO;
% 
% % Save observers' data
% save linObs.mat

%% Simulation parameters
Time = 60;                              % Simulation end time 
Nsim = Time/Ts;                     % Simulation steps
t = 0:Ts:Time-Ts;                    % Simulation time

Fail_Q1 = -5; Fail_Q2 = 5;      % Actuator fault magnitude [5%, 5%]
Fail_S1 = 2; Fail_S3 = -1.5;    % Sensor fault magnitude [2% 0.5%]
x0 = [95; 0.105; 437];            % Start-point
xsp = [Vr; Ca; Tr];                   % Set-point

%% MPC controller
% Constraints
xmin = [90; 0.07; 434];
xmax = [105; 0.11; 444];
umin = [90; 90];
umax = [110; 110];

% Wheight matrix
Qx = eye(nx);
% Qx = diag([1 1 1]);
Rx = eye(nu);
% Rx = diag([1 1]);
gamma = 1e2*diag([1 1e4 1e-1]);

% Controller
Np = 5;                                   % Prediction horizon
run MPC;

%% Noise
sig = 5e-15*([1 5 2])';         % Ouput noise sigma

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
mag_1 = 2e-2;     % Value Q1
mag_2 = 2e-4;     % Value Q2
mag_3 = 3e-6;     % Value O1
mag_4 = 1e-3;     % Value O2

threshold = zeros(4, Nsim);

for k = 1:Nsim
    threshold(1, k) = mag_1 + 500*exp(-(k-1)/Tau);  % Q1
    threshold(2, k) = mag_2 + 100*exp(-(k-1)/Tau);    % Q2
    threshold(3, k) = mag_3 + 0.6*exp(-(k-1)/Tau);    % O1
    threshold(4, k) = mag_4 + 400*exp(-(k-1)/Tau);  % O2
end

%% Simulation Setup
% Fault Tolerant Control System (1 = disable; 2 = enable).Matrix 
FTCS = struct;

%% Simulation
disp('Simulating...')
FTC_OFF = 1; FTC_ON = 2;
for FT = 1:2    % 1 - FT is off; 2 -  FT is on
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
    
    % FTCS matrices
    FTCS(FT).U = zeros(nu, Nsim);                   % Control Input
    FTCS(FT).Ufail = zeros(nu, Nsim);              % Faulty control Input
    FTCS(FT).Uff = zeros(nu, Nsim+1);            % Feedforward control Input
    FTCS(FT).Ufails = zeros(nu, Nsim);             % Fails of control inputs
    FTCS(FT).X = zeros(nx, Nsim+1);               % States
    FTCS(FT).Y = zeros(ny, Nsim);                    % Measure outputs
    FTCS(FT).Y_hat = zeros(ny, Nsim);             % Estimated outputs
    FTCS(FT).Yfail = zeros(ny, Nsim);               % Faulty measure outputs
    FTCS(FT).Obj = zeros(1, Nsim);                  % Objective cost

    % Initial states and inputs
    FTCS(FT).X(:, 1) = x0;
    FTCS(FT).Y(:, 1) = C*x0;
    FTCS(FT).Y_hat(:, 1) = C*x0;
    FTCS(FT).U(:, 1) = U_lin;
    RUIO(FT).X(:, 1) = x0;
    UIOO(FT).X(:, 1) = x0;
    
    % Vector initialization for plots
    FTCS(FT).Yc_sim = C*x0;
    
    if FT == FTC_ON
        disp('Fault tolerant = ON')
    else
        disp('Fault tolerant = OFF')
    end

    for k = 1:Nsim
        tk = k*Ts; % Simulation time
        
        % Run MPC controller
        [sol, diag] = mpc{FTCS(FT).Y_hat(:, k), xsp, FTCS(FT).Uff(:, k)};
        if diag
            msg = ['Infeasible problem at t = ', num2str(k*Ts)];
            disp(msg)
            return;
        end
        FTCS(FT).U(:, k) = sol{1}; FTCS(FT).Obj(:, k) = sol{2};

        %% Actuator fault income
        % No fault
        FTCS(FT).Ufails(:, k) = [0; 0];
        FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k);      
        FTCS(FT).Umax(:, k) = umax;
        FTCS(FT).Umin(:, k) = umin;
 
        % Q1 fault
        if tk > 5 && tk < 15
            FTCS(FT).Ufails(:, k) = [Fail_Q1; 0];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
            % TODO: Check if the constraints may be modified
%             FTCS(FT).Umax(:, k) = umax + FTCS(FT).Ufails(:, k);
%             FTCS(FT).Umin(:, k) = umin + FTCS(FT).Ufails(:, k);
        end

        % Q2 fault
        if tk > 50 && tk < 60
            FTCS(FT).Ufails(:, k) = [0; +Fail_Q2-Fail_Q2*(exp(-2*(tk-50)/1))];
            FTCS(FT).Ufail(:, k) = FTCS(FT).U(:, k) + FTCS(FT).Ufails(:, k);
%             FTCS(FT).Umax(:, k) = umax + FTCS(FT).Ufails(:, k);
%             FTCS(FT).Umin(:, k) = umin + FTCS(FT).Ufails(:, k);
        elseif tk >= 60 && tk < 62
            FTCS(FT).Ufails(:, k) = [0; Fail_Q2*(exp(-8*(tk-60)/1))];
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
        Xc_sim = lsim(sys, u, Dt, FTCS(FT).X(:, k)-X_lin);   % Nominal continuous system
        Xc_sim = Xc_sim' + X_lin;                                      % (with state offset)
        FTCS(FT).Yc_sim = [FTCS(FT).Yc_sim C*Xc_sim];  % System output with offset (Continuous time)
        FTCS(FT).X(:, k+1) = Xc_sim(:, end) + v(:, k);       % Discrete-time state + noise

        % Natural state limits (maybe should be erased)
        for j = 1:nu
            if FTCS(FT).X(j, k) >= xmax(j)
                FTCS(FT).X(j, k) = xmax(j);
            elseif FTCS(FT).X(j, k) <= xmin(j)
                FTCS(FT).X(j, k) = xmin(j);
            end
        end
        
        FTCS(FT).Y(:, k) = C*FTCS(FT).X(:, k);                   % Discrete-time output

        %% Sensor fault income
        FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k);
        
        if tk > 20 && tk < 30
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; 0; Fail_S3-Fail_S3*(exp(-3*(tk-20)/1))];
        elseif tk >= 30 && tk < 32
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [0; 0; Fail_S3*(exp(-5*(tk-30)/1))];
        end

        if tk > 35 && tk < 45
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [Fail_S1-Fail_S1*(exp(-3*(tk-35)/1)); 0; 0];
        elseif tk >= 45 && tk < 47
            FTCS(FT).Yfail(:, k) = FTCS(FT).Y(:, k) + [Fail_S1*(exp(-7*(tk-45)/1)); 0; 0];
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
        UIOO(2).error(k) = sqrt(UIOO(2).res(3, k)^2);

        if UIOO(2).error(k) > threshold(4, k)
            UIOO(2).FO(k) = true;
        else
            UIOO(2).FO(k) = false;
        end
        
        %% Actuator fault estimation
        if RUIO(1).FQ(k) && ~RUIO(2).FQ(k) && UIOO(1).FO(k) && UIOO(2).FO(k) % Actuator fault 1
            if RUIO(2).delay > 1
                RUIO(2).Fact(k) = RUIO(2).Fact(k);
            else
                RUIO(2).delay = RUIO(2).delay + 1;
                RUIO(2).Fact(k) = 0;
            end
        else
            RUIO(2).delay = 0;
            RUIO(2).Fact(k) = 0;
        end
        if ~RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k) % Actuator fault 2
            if RUIO(1).delay > 1
                RUIO(1).Fact(k) = RUIO(1).Fact(k);
            else
                RUIO(1).delay = RUIO(1).delay + 1;
                RUIO(1).Fact(k) = 0;
            end
        else
            RUIO(1).delay = 0;
            RUIO(1).Fact(k) = 0;       
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
            UIOO(2).Fsen(k) = UIOO(1).res(3, k);
        else
            UIOO(2).Fsen(k) = zeros(size(UIOO(1).res(3, k)));
        end
        
        % If FT-MPC is enabled
        if FT == FTC_ON
            FTCS(FT).Uff(:, k+1) = [RUIO(1).Fact(k); RUIO(2).Fact(k)];
            FTCS(FT).Y_hat(:, k+1) = [FTCS(FT).Yfail(1, k)-UIOO(1).Fsen(k); FTCS(FT).Yfail(2, k); FTCS(FT).Yfail(3, k)-UIOO(2).Fsen(k)];
        else
            FTCS(FT).Uff(:, k+1) = [0; 0];
            FTCS(FT).Y_hat(:, k+1) = FTCS(FT).Yfail(:, k);
        end
        
%         FTCS(FT).Uff(:, k+1) = [0; 0];
%         FTCS(FT).Y_hat(:, k+1) = FTCS(FT).Yfail(:, k);
      
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

yalmip('clear')
disp('Saving...')
save FTCS.mat

run enPlotCSTR
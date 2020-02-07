%% MPC on HE
clc; clear; yalmip('clear');
close all;

UIOO = struct;
RUIO = struct;

% When generates flat figures
set(0, 'DefaultFigureRenderer', 'painters');

%% LTI model
Ts = 0.05;               % Sample time [min]
Theta_1s = 498;     % Temperatura de salida de fluido 1 (K) (set-point) (498)
Theta_2s = 690;     % Temperatura de salida de fluido 2 (K)

run HE;

%% Parameters
Np = 5;                                   % Prediction horizon
Time = 50;                              % Simulation end time 
Nsim = Time/Ts;                     % Simulation steps
t = 0:Ts:Time-Ts;                    % Simulation time

x0 = [495; 680; 570];             % Start-point
f1 = [0; 1]; f2 = [-8; 0];          % Fault magnitude
xsp = [Theta_1s; Theta_2s; Theta_p];% Set-point

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

%% Reduced-order unknown input observer
N = 2;
run HE_RUIO;

%% Unknown input output observer
run HE_UIOO;

%% Noise
sig = 3e-3*([1 1 1])';         % Ouput noise sigma

rng default;                        % Random seed start
v = sig*randn(1, Nsim);    % Measurement noise v~N(0, sig)

%% Error detection threshold
Tau = 2;                % Convergence period
mag_1 = 1e-1;     % Value Q1
mag_2 = 2e-1;     % Value Q2
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
U = zeros(nu, Nsim);                   % Control Input
Ufail = zeros(nu, Nsim);              % Faulty control Input
Uff = zeros(nu, Nsim+1);            % Feedforward control Input
X = zeros(nx, Nsim+1);               % States
Y = zeros(ny, Nsim);                    % Measure outputs
Yfail = zeros(ny, Nsim);                % Faulty measure outputs
Obj = zeros(1, Nsim);                   % Objective cost

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
X(:, 1) = x0;
Y(:, 1) = C*x0;
U(:, 1) = U_lin;
RUIO(1).X(:, 1) = x0;
RUIO(2).X(:, 1) = x0;
UIOO(1).X(:, 1) = x0;
UIOO(2).X(:, 1) = x0;
              
%% Simulation
for FTC = 0:1 % 0 - FTC is off; 1 - FTC is on
    
    % Vector initialization for plots
    Y_sim = C*x0; Umin = umin; Umax = umax;
    
    for k = 1:Nsim
        tk = k*Ts; % Simulation time
        
        % Run MPC controller
        [sol, diag] = mpc{X(:, k), xsp, Uff(:, k)};
        if diag
            msg = ['Infeasible problem at t = ', num2str(k*Ts)];
            disp(msg)
            return;
        end
        U(:, k) = sol{1}; Obj(:, k) = sol{2};

        %% Actuator fault income
        if k*Ts > 5 && k*Ts < 15
            Ufail(:, k) = U(:, k) + f1;
            Umax(:, k) = umax + f1;
            Umin(:, k) = umin + f1;
        elseif k*Ts > 20 && k*Ts < 25
            Ufail(:, k) = U(:, k) + f2;
            Umax(:, k) = umax + f2;
            Umin(:, k) = umin + f2;
        else
            Ufail(:, k) = U(:, k);
            Umax(:, k) = umax;
            Umin(:, k) = umin;
        end

        % Natural system saturation
        for j = 1:nu
            if Ufail(j, k) >= Umax(j, k)
                Ufail(j, k) = Umax(j, k);
            elseif Ufail(j, k) <= Umin(j, k)
                Ufail(j, k) = Umin(j, k);
            end
        end
        
        % Continuous-time simulation (reality)
        %TODO: Evaluate use ODE tool
        Dt  = linspace(0, Ts, 10)';
        u = (Ufail(:, k) - U_lin)*ones(1, numel(Dt)); % Zero-order hold input (with Input offset)
        X_sim = lsim(sys, u, Dt, X(:, k)-X_lin);   % Nominal continuous system
        X_sim = X_sim' + X_lin;                        % (with state offset)
        Y_sim = [Y_sim C*X_sim];                    % System output with offset (Continuous time)
        X(:, k+1) = X_sim(:, end);                     % Discrete-time state
        
        % Natural state limits (maybe should be erased)
        for j = 1:nu
            if X(j, k) >= xmax(j)
                X(j, k) = xmax(j);
            elseif X(j, k) <= xmin(j)
                X(j, k) = xmin(j);
            end
        end
        
        Y(:, k) = C*X(:, k);                                  % Discrete-time output

        %% Sensor fault income
        %TODO: Add output fails
        Yfail(:, k) = Y(:, k);
        
%         if k*Ts > 28 && k*Ts < 38
%             Yfail(:, k) = Y(:, k) + [0; -1; 0];
%         end
% 
%         if k*Ts >40 && k*Ts < 50
%             Yfail(:, k) = Y(:, k) + [1; 0; 0];
%         end
        
        %% RUIO 1
        RUIO(1).Phi(:, k+1) = RUIO(1).K*RUIO(1).Phi(:, k) + RUIO(1).L_ast*Yfail(:, k) + RUIO(1).B_bar_1*U(:, k) + RUIO(1).delta_bar_1;
        RUIO(1).X(:, k) = RUIO(1).T*[RUIO(1).Phi(:, k); RUIO(1).U_1*Yfail(:, k)-RUIO(1).U_1*RUIO(1).C_tilde_1*RUIO(1).Phi(:, k)];

        RUIO(1).Fact(k) = RUIO(1).U_1*(X(:, k+1) - RUIO(1).C_tilde_1*RUIO(1).Phi(:, k+1)) + RUIO(1).A_bar_22*RUIO(1).U_1*(RUIO(1).C_tilde_1*RUIO(1).Phi(:, k) - Yfail(:, k)) - RUIO(1).A_bar_21*RUIO(1).Phi(:, k) - RUIO(1).B_bar_2*U(:, k) - RUIO(1).delta_bar_2;

        % Error norm 1
        RUIO(1).error(k) = sqrt((RUIO(1).X(1, k)-Yfail(1, k))^2 + (RUIO(1).X(2, k)-Yfail(2, k))^2 + (RUIO(1).X(3, k)-Yfail(3, k))^2);
        
        if RUIO(1).error(k) > threshold(1, k)
            RUIO(1).FQ(k) = true;
        else
            RUIO(1).FQ(k) = false;
        end        

        %% RUIO 2
        RUIO(2).Phi(:, k+1) = RUIO(2).K*RUIO(2).Phi(:, k) + RUIO(2).L_ast*Yfail(:, k) + RUIO(2).B_bar_1*U(:, k) + RUIO(2).delta_bar_1;
        RUIO(2).X(:, k) = RUIO(2).T*[RUIO(2).Phi(:, k); RUIO(2).U_1*Yfail(:, k)-RUIO(2).U_1*RUIO(2).C_tilde_1*RUIO(2).Phi(:, k)];

        RUIO(2).Fact(k) = RUIO(2).U_1*(X(:, k+1) - RUIO(2).C_tilde_1*RUIO(2).Phi(:, k+1)) + RUIO(2).A_bar_22*RUIO(2).U_1*(RUIO(2).C_tilde_1*RUIO(2).Phi(:, k) - Yfail(:, k)) - RUIO(2).A_bar_21*RUIO(2).Phi(:, k) - RUIO(2).B_bar_2*U(:, k) - RUIO(2).delta_bar_2;

        % Error norm 2
        RUIO(2).error(k) = sqrt((RUIO(2).X(1, k)-Yfail(1, k))^2 + (RUIO(2).X(2, k)-Yfail(2, k))^2 + (RUIO(2).X(3, k)-Yfail(3, k))^2);
        
        if RUIO(2).error(k) > threshold(2, k)
            RUIO(2).FQ(k) = true;
        else
            RUIO(2).FQ(k) = false;
        end
        
        %% UIOO 1
        UIOO(1).Ymon(:, k) = UIOO(1).T2*Yfail(:, k);
        UIOO(1).Z(:, k+1) = zeros(nx, 1);      

        UIOO(1).Z(:, k+1) = UIOO(1).N*UIOO(1).Z(:, k) + UIOO(1).L*UIOO(1).Ymon(:, k) + UIOO(1).G*U(:, k) + UIOO(1).Tg;

        UIOO(1).X(:, k) = UIOO(1).Z(:, k) - UIOO(1).E*UIOO(1).Ymon(:, k);
        UIOO(1).Y(:, k) = Cd*UIOO(1).X(:, k);

        % Residue 1
        UIOO(1).res(:, k) = Yfail(:, k) - UIOO(1).Y(:, k);

        % Error norm 1
        UIOO(1).error(k) = sqrt(UIOO(1).res(1, k)^2);

        if UIOO(1).error(k) > threshold(3, k)
            UIOO(1).FO(k) = true;
        else
            UIOO(1).FO(k) = false;
        end        

        %% UIOO 2
        UIOO(2).Ymon(:, k) = UIOO(2).T2*Yfail(:, k);
        UIOO(2).Z(:, k+1) = zeros(nx, 1);      

        UIOO(2).Z(:, k+1) = UIOO(2).N*UIOO(2).Z(:, k) + UIOO(2).L*UIOO(2).Ymon(:, k) + UIOO(2).G*U(:, k) + UIOO(2).Tg;

        UIOO(2).X(:, k) = UIOO(2).Z(:, k) - UIOO(2).E*UIOO(2).Ymon(:, k);
        UIOO(2).Y(:, k) = Cd*UIOO(2).X(:, k);

        % Residue 2
        UIOO(2).res(:, k) = Yfail(:, k) - UIOO(2).Y(:, k);

        % Error norm 2
        UIOO(2).error(k) = sqrt(UIOO(2).res(2, k)^2);

        if UIOO(2).error(k) > threshold(4, k)
            UIOO(2).FO(k) = true;
        else
            UIOO(2).FO(k) = false;
        end
        
        %% Actuator fault estimation
        if RUIO(1).FQ(k) && ~RUIO(2).FQ(k)% && ~UIOO(1).FO(k) && UIOO(2).FO(k) % Actuator fault 1
            if RUIO(2).delay
                RUIO(2).Fact(k) = RUIO(2).Fact(k);
            else
                RUIO(2).delay = 1;
                RUIO(2).Fact(k) = 0;
            end
        elseif ~RUIO(1).FQ(k) && RUIO(2).FQ(k)% && UIOO(1).FO(k) && ~UIOO(2).FO(k) % Actuator fault 2
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

%         %% Next step fault compensation (feedforward)
%         %TODO: Add exponential threshold
%         if  k*Ts>= 2 % Estimation is ok
%             if RUIO(1).error(k) < threshold(1, k)  % Error 1 threshold
%                 RUIO(1).FQ(k) = true;
%                 RUIO(2).Fact(k) = 0;
%             else
%                 RUIO(1).FQ(k) = false;
%             end
%             if RUIO(2).error(k) < threshold(2, k)  % Error 2 threshold
%                 RUIO(2).FQ(k) = true;
%                 RUIO(1).Fact(k) = 0;
%             else
%                 RUIO(2).FQ(k) = false;
%             end
%         else
%             RUIO(1).Fact(k) = 0; RUIO(2).Fact(k) = 0;
%         end
% 
%         %% Sensor fault estimation
%         % Sensor fault 1
%         if RUIO(1).FQ(k) && RUIO(2).FQ(k) && ~UIOO(1).FO(k) && UIOO(2).FO(k)
%             UIOO(1).Fsen(k) = UIOO(2).res(1, k);
%         else
%             UIOO(1).Fsen(k) = zeros(size(UIOO(2).res(1, k)));
%         end
% 
%         % Sensor fault 2
%         if RUIO(1).FQ(k) && RUIO(2).FQ(k) && UIOO(1).FO(k) && ~UIOO(2).FO(k)
%             UIOO(2).Fsen(k) = UIOO(1).res(2, k);
%         else
%             UIOO(2).Fsen(k) = zeros(size(UIOO(1).res(2, k)));
%         end
        
        % If FT-MPC is enabled
        if FTC == 1
            Uff(:, k+1) = [RUIO(1).Fact(k); RUIO(2).Fact(k)];
        else
            Uff(:, k+1) = [0; 0];
        end
        
    end

    %% Set Plots
    vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];
    tc = 0:Ts/10:Nsim/20;
    
    % Outputs
    figure(1)
    if FTC == 0
        subplot(311)
%         plot(t, state(1, :), 'b.', tc, y(1, :), 'y-.', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
        plot(t, xsp(1)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
        plot(tc, Y_sim(1, :), 'b', 'LineWidth', 1.5);
        % plot(t, xmin(1)*ones(length(t)), 'y--')
        % plot(t, xmax(1)*ones(length(t)), 'y--')
        % plot(t, state(1, :), 'b-.', t, x1_hat(1, :), 'g:', t, x2_hat(1, :), 'y:', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
        xlabel('Time [min]'); ylabel('\theta_1 [K]'); grid on
        subplot(312)
%         plot(t, state(2, :), 'b.', tc, y(2, :), 'y-.', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(2)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
        h(2) = plot(tc, Y_sim(2, :), 'b', 'LineWidth', 1.5);
        % plot(t, xmin(2)*ones(length(t)), 'y--')
        % plot(t, xmax(2)*ones(length(t)), 'y--')
        % plot(t, state(2, :), 'b-.', t, x1_hat(2, :), 'g:', t, x2_hat(2, :), 'y:', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_2 [K]'); grid on
        subplot(313)
%         plot(t, state(3, :), 'b.', tc, y(3, :), 'y-.', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(3)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
        plot(tc, Y_sim(3, :), 'b', 'LineWidth', 1.5);   
        % plot(t, xmin(3)*ones(length(t)), 'y--')
        % plot(t, xmax(3)*ones(length(t)), 'y--')
        % plot(t, state(3, :), 'b-.', t, x1_hat(3, :), 'g:', t, x2_hat(3, :), 'y:', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_p [K]'); grid on
    else
        subplot(311)
%         plot(t, state(1, :), 'k.', tc, y(1, :), 'g-.', 'LineWidth', 1.5);
        plot(tc, Y_sim(1, :), 'k-.', 'LineWidth', 1.5);
        hold off
        legend('x_s', 'Location', 'SouthEast');
        legend boxoff
        subplot(312)
%         plot(t, state(2, :), 'k.', tc, y(2, :), 'g-.', 'LineWidth', 1.5);
        plot(tc, Y_sim(2, :), 'k-.', 'LineWidth', 1.5);        
        hold off
        legend(h(2), 'MPC', 'Location', 'SouthEast');
        legend boxoff
        subplot(313)
%         plot(t, state(3, :), 'k.', tc, y(3, :), 'g-.', 'LineWidth', 1.5);
        h(3) = plot(tc, Y_sim(3, :), 'k-.', 'LineWidth', 1.5);
        hold off
        legend(h(3), 'FTMPC', 'Location', 'NorthEast');
        legend boxoff
        print -dsvg figs/outputHE.svg
    end

% 	%% Manipulated variables
%     figure(2)
%     if FTC == 0
%         subplot(211)
%         stairs(t, U(1, 1:end-1), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_1 [l/min]'); grid on
%         subplot(212)
%         stairs(t, U(2, 1:end-1), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_2 [l/min]'); grid on
% 	else
%         subplot(211)
%         stairs(t, U(1, 1:end-1), 'k-.', 'LineWidth', 1.5)
%         plot(t, umin(1)*ones(length(t)), 'r--')
%         plot(t, umax(1)*ones(length(t)), 'r--')
%         hold off
%         legend('MPC', 'FTMPC', 'Location', 'NorthWest');
%         legend boxoff
%         subplot(212)
%         stairs(t, U(2, 1:end-1), 'k-.', 'LineWidth', 1.5)
%         plot(t, umin(2)*ones(length(t)), 'r--')
%         plot(t, umax(2)*ones(length(t)), 'r--')
%         hold off
%         print -dsvg figs/inputHE.svg
%     end
% 
%     % Failure inputs
%     figure(3)
%     if FTC == 0
%         subplot(211)
%         plot(t, Umin(1, :), 'r--')
%         hold on
%         plot(t, Umax(1, :), 'r--')
%         stairs(t, Ufail(1, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
%         subplot(212)
%         plot(t, Umin(2, :), 'r--')
%         hold on
%         plot(t, Umax(2, :), 'r--')
%         stairs(t, Ufail(2, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
% 	else
%         subplot(211)
%         stairs(t, Ufail(1, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, Ufail(2, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%         print -dsvg figs/inputfailHE.svg
%     end

    % RUIO error detection
    figure(4)
    if FTC == 0
        subplot(211)
        plot(t, RUIO(1).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        plot(t, threshold(1, :),  'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('|e_x|');
        axis([0 inf 0 6.5])
        subplot(212)
        plot(t, RUIO(2).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        plot(t, threshold(2, :),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.6])
        xlabel('Time [min]'); ylabel('|e_x|');
	else
        subplot(211)
        plot(t, RUIO(1).error, 'k-.', 'LineWidth', 1.5)
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
        legend boxoff
        subplot(212)
        plot(t, RUIO(2).error, 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/RUIOerrorHE.svg
    end
    
	% UIOO error detection
    figure(5)
    if FTC == 0
        subplot(211)
        plot(t, UIOO(1).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        plot(t, threshold(3, :),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.04])
        xlabel('Time [min]'); ylabel('|e_x|');
        subplot(212)
        plot(t, UIOO(2).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        plot(t, threshold(4, :),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.03])
        xlabel('Time [min]'); ylabel('|e_x|');
	else
        subplot(211)
        plot(t, UIOO(1).error, 'k-.', 'LineWidth', 1.5)
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
        legend boxoff
        subplot(212)
        plot(t, UIOO(2).error, 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/UIOOerrorHE.svg
    end

%     % Fault estimation
%     figure(6)
%     if FTC == 0
%         subplot(211)
%         stairs(t, RUIO(1).Fact, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, Ufail(1, :) - U(1, 1:end-1), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
%         subplot(212)
%         stairs(t, RUIO(2).Fact, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, Ufail(2, :) - U(2, 1:end-1), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
%     else
%         subplot(211)
%         stairs(t, RUIO(1).Fact, 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, RUIO(2).Fact, 'k-.', 'LineWidth', 1.5)
%         hold off
%         print -dsvg figs/estimationHE.svg
%     end
% 
%     % Objective
%     figure(7)
%     if FTC == 0
%         plot(t, Obj, 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('Objective'); grid on
%     else
%         plot(t, Obj, 'k-.', 'LineWidth', 1.5)
%         hold off
%         axis([0 inf 0 400])
%         print -dsvg figs/objectiveHE.svg
%     end
%     
%     % State evolution
%     figure(8)
%     if FTC == 0
%         plot3(x0(1), x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
%         plot3(Y(1, :), Y(2, :), Y(3, :), 'b', 'LineWidth', 1.5)
%         plot3(Y(1, end), Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
%         plot3(xsp(1), xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(Xpoly+X_lin, 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xs+X_lin, 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); zlabel('\theta_3 [K]'); grid on;
%     else
%         plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)        
%         plot3(Y(1, :), Y(2, :), Y(3, :), 'k-.', 'LineWidth', 1.5)
%         plot3(Y(1, end), Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%         print -dsvg figs/stateHE.svg
%     end
%     
%     % State evolution
%     Xx = Xpoly.projection(1:2).minHRep();
%     Xxs = Xs.projection(1:2).minHRep();
% 
%     figure(9)
%     if FTC == 0
%         plot(x0(1), x0(2), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot(Y_sim(1, :), Y_sim(2, :), 'y', 'LineWidth', 1.5)        
%         plot(Y(1, :), Y(2, :), 'b.', 'LineWidth', 1.5)
%         plot(Y(1, end), Y(2, end), 'mo', 'LineWidth', 1.5)
%         plot(xsp(1), xsp(2), 'rp', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(1:2), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xx+X_lin(1:2), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); grid on;
%     else
%         plot(Y_sim(1, :), Y_sim(2, :), 'g', 'LineWidth', 1.5)        
%         plot(Y(1, :), Y(2, :), 'k.', 'LineWidth', 1.5)
%         plot(Y(1, end), Y(2, end), 'ro', 'LineWidth', 1.5)
%         hold off
%         print -dsvg figs/state1-2HE.svg
%     end
%     
%     % State evolution
%     Xx = Xpoly.projection(1:2:3).minHRep();
%     Xxs = Xs.projection(1:2:3).minHRep();
% 
%     figure(10)
%     if FTC == 0
%         plot(x0(1), x0(3), 'gd', 'LineWidth', 1.5);
%         hold on
%         plot(Y_sim(1, :), Y_sim(3, :), 'b', 'LineWidth', 1.5)
%         plot(Y(1, end), Y(3, end), 'mo', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_3 [K]'); grid on;
%     else
%         plot(Y_sim(1, :), Y_sim(3, :), 'k--', 'LineWidth', 1.5)
%         plot(Y(1, end), Y(3, end), 'y*', 'LineWidth', 1.5)
%         plot(xsp(1), xsp(3), 'ro', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(1:2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'LineWidth', 1.5)
%         plot(Xx+X_lin(1:2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'LineWidth', 1.5)
%         hold off
%         leg = legend('$x(0)$', '$x_{MPC}$', '$x_{MPC}(end)$', '$x_{FTMPC}$', '$x_{FTMPC}(end)$', '$x_s$', '$\bf{X}_s$', '$\bf{X}$', 'Location', 'SouthEast');
%         set(leg, 'Interpreter', 'latex');
%         print -dsvg figs/state1-3HE.svg
%     end
%     
%     % State evolution
%     Xx = Xpoly.projection(2:3).minHRep();
%     Xxs = Xs.projection(2:3).minHRep();
% 
%     figure(11)
%     if FTC == 0
%         plot(x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot(Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
%         plot(Y(2, :), Y(3, :), 'b.', 'LineWidth', 1.5)
%         plot(Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
%         plot(xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xx+X_lin(2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_2 [K]'); ylabel('\theta_3 [K]'); grid on;
%     else
%         plot(Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)
%         plot(Y(2, :), Y(3, :), 'k.', 'LineWidth', 1.5)
%         plot(Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%         print -dsvg figs/state2-3HE.svg    
%     end
end
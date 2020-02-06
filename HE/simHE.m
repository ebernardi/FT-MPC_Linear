%% MPC on HE
clc; clear; yalmip('clear');
close all;

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
Time = 25;                              % Simulation end time 
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

        % Run MPC controller       
        [sol, diag] = mpc{X(:,k), xsp, Uff(:, k)};
        if diag
            msg = ['Infeasible problem at t = ', num2str(k*Ts)];
            disp(msg)
            return;
        end
        U(:, k) = sol{1}; Obj(:, k) = sol{2};

        % Input limits
        for j = 1:nu
            if U(j, k) >= umax(j)
                U(j, k) = umax(j);
            elseif U(j, k) <= umin(j)
                U(j, k) = umin(j);
            end
        end

        % Add actuator fails
        if k*Ts > 5 && k*Ts < 15
            Ufail(:, k) = U(:, k) + f1;
            Umax(:, k) = umax + f1;
            Umin(:, k) = umin + f1;
        elseif k*Ts > 20
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
        Y(:, k) = C*X(:, k);                                  % Discrete-time output

        %% Sensor fault income
        %TODO: Add output fails
        Yfail(:, k) = Y(:, k);
        
        %% RUIO 1
        RUIO(1).Phi(:, k+1) = RUIO(1).K*RUIO(1).Phi(:, k) + RUIO(1).L_ast*Yfail(:, k) + RUIO(1).B_bar_1*U(:, k) + RUIO(1).delta_bar_1;
        RUIO(1).X(:, k) = RUIO(1).T*[RUIO(1).Phi(:, k); RUIO(1).U_1*Yfail(:, k)-RUIO(1).U_1*RUIO(1).C_tilde_1*RUIO(1).Phi(:, k)];

        RUIO(1).Fact(k) = RUIO(1).U_1*(X(:, k+1) - RUIO(1).C_tilde_1*RUIO(1).Phi(:, k+1)) + RUIO(1).A_bar_22*RUIO(1).U_1*(RUIO(1).C_tilde_1*RUIO(1).Phi(:, k) - Yfail(:, k)) - RUIO(1).A_bar_21*RUIO(1).Phi(:, k) - RUIO(1).B_bar_2*U(:, k) - RUIO(1).delta_bar_2;

        % Error norm 1
        RUIO(1).error(k) = sqrt((RUIO(1).X(1, k)-Yfail(1, k))^2 + (RUIO(1).X(2, k)-Yfail(2, k))^2 + (RUIO(1).X(3, k)-Yfail(3, k))^2);

        %% RUIO 2
        RUIO(2).Phi(:, k+1) = RUIO(2).K*RUIO(2).Phi(:, k) + RUIO(2).L_ast*Yfail(:, k) + RUIO(2).B_bar_1*U(:, k) + RUIO(2).delta_bar_1;
        RUIO(2).X(:, k) = RUIO(2).T*[RUIO(2).Phi(:, k); RUIO(2).U_1*Yfail(:, k)-RUIO(2).U_1*RUIO(2).C_tilde_1*RUIO(2).Phi(:, k)];

        RUIO(2).Fact(k) = RUIO(2).U_1*(X(:, k+1) - RUIO(2).C_tilde_1*RUIO(2).Phi(:, k+1)) + RUIO(2).A_bar_22*RUIO(2).U_1*(RUIO(2).C_tilde_1*RUIO(2).Phi(:, k) - Yfail(:, k)) - RUIO(2).A_bar_21*RUIO(2).Phi(:, k) - RUIO(2).B_bar_2*U(:, k) - RUIO(2).delta_bar_2;

        % Error norm 2
        RUIO(2).error(k) = sqrt((RUIO(2).X(1, k)-Yfail(1, k))^2 + (RUIO(2).X(2, k)-Yfail(2, k))^2 + (RUIO(2).X(3, k)-Yfail(3, k))^2);

        % State limits
        for j = 1:nu
            if X(j, k) >= xmax(j)
                X(j, k) = xmax(j);
            elseif X(j, k) <= xmin(j)
                X(j, k) = xmin(j);
            end
        end
        
        % Next step fault compensation
        %TODO: Add exponential threshold
        if  k*Ts>= 2 % Estimation is ok
            if RUIO(1).error(k) < 0.1  % Error 1 threshold
                RUIO(1).FQ(k) = true;
                RUIO(2).Fact(k) = 0;
            else
                RUIO(1).FQ(k) = false;
            end
            if RUIO(2).error(k) < 0.16  % Error 2 threshold
                RUIO(2).FQ(k) = true;
                RUIO(1).Fact(k) = 0;
            else
                RUIO(2).FQ(k) = false;
            end
        else
            RUIO(1).FQ(k) = false; RUIO(2).FQ(k) = false;
            RUIO(1).Fact(k) = 0; RUIO(2).Fact(k) = 0;
        end
        
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

    % Inputs
    figure(2)
    if FTC == 0
        subplot(211)
        stairs(t, U(1, :), 'b', 'LineWidth', 1.5)
        hold on
        xlabel('Time [min]'); ylabel('q_1 [l/min]'); grid on
        subplot(212)
        stairs(t, U(2, :), 'b', 'LineWidth', 1.5)
        hold on
        xlabel('Time [min]'); ylabel('q_2 [l/min]'); grid on
	else
        subplot(211)
        stairs(t, U(1, :), 'k-.', 'LineWidth', 1.5)
        plot(t, umin(1)*ones(length(t)), 'r--')
        plot(t, umax(1)*ones(length(t)), 'r--')
        hold off
        legend('MPC', 'FTMPC', 'Location', 'NorthWest');
        legend boxoff
        subplot(212)
        stairs(t, U(2, :), 'k-.', 'LineWidth', 1.5)
        plot(t, umin(2)*ones(length(t)), 'r--')
        plot(t, umax(2)*ones(length(t)), 'r--')
        hold off
        print -dsvg figs/inputHE.svg
    end

    % Failure inputs
    figure(3)
    if FTC == 0
        subplot(211)
        plot(t, Umin(1, :), 'r--')
        hold on
        plot(t, Umax(1, :), 'r--')
        stairs(t, Ufail(1, :), 'b', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
        subplot(212)
        plot(t, Umin(2, :), 'r--')
        hold on
        plot(t, Umax(2, :), 'r--')
        stairs(t, Ufail(2, :), 'b', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
	else
        subplot(211)
        stairs(t, Ufail(1, :), 'k-.', 'LineWidth', 1.5)
        hold off
        subplot(212)
        stairs(t, Ufail(2, :), 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/inputfailHE.svg
    end

    % Estimation error
    figure(4)
    if FTC == 0
        subplot(211)
        plot(t, RUIO(1).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        xlabel('Time [min]'); ylabel('|e_x|');
        subplot(212)
        plot(t, RUIO(2).error, 'b', 'LineWidth', 1.5)
        hold on; grid on        
        xlabel('Time [min]'); ylabel('|e_x|');
	else
        subplot(211)
        plot(t, RUIO(1).error, 'k-.', 'LineWidth', 1.5)
        plot(t, 0.1*ones(length(t)),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 6.5])
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
        legend boxoff
        subplot(212)
        plot(t, RUIO(2).error, 'k-.', 'LineWidth', 1.5)
        plot(t, 0.16*ones(length(t)),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.6])
        hold off
        print -dsvg figs/errorHE.svg
    end

    % Fault estimation
    figure(5)
    if FTC == 0
        subplot(211)
        stairs(t, RUIO(1).Fact, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, Ufail(1, :) - U(1, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
        subplot(212)
        stairs(t, RUIO(2).Fact, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, Ufail(2, :) - U(2, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
    else
        subplot(211)
        stairs(t, RUIO(1).Fact, 'k-.', 'LineWidth', 1.5)
        hold off
        subplot(212)
        stairs(t, RUIO(2).Fact, 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/estimationHE.svg
    end

    % Objective
    figure(6)
    if FTC == 0
        plot(t, Obj, 'b', 'LineWidth', 1.5)
        hold on
        xlabel('Time [min]'); ylabel('Objective'); grid on
    else
        plot(t, Obj, 'k-.', 'LineWidth', 1.5)
        hold off
        axis([0 inf 0 400])
        print -dsvg figs/objectiveHE.svg
    end
    
    % State evolution
    figure(7)
    if FTC == 0
        plot3(x0(1), x0(2), x0(3), 'g*', 'LineWidth', 1.5);
        hold on
        plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
        plot3(Y(1, :), Y(2, :), Y(3, :), 'b', 'LineWidth', 1.5)
        plot3(Y(1, end), Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
        plot3(xsp(1), xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
        plot(Xpoly+X_lin, 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
        plot(Xs+X_lin, 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
        xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); zlabel('\theta_3 [K]'); grid on;
    else
        plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)        
        plot3(Y(1, :), Y(2, :), Y(3, :), 'k-.', 'LineWidth', 1.5)
        plot3(Y(1, end), Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/stateHE.svg
    end
    
    % State evolution
    Xx = Xpoly.projection(1:2).minHRep();
    Xxs = Xs.projection(1:2).minHRep();

    figure(8)
    if FTC == 0
        plot(x0(1), x0(2), 'g*', 'LineWidth', 1.5);
        hold on
        plot(Y_sim(1, :), Y_sim(2, :), 'y', 'LineWidth', 1.5)        
        plot(Y(1, :), Y(2, :), 'b.', 'LineWidth', 1.5)
        plot(Y(1, end), Y(2, end), 'mo', 'LineWidth', 1.5)
        plot(xsp(1), xsp(2), 'rp', 'LineWidth', 1.5)
        plot(Xxs+X_lin(1:2), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
        plot(Xx+X_lin(1:2), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
        xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); grid on;
    else
        plot(Y_sim(1, :), Y_sim(2, :), 'g', 'LineWidth', 1.5)        
        plot(Y(1, :), Y(2, :), 'k.', 'LineWidth', 1.5)
        plot(Y(1, end), Y(2, end), 'ro', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/state1-2HE.svg
    end
    
    % State evolution
    Xx = Xpoly.projection(1:2:3).minHRep();
    Xxs = Xs.projection(1:2:3).minHRep();

    figure(9)
    if FTC == 0
        plot(x0(1), x0(3), 'gd', 'LineWidth', 1.5);
        hold on
        plot(Y_sim(1, :), Y_sim(3, :), 'b', 'LineWidth', 1.5)
        plot(Y(1, end), Y(3, end), 'mo', 'LineWidth', 1.5)
        xlabel('\theta_1 [K]'); ylabel('\theta_3 [K]'); grid on;
    else
        plot(Y_sim(1, :), Y_sim(3, :), 'k--', 'LineWidth', 1.5)
        plot(Y(1, end), Y(3, end), 'y*', 'LineWidth', 1.5)
        plot(xsp(1), xsp(3), 'ro', 'LineWidth', 1.5)
        plot(Xxs+X_lin(1:2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'LineWidth', 1.5)
        plot(Xx+X_lin(1:2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'LineWidth', 1.5)
        hold off
        leg = legend('$x(0)$', '$x_{MPC}$', '$x_{MPC}(end)$', '$x_{FTMPC}$', '$x_{FTMPC}(end)$', '$x_s$', '$\bf{X}_s$', '$\bf{X}$', 'Location', 'SouthEast');
        set(leg, 'Interpreter', 'latex');
        print -dsvg figs/state1-3HE.svg
    end
    
    % State evolution
    Xx = Xpoly.projection(2:3).minHRep();
    Xxs = Xs.projection(2:3).minHRep();

    figure(10)
    if FTC == 0
        plot(x0(2), x0(3), 'g*', 'LineWidth', 1.5);
        hold on
        plot(Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
        plot(Y(2, :), Y(3, :), 'b.', 'LineWidth', 1.5)
        plot(Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
        plot(xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
        plot(Xxs+X_lin(2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
        plot(Xx+X_lin(2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
        xlabel('\theta_2 [K]'); ylabel('\theta_3 [K]'); grid on;
    else
        plot(Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)
        plot(Y(2, :), Y(3, :), 'k.', 'LineWidth', 1.5)
        plot(Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/state2-3HE.svg    
    end
end
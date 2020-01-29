%% MPC on CSTR
clc; clear; yalmip('clear');
close all;

RUIO = struct;

%% LTI model
Ts = 0.05;               % Sample time [min]
Theta_1s = 498;     % Temperatura de salida de fluido 1 (K) (set-point) (460)
Theta_2s = 690;     % Temperatura de salida de fluido 2 (K)

run HE;

%% Parameters
Np = 5;                                   % Prediction horizon
Time = 25;                              % Simulation end time 
Nsim = Time/Ts;                     % Simulation steps
t = 0:Ts:Time;                         % Simulation time
x0 = [495; 680; 570];          % Start-point
f1 = [0; 1]; f2 = [-8; 0];      % Fault magnitude

% Set-point
xsp = [Theta_1s; Theta_2s; Theta_p]; % Set-point

% Constraints
xmin = [495; 650; 530];
xmax = [500; 750; 590];
umin = [90; 6];
umax = [110; 10];

% Wheight matrix
Qx = eye(nx);
% Qx = diag([1 1 1]);
%Rx = eye(nu);
Rx = diag([1 0.1]);
gamma = 1e2*diag([1e6 1 1]); 

% %% Max reachable set
% X = Polyhedron('lb', xmin-X_lin, 'ub', xmax-X_lin);     % State polyhedron
% U = Polyhedron('lb', umin-U_lin, 'ub', umax-U_lin);    % Input polyhedron
% Acl = [Ad, Bd; zeros(nu, nx) eye(nu)];
% Maxiter = 50;
% Xs = X;
% for k = 1:Maxiter
%     Xo = Xs;
%     Z = Xo*U;   % Extended Set
%     
%     % Next set
%     S = Z.invAffineMap(Acl);
%     S = S.intersect(Z);
%     Xs = S.projection(1:nx).minHRep();
%     Xs = Xs.intersect(Xo).minHRep();
%     
%     % Check Invariant
% 	if Xs == Xo
%         break;
% 	end
% end

%% Reduced-order unknown input observer
N = 2;
run HE_RUIO;

%% MPC controller
run MPC;

%% Simulation Setup
U = zeros(nu, Nsim);                   % Control Input
Ufail = zeros(nu, Nsim);               % Fault control Input
X = zeros(nx, Nsim+1);               % States
Y = zeros(ny, Nsim);                    % Measure outputs
Yfail = zeros(ny, Nsim);                % Faulty measure outputs

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
U(:, 1) = U_lin;
RUIO(1).X(:, 1) = x0;
RUIO(2).X(:, 1) = x0;
UIOO(1).X(:, 1) = x0;
UIOO(2).X(:, 1) = x0;
              
%% Simulation
for FTC = 0:0 % 0 - FTC is off; 1 - FTC is on
    
    % Vector initialization for plots
    obj = NaN; state = x0; input = [NaN; NaN]; u_max = umax; u_min = umin;
    error1 = NaN; x1_hat = x0; f1est = 0; error2 = NaN; x2_hat = x0; f2est = 0;
    ufail = [NaN; NaN]; ufail_max = umax; ufail_min = umin; y = x0;

    % Vector initialization for lsim
    x = x0; Y = x0;
    umin0 = umin; umax0 = umax;
    Umin = umin; Umax = umax;
    uf = [0; 0];
    for k = 1:Nsim

        % Run MPC controller       
        [sol, diag] = mpc{x, xsp, uf};
        if diag
            msg = ['Infeasible problem at t = ', num2str(k*Ts)];
            disp(msg)
            return;
        end
        umpc = sol{1}; Obj = sol{2};

        % Input limits
        for j = 1:nu
            if umpc(j) >= umax(j)
                umpc(j) = umax(j);
            elseif umpc(j) <= umin(j)
                umpc(j) = umin(j);
            end
        end

        % Add fails
        if k*Ts >= 5 && k*Ts < 15
            Ufail = umpc + f1;
            Umax = umax0 + f1;
            Umin = umin0 + f1;
        elseif k*Ts >= 20
            Ufail = umpc + f2;
            Umax = umax0 + f2;
            Umin = umin0 + f2;
        else
            Ufail = umpc;
            Umax = umax0;
            Umin = umin0;
        end

        % Natural system saturation
        for j = 1:nu
            if Ufail(j) >= Umax(j)
                Ufail(j) = Umax(j);
            elseif Ufail(j) <= Umin(j)
                Ufail(j) = Umin(j);
            end
        end
        
        U(:, k) = umpc;
        
        % Continuous-time simulation (reality)
        Dt  = linspace(0, Ts, 10)';
        u = (Ufail - U_lin)*ones(1, numel(Dt));% Zero-order hold input (with Input offset)
        X_sim = lsim(sys, u, Dt, X(:, k)-X_lin);  % Nominal system (with state offset)
        X_sim = X_sim' + X_lin;
        X(:, k+1) = X_sim(:, end);
        Y(:, k) = C*X(:, k);

        %% Sensor fault income
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

        x = X_sim(:, end);                                              % System output with offset in discrete time
        % state limits
        for j = 1:nu
            if x(j) >= xmax(j)
                x(j) = xmax(j);
            elseif x(j) <= xmin(j)
                x(j) = xmin(j);
            end
        end
        
%         % Next Step Fault compensation
%         if  i*Ts>= 2 % Estimation is ok
%             if Error1 < 0.1  % Error 1 Threshold
%                 F2est = 0;
%             end
%             if Error2 < 0.1  % Error 2 Threshold
%                 F1est = 0;
%             end
%         else
%             F1est = 0; F2est = 0;
%         end
%         if FTC == 1
%             uf = [F1est; F2est];
%         else
%             uf = [0; 0];
%         end

        y = [y X_sim];                                    % System output with offset
        ufail = [ufail Ufail];                             % Failure input sequence    
        input = [input umpc];                        % Predicted input sequence    
        state = [state x];                               % System state
        obj = [obj Obj];                                  % Objective function
        u_max = [u_max umax];
        u_min = [u_min umin];
        ufail_max = [ufail_max Umax];
        ufail_min = [ufail_min Umin];

%         x1_hat = [x1_hat X1_hat];
%         error1 = [error1 Error1];
%         f1est = [f1est F1est];
% 
%         x2_hat = [x2_hat X2_hat];
%         error2 = [error2 Error2];
%         f2est = [f2est F2est];
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
        plot(tc, y(1, :), 'b', 'LineWidth', 1.5);
        % plot(t, xmin(1)*ones(length(t)), 'y--')
        % plot(t, xmax(1)*ones(length(t)), 'y--')
        % plot(t, state(1, :), 'b-.', t, x1_hat(1, :), 'g:', t, x2_hat(1, :), 'y:', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
        xlabel('Time [min]'); ylabel('\theta_1 [K]'); grid on
        subplot(312)
%         plot(t, state(2, :), 'b.', tc, y(2, :), 'y-.', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(2)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
        h(2) = plot(tc, y(2, :), 'b', 'LineWidth', 1.5);
        % plot(t, xmin(2)*ones(length(t)), 'y--')
        % plot(t, xmax(2)*ones(length(t)), 'y--')
        % plot(t, state(2, :), 'b-.', t, x1_hat(2, :), 'g:', t, x2_hat(2, :), 'y:', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_2 [K]'); grid on
        subplot(313)
%         plot(t, state(3, :), 'b.', tc, y(3, :), 'y-.', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(3)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
        plot(tc, y(3, :), 'b', 'LineWidth', 1.5);   
        % plot(t, xmin(3)*ones(length(t)), 'y--')
        % plot(t, xmax(3)*ones(length(t)), 'y--')
        % plot(t, state(3, :), 'b-.', t, x1_hat(3, :), 'g:', t, x2_hat(3, :), 'y:', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_p [K]'); grid on
    else
        subplot(311)
%         plot(t, state(1, :), 'k.', tc, y(1, :), 'g-.', 'LineWidth', 1.5);
        plot(tc, y(1, :), 'k-.', 'LineWidth', 1.5);
        hold off
        legend('x_s', 'Location', 'SouthEast');
        legend boxoff
        subplot(312)
%         plot(t, state(2, :), 'k.', tc, y(2, :), 'g-.', 'LineWidth', 1.5);
        plot(tc, y(2, :), 'k-.', 'LineWidth', 1.5);        
        hold off
        legend(h(2), 'MPC', 'Location', 'SouthEast');
        legend boxoff
        subplot(313)
%         plot(t, state(3, :), 'k.', tc, y(3, :), 'g-.', 'LineWidth', 1.5);
        h(3) = plot(tc, y(3, :), 'k-.', 'LineWidth', 1.5);
        hold off
        legend(h(3), 'FTMPC', 'Location', 'NorthEast');
        legend boxoff
        pause(2)
        print -dsvg figs/outputHE.svg
    end
    % print -dsvg outputCSTR.svg

%     % Inputs
%     figure(2)
%     if FTC == 0
%         subplot(211)
%         stairs(t, input(1, :), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_1 [l/min]'); grid on
%         subplot(212)
%         stairs(t, input(2, :), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_2 [l/min]'); grid on
% 	else
%         subplot(211)
%         stairs(t, input(1, :), 'k-.', 'LineWidth', 1.5)
%         plot(t, u_min(1, :), 'r--')
%         plot(t, u_max(1, :), 'r--')
%         hold off
%         legend('MPC', 'FTMPC', 'Location', 'NorthWest');
%         legend boxoff
%         subplot(212)
%         stairs(t, input(2, :), 'k-.', 'LineWidth', 1.5)
%         plot(t, u_min(2, :), 'r--')
%         plot(t, u_max(2, :), 'r--')
%         hold off
%         pause(2)
%         print -dsvg inputHE.svg
%     end
    % print -dsvg inputCSTR.svg

%     % Failure inputs
%     figure(3)
%     if FTC == 0
%         subplot(211)
%         plot(t, ufail_min(1, :), 'r--')
%         hold on
%         plot(t, ufail_max(1, :), 'r--')
%         stairs(t, ufail(1, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
%         subplot(212)
%         plot(t, ufail_min(2, :), 'r--')
%         hold on
%         plot(t, ufail_max(2, :), 'r--')
%         stairs(t, ufail(2, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
% 	else
%         subplot(211)
%         stairs(t, ufail(1, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, ufail(2, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg inputfailCSTR.svg    
% 
%     % Estimation error
%     figure(4)
%     if FTC == 0
%         subplot(211)
%         plot(t, error1, 'b', 'LineWidth', 1.5)
%         hold on; grid on
%         xlabel('Time [min]'); ylabel('|e_x|');
%         subplot(212)
%         plot(t, error2, 'b', 'LineWidth', 1.5)
%         hold on; grid on        
%         xlabel('Time [min]'); ylabel('|e_x|');
% 	else
%         subplot(211)
%         plot(t, error1, 'k-.', 'LineWidth', 1.5)
%         plot(t, 0.1*ones(length(t)),  'r--', 'LineWidth', 1.5)
%         hold off
%         legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
%         legend boxoff
%         subplot(212)
%         plot(t, error2, 'k-.', 'LineWidth', 1.5)
%         plot(t, 0.1*ones(length(t)),  'r--', 'LineWidth', 1.5)
%         hold off
%         pause(2)
%         print -dsvg errorHE.svg
%     end
    % print -dsvg errorCSTR.svg

%     % Fault estimation
%     figure(5)
%     if FTC == 0
%         subplot(211)
%         stairs(t, f1est, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, ufail(1, :) - input(1, :), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
%         subplot(212)
%         stairs(t, f2est, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, ufail(2, :) - input(2, :), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
%     else
%         subplot(211)
%         stairs(t, f1est, 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, f2est, 'k-.', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg figs/estimationCSTR.svg


%     % Objective
%     figure(6)
%     if FTC == 0
%         plot(t, obj, 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('Objective'); grid on
%     else
%         plot(t, obj, 'k-.', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg objectiveCSTR.svg
    
%     % State evolution
%     figure(7)
%     if FTC == 0
%         plot3(x0(1), x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot3(y(1, :), y(2, :), y(3, :), 'y', 'LineWidth', 1.5)
%         plot3(state(1, :), state(2, :), state(3, :), 'b', 'LineWidth', 1.5)
%         plot3(state(1, end), state(2, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         plot3(xsp(1), xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(X+X_lin, 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xs+X_lin, 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); zlabel('\theta_3 [K]'); grid on;
%     else
%         plot3(y(1, :), y(2, :), y(3, :), 'g', 'LineWidth', 1.5)        
%         plot3(state(1, :), state(2, :), state(3, :), 'k-.', 'LineWidth', 1.5)
%         plot3(state(1, end), state(2, end), state(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg stateCSTR.svg
    
%     % State evolution
%     Xx = X.projection(1:2).minHRep();
%     Xxs = Xs.projection(1:2).minHRep();
% 
%     figure(8)
%     if FTC == 0
%         plot(x0(1), x0(2), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot(y(1, :), y(2, :), 'y', 'LineWidth', 1.5)        
%         plot(state(1, :), state(2, :), 'b.', 'LineWidth', 1.5)
%         plot(state(1, end), state(2, end), 'mo', 'LineWidth', 1.5)
%         plot(xsp(1), xsp(2), 'rp', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(1:2), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xx+X_lin(1:2), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_2 [K]'); grid on;
%     else
%         plot(y(1, :), y(2, :), 'g', 'LineWidth', 1.5)        
%         plot(state(1, :), state(2, :), 'k.', 'LineWidth', 1.5)
%         plot(state(1, end), state(2, end), 'ro', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg stateCSTR.svg
%     
%     % State evolution
%     Xx = X.projection(1:2:3).minHRep();
%     Xxs = Xs.projection(1:2:3).minHRep();
% 
%     figure(9)
%     if FTC == 0
%         plot(x0(1), x0(3), 'gd', 'LineWidth', 1.5);
%         hold on
%         plot(y(1, :), y(3, :), 'b', 'LineWidth', 1.5)
% %         plot(state(1, :), state(3, :), 'b.', 'LineWidth', 1.5)
%         plot(state(1, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         xlabel('\theta_1 [K]'); ylabel('\theta_3 [K]'); grid on;
%     else
%         plot(y(1, :), y(3, :), 'k--', 'LineWidth', 1.5)
% %         plot(state(1, :), state(3, :), 'k.', 'LineWidth', 1.5)
%         plot(state(1, end), state(3, end), 'y*', 'LineWidth', 1.5)
%         plot(xsp(1), xsp(3), 'ro', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(1:2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'LineWidth', 1.5)
%         plot(Xx+X_lin(1:2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'LineWidth', 1.5)
%         hold off
%         leg = legend('$x(0)$', '$x_{MPC}$', '$x_{MPC}(end)$', '$x_{FTMPC}$', '$x_{FTMPC}(end)$', '$x_s$', '$\bf{X}_s$', '$\bf{X}$', 'Location', 'SouthEast');
%         set(leg, 'Interpreter', 'latex');
%         pause(2)
%         print -dsvg stateHE.svg
%     end
%     % print -dsvg stateCSTR.svg
%     
%     % State evolution
%     Xx = X.projection(2:3).minHRep();
%     Xxs = Xs.projection(2:3).minHRep();
% 
%     figure(10)
%     if FTC == 0
%         plot(x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot(y(2, :), y(3, :), 'y', 'LineWidth', 1.5)
%         plot(state(2, :), state(3, :), 'b.', 'LineWidth', 1.5)
%         plot(state(2, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         plot(xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xx+X_lin(2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('\theta_2 [K]'); ylabel('\theta_3 [K]'); grid on;
%     else
%         plot(y(2, :), y(3, :), 'g', 'LineWidth', 1.5)
%         plot(state(2, :), state(3, :), 'k.', 'LineWidth', 1.5)
%         plot(state(2, end), state(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg stateCSTR.svg    
end
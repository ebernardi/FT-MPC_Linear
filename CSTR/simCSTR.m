%% MPC on CSTR
clc; clear; yalmip('clear');
close all;

%% LTI model
Ts = 0.05;                       % Sample time [min]
Vr = 100;                        % Reactor volume [l]
Tr = 440;                         % Temperature [K]

run CSTR;

% Continuous-time system
sys = ss(A, B, C, D);

% Euler discretization method
Ad = (A*Ts) + eye(nx); Bd = B*Ts; Cd = C; Dd = D; deltad = delta*Ts;

%% Parameters
N = 5;                                     % Prediction horizon
Time = 20;                              % Simulation end time 
Nsim = Time/Ts;                     % Simulation steps
t = 0:Ts:Time;                         % Simulation time
x0 = [95; 0.105; 437];            % Start-point
f1 = [0; -8]; f2 = [8; 0];          % Fault magnitude

% Set-point
% Vr = 103;                           % Reactor volume [l]
% Tr = 450;                            % Temperature [K]
% Ca = CAe/(1+(k0*(Vr/q)*exp(-E_R/Tr)));% Concentration [mol/l]
xsp = [Vr; Ca; Tr];  % Set-point

% Constraints
xmin = [90; 0.08; 434];
xmax = [105; 0.11; 444];
umin = [90; 90];
umax = [110; 110];

% Wheight matrix
Qx = eye(nx);
% Qx = diag([1 1 1]);
Rx = eye(nu);
% Rx = diag([1 1]);
gamma = 1e2*diag([1 1e4 1e2]); 

%% Max reachable set
X = Polyhedron('lb', xmin-X_lin, 'ub', xmax-X_lin);     % State polyhedron
U = Polyhedron('lb', umin-U_lin, 'ub', umax-U_lin);    % Input polyhedron
Acl = [Ad, Bd; zeros(nu, nx) eye(nu)];
Maxiter = 50;
Xs = X;
for i = 1:Maxiter
    Xo = Xs;
    Z = Xo*U;   % Extended Set
    
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
run RUIO;
nobs = 2;

%% MPC controller
run MPC;

Aso = [Ad zeros(nx, nobs+nobs); L_ast_1 K_1 zeros(nobs); L_ast_2 zeros(nobs) K_2];
Bso = [Bd zeros(nx, nu); zeros(nobs, nu) B1_bar_1; zeros(nobs, nu) B2_bar_1];
Cso = [T_1*[zeros(nobs, nx) eye(nobs) zeros(nobs); U1_1 -U1_1*Cd*N_1 zeros(1, nobs)]; ...
            T_2*[zeros(nobs, nx+nobs) eye(nobs); U1_2 zeros(1, nobs) -U1_2*Cd*N_2]];
Dso = zeros(nx+nx, nu+nu);
sys_obs = ss(Aso, Bso, Cso, Dso, Ts);
            
%% Simulation
for FTC = 0:1 % 0 - FTC is off; 1 - FTC is on
    
    % Vector initialization for plots
    obj = NaN; state = x0; input = [NaN; NaN]; u_max = umax; u_min = umin;
    error1 = NaN; x1_hat = x0; f1est = 0; error2 = NaN; x2_hat = x0; f2est = 0;
    ufail = [NaN; NaN]; ufail_max = umax; ufail_min = umin; y = x0;

    % Vector initialization for lsim
    theta1 = [0; 0]; theta2 = [0; 0];
    x = x0; Y = x0;
    umin0 = umin; umax0 = umax;
    Umin = umin; Umax = umax;
    F1est = 0; F2est = 0; Error1 = 0; Error2 = 0;
    uf = [0; 0];
    xk = x0-X_lin;
    uk = umin-U_lin;
    theta1k = theta1; theta2k = theta2;
    
    for i = 1:Nsim
        
        % Fault of first observer
        F1est = U1_1*((x-X_lin) - Cd*N_1*theta1) + A1_bar_22*U1_1*(Cd*N_1*theta1k - xk) - (A1_bar_21*theta1k) - (B1_bar_2*uk);
        % Fault of second observer
        F2est = U1_2*((x-X_lin) - Cd*N_2*theta2) + A2_bar_22*U1_2*(Cd*N_2*theta2k - xk) - (A2_bar_21*theta2k) - (B2_bar_2*uk);
        
        % Next Step Fault compensation
        if  i*Ts>= 2 % Estimation is ok
            if Error1 < 0.1  % Error 1 Threshold
                F2est = 0;
            end
            if Error2 < 0.1  % Error 2 Threshold
                F1est = 0;
            end
        else
            F1est = 0; F2est = 0;
        end
        if FTC == 1
            uf = [F1est; F2est];
        else
            uf = [0; 0];
        end
        
        % Run MPC controller
        [sol, diag] = mpc{x, xsp, uf};
        if diag
            msg = ['Infeasible problem at t = ', num2str(i*Ts)];
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
        if i*Ts >= 5 && i*Ts < 10
            Ufail = umpc + f2;
            Umax = umax0 + f2;
            Umin = umin0 + f2;
        elseif i*Ts >= 12
            Ufail = umpc + f1;
            Umax = umax0 + f1;
            Umin = umin0 + f1;
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

        % Past values for Fault-estimation
        xk = x-X_lin;
        uk = umpc-U_lin;
        theta1k = theta1;
        theta2k = theta2;
        
        % Continuous-time simulation (reality)
        Dt  = linspace(0, Ts, 10)';
        u = (Ufail - U_lin)*ones(1, numel(Dt));% Zero-order hold input (with Input offset)
        Y = lsim(sys, u, Dt, x-X_lin);  % Nominal system (with state offset)
        Y = Y'+X_lin;
        y = [y Y];                                    % System output with offset
        
        % Discrete-time Observer simulation
        Dt  =0:Ts:Ts;                        % Vector of instants to be passed to lsim
        uzoh = [(Ufail - U_lin); (umpc - U_lin)]*ones(1, numel(Dt));  % Zero-order hold input (with Input offset)
        [yout, tout, xout] = lsim(sys_obs, uzoh, Dt, [x-X_lin; theta1; theta2]);  % Nominal system (with state offset)
%         x = xout(end, 1:3)' + X_lin;                    % System output with offset in discrete time
        x = Y(:, end);                                              % System output with offset in continuos time
        theta1 = xout(end, 4:5)';                           % States of first observer
        theta2 = xout(end, 6:7)';                           % States of second observer    
        X1_hat = yout(end, 1:3)' + X_lin;               % Estimation of first observer
        X2_hat = yout(end, 4:6)' + X_lin;               % Estimation of second observer    
        Error1 = norm(X1_hat-x, 2);                      % Error 1 detection 
        Error2 = norm(X2_hat-x, 2);                      % Error 2 detection 

        % State limits
        for j = 1:nu
            if x(j) >= xmax(j)
                x(j) = xmax(j);
            elseif x(j) <= xmin(j)
                x(j) = xmin(j);
            end
        end
        
        ufail = [ufail Ufail];                             % Failure input sequence    
        input = [input umpc];                        % Predicted input sequence    
        state = [state x];                               % System state
        obj = [obj Obj];                                  % Objective function
        u_max = [u_max umax];
        u_min = [u_min umin];
        ufail_max = [ufail_max Umax];
        ufail_min = [ufail_min Umin];

        x1_hat = [x1_hat X1_hat];
        error1 = [error1 Error1];
        f1est = [f1est F1est];

        x2_hat = [x2_hat X2_hat];
        error2 = [error2 Error2];
        f2est = [f2est F2est];
    end

    %% Set Plots
    vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];
    tc = 0:Ts/10:Nsim/20;
    
%     % Outputs
%     figure(1)
%     if FTC == 0
%         subplot(311)
% %         plot(t, state(1, :), 'b.', tc, y(1, :), 'y-.', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
%         plot(t, xsp(1)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
%         hold on
%         plot(tc, y(1, :), 'b', 'LineWidth', 1.5);
%         % plot(t, xmin(1)*ones(length(t)), 'y--')
%         % plot(t, xmax(1)*ones(length(t)), 'y--')
%         % plot(t, state(1, :), 'b-.', t, x1_hat(1, :), 'g:', t, x2_hat(1, :), 'y:', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
%         xlabel('Time [min]'); ylabel('Vol [l]'); grid on; axis([-inf inf 95 102]); yticks([96 98 100 102]);
%         subplot(312)
% %         plot(t, state(2, :), 'b.', tc, y(2, :), 'y-.', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
%         plot(t, xsp(2)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
%         hold on
%         h(2) = plot(tc, y(2, :), 'b', 'LineWidth', 1.5);
%         % plot(t, xmin(2)*ones(length(t)), 'y--')
%         % plot(t, xmax(2)*ones(length(t)), 'y--')
%         % plot(t, state(2, :), 'b-.', t, x1_hat(2, :), 'g:', t, x2_hat(2, :), 'y:', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('C_A [mol/l]'); grid on; axis([-inf inf 0.08 0.11]); yticks([0.08 0.09 0.1 0.11]);
%         subplot(313)
% %         plot(t, state(3, :), 'b.', tc, y(3, :), 'y-.', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
%         plot(t, xsp(3)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
%         hold on
%         plot(tc, y(3, :), 'b', 'LineWidth', 1.5);        
%         % plot(t, xmin(3)*ones(length(t)), 'y--')
%         % plot(t, xmax(3)*ones(length(t)), 'y--')
%         % plot(t, state(3, :), 'b-.', t, x1_hat(3, :), 'g:', t, x2_hat(3, :), 'y:', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Temp [K]'); grid on; axis([-inf inf 437 442]); yticks([438 440 442]);
%     else
%         subplot(311)
% %         plot(t, state(1, :), 'k.', tc, y(1, :), 'g-.', 'LineWidth', 1.5);
%         plot(tc, y(1, :), 'k-.', 'LineWidth', 1.5);
%         hold off
%         legend('x_s', 'Location', 'SouthEast');
%         legend boxoff
%         subplot(312)
% %         plot(t, state(2, :), 'k.', tc, y(2, :), 'g-.', 'LineWidth', 1.5);
%         plot(tc, y(2, :), 'k-.', 'LineWidth', 1.5);        
%         hold off
%         legend(h(2), 'MPC', 'Location', 'NorthEast');
%         legend boxoff
%         subplot(313)
% %         plot(t, state(3, :), 'k.', tc, y(3, :), 'g-.', 'LineWidth', 1.5);
%         h(3) = plot(tc, y(3, :), 'k-.', 'LineWidth', 1.5);
%         hold off
%         legend(h(3), 'FTMPC', 'Location', 'SouthEast');
%         legend boxoff
%         pause(2)
%         print -dsvg outputCSTR.svg
%     end

%     % Inputs
%     figure(2)
%     if FTC == 0
%         subplot(211)
%         stairs(t, input(1, :), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_s [l/min]'); grid on
%         subplot(212)
%         stairs(t, input(2, :), 'b', 'LineWidth', 1.5)
%         hold on
%         xlabel('Time [min]'); ylabel('q_c [l/min]'); grid on
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
%         print -dsvg inputCSTR.svg
%     end

%     % Failure inputs
%     figure(3)
%     if FTC == 0
%         subplot(211)
%         plot(t, ufail_min(1, :), 'r--')
%         hold on
%         plot(t, ufail_max(1, :), 'r--')
%         stairs(t, ufail(1, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('q [l/min]'); grid on
%         subplot(212)
%         plot(t, ufail_min(2, :), 'r--')
%         hold on
%         plot(t, ufail_max(2, :), 'r--')
%         stairs(t, ufail(2, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('qc [l/min]'); grid on
% 	else
%         subplot(211)
%         stairs(t, ufail(1, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, ufail(2, :), 'k-.', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg inputfailCSTR.svg    

    % Estimation error
    figure(4)
    if FTC == 0
        subplot(211)
        plot(t, error1, 'b', 'LineWidth', 1.5)
        hold on; grid on
        xlabel('Time [min]'); ylabel('r_{q_c}'); axis([-inf inf 0 0.8]);
        subplot(212)
        plot(t, error2, 'b', 'LineWidth', 1.5)
        hold on; grid on        
        xlabel('Time [min]'); ylabel('r_{q_s}'); axis([-inf inf 0 0.8]);
	else
        subplot(211)
        plot(t, error1, 'k-.', 'LineWidth', 1.5)
        plot(t, 0.2*ones(length(t)),  'r--', 'LineWidth', 1.5)
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');
        legend boxoff
        subplot(212)
        plot(t, error2, 'k-.', 'LineWidth', 1.5)
        plot(t, 0.2*ones(length(t)),  'r--', 'LineWidth', 1.5)
        hold off
        pause(2)
        print -dsvg errorCSTR.svg
    end

%     % Fault estimation
%     figure(5)
%     if FTC == 0
%         subplot(211)
%         stairs(t, f1est, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, ufail(1, :) - input(1, :), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('qc [l/min]'); grid on
%         subplot(212)
%         stairs(t, f2est, 'b', 'LineWidth', 1.5)
%         hold on
%         stairs(t, ufail(2, :) - input(2, :), 'm--', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('qs [l/min]'); grid on
%     else
%         subplot(211)
%         stairs(t, f1est, 'k-.', 'LineWidth', 1.5)
%         hold off
%         subplot(212)
%         stairs(t, f2est, 'k-.', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg estimationCSTR.svg
% 
% 
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
%     
%     % State evolution
%     figure(7)
%     if FTC == 0
%         plot3(x0(1), x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot3(y(1, :), y(2, :), y(3, :), 'y', 'LineWidth', 1.5)
%         plot3(state(1, :), state(2, :), state(3, :), 'b.', 'LineWidth', 1.5)
%         plot3(state(1, end), state(2, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         plot3(xsp(1), xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(Xs+X_lin, 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(X+X_lin, 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('Vol [l]'); ylabel('C_A [mol/l]'); zlabel('Temp [K]'); grid on;
%     else
%         plot3(y(1, :), y(2, :), y(3, :), 'g', 'LineWidth', 1.5)
%         plot3(state(1, :), state(2, :), state(3, :), 'k.', 'LineWidth', 1.5)
%         plot3(state(1, end), state(2, end), state(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg stateCSTR.svg
%     
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
%         xlabel('Vol [l]'); ylabel('C_A [mol/l]'); grid on;
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
%         plot(x0(1), x0(3), 'g*', 'LineWidth', 1.5);
%         hold on
%         plot(y(1, :), y(3, :), 'y', 'LineWidth', 1.5)
%         plot(state(1, :), state(3, :), 'b.', 'LineWidth', 1.5)
%         plot(state(1, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         plot(xsp(1), xsp(3), 'rp', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(1:2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%         plot(Xx+X_lin(1:2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%         xlabel('Vol [l]'); ylabel('Temp [K]'); grid on;
%     else
%         plot(y(1, :), y(3, :), 'g', 'LineWidth', 1.5)
%         plot(state(1, :), state(3, :), 'k.', 'LineWidth', 1.5)
%         plot(state(1, end), state(3, end), 'ro', 'LineWidth', 1.5)
%         hold off
%     end
%     % print -dsvg stateCSTR.svg
    
%     % State evolution
%     Xx = X.projection(2:3).minHRep();
%     Xxs = Xs.projection(2:3).minHRep();
% 
%     figure(10)
%     if FTC == 0
%         plot(x0(2), x0(3), 'gd', 'LineWidth', 1.5);
%         hold on
%         plot(y(2, :), y(3, :), 'b', 'LineWidth', 1.5)
% %         plot(state(2, :), state(3, :), 'b.', 'LineWidth', 1.5)
%         plot(state(2, end), state(3, end), 'mo', 'LineWidth', 1.5)
%         xlabel('C_A [mol/l]'); ylabel('Temp [K]'); grid on;
%     else
%         plot(y(2, :), y(3, :), 'k--', 'LineWidth', 1.5)
% %         plot(state(2, :), state(3, :), 'k.', 'LineWidth', 1.5)
%         plot(state(2, end), state(3, end), 'y*', 'LineWidth', 1.5)
%         plot(xsp(2), xsp(3), 'ro', 'LineWidth', 1.5)
%         plot(Xxs+X_lin(2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'LineWidth', 1.5)
%         plot(Xx+X_lin(2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'LineWidth', 1.5)
%         hold off
%         leg = legend('$x(0)$', '$x_{MPC}$', '$x_{MPC}(end)$', '$x_{FTMPC}$', '$x_{FTMPC}(end)$', '$x_s$', '$\bf{X}_s$', '$\bf{X}$', 'Location', 'SouthWest');
%         set(leg, 'Interpreter', 'latex');
%         pause(2)
%         print -dsvg stateCSTR.svg
%     end  
end
%% Set Plots
clear all; close all; clc;

load FTCS

vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];
tc = 0:Ts/10:Nsim/20;

for FT = 2:-1:1
    % Outputs
    figure(1)
    if FT == 2
        subplot(311)
%         plot(t, state(1, :), 'b.', tc, y(1, :), 'y-.', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
        plot(t, xsp(1)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
%         plot(tc, Y_sim(1, :), 'b', 'LineWidth', 1.5);
%         plot(t, Y(1, :), 'g:', t, Yfail(1, :), 'y-.', tc, Y_sim(1, :), 'b', 'LineWidth', 1.5);
        plot(t, FTCS(FT).Y(1, :), 'b', t, FTCS(FT).Yfail(1, :), 'y-.', 'LineWidth', 1.5);
        % plot(t, xmin(1)*ones(length(t)), 'y--')
        % plot(t, xmax(1)*ones(length(t)), 'y--')
        % plot(t, state(1, :), 'b-.', t, x1_hat(1, :), 'g:', t, x2_hat(1, :), 'y:', t, xsp(1)*ones(length(t)), 'r--', 'LineWidth', 1.5);
        xlabel('Time [min]'); ylabel('\theta_1 [K]'); grid on
        subplot(312)
%         plot(t, state(2, :), 'b.', tc, y(2, :), 'y-.', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(2)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
%         h(2) = plot(tc, Y_sim(2, :), 'b', 'LineWidth', 1.5);
        plot(t, FTCS(FT).Y(2, :), 'b', t, FTCS(FT).Yfail(2, :), 'y-.', 'LineWidth', 1.5);
        legend('Actual-MPC', 'Measured-MPC', 'Location', 'NorthEast');
%         plot(t, Y(2, :), 'g:', t, Yfail(2, :), 'y-.', tc, Y_sim(2, :), 'b', 'LineWidth', 1.5);
        % plot(t, xmin(2)*ones(length(t)), 'y--')
        % plot(t, xmax(2)*ones(length(t)), 'y--')
        % plot(t, state(2, :), 'b-.', t, x1_hat(2, :), 'g:', t, x2_hat(2, :), 'y:', t, xsp(2)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_2 [K]'); grid on
        subplot(313)
%         plot(t, state(3, :), 'b.', tc, y(3, :), 'y-.', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        plot(t, xsp(3)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
        hold on
%         plot(tc, Y_sim(3, :), 'b', 'LineWidth', 1.5);   
        plot(t, FTCS(FT).Y(3, :), 'b', t, FTCS(FT).Yfail(3, :), 'y-.', 'LineWidth', 1.5);
        % plot(t, xmin(3)*ones(length(t)), 'y--')
        % plot(t, xmax(3)*ones(length(t)), 'y--')
        % plot(t, state(3, :), 'b-.', t, x1_hat(3, :), 'g:', t, x2_hat(3, :), 'y:', t, xsp(3)*ones(length(t)), 'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_p [K]'); grid on
    else
        subplot(311)
%         plot(t, state(1, :), 'k.', tc, y(1, :), 'g-.', 'LineWidth', 1.5);
        plot(t, FTCS(FT).Y(1, :), 'k-.', t, FTCS(FT).Yfail(1, :), 'g:', 'LineWidth', 1.5);
%         plot(t, Y(1, :), 'g:', t, Yfail(1, :), 'y-.', tc, Y_sim(1, :), 'b', 'LineWidth', 1.5);
        hold off
        legend('x_s', 'Location', 'SouthEast');
        legend boxoff
        subplot(312)
%         plot(t, state(2, :), 'k.', tc, y(2, :), 'g-.', 'LineWidth', 1.5);
        plot(t, FTCS(FT).Y(2, :), 'k-.', t, FTCS(FT).Yfail(2, :), 'g:', 'LineWidth', 1.5);
%         plot(t, Y(2, :), 'g:', t, Yfail(2, :), 'y-.', tc, Y_sim(2, :), 'b', 'LineWidth', 1.5);
        legend('', 'Actual-MPC', 'Measured-MPC', 'Location', 'SouthEast');
        hold off
%         legend(h(2), 'MPC', 'Location', 'SouthEast');
        legend boxoff
        subplot(313)
%         plot(t, state(3, :), 'k.', tc, y(3, :), 'g-.', 'LineWidth', 1.5);
        plot(t, FTCS(FT).Y(3, :), 'k-.', t, FTCS(FT).Yfail(3, :), 'g:', 'LineWidth', 1.5);
%         h(3) = plot(t, Y(3, :), 'g:', t, Yfail(3, :), 'y-.', tc, Y_sim(3, :), 'b', 'LineWidth', 1.5);
%         plot(t, Y(3, :), 'g:', t, Yfail(3, :), 'y-.', tc, Y_sim(3, :), 'b', 'LineWidth', 1.5);
        hold off
        legend('', '', '', 'Actual-FTMPC', 'Measured-FTMPC', 'Location', 'NorthEast');
        legend boxoff
        print -dsvg figs/outputHE.svg
    end

	%% Manipulated variables
    figure(2)
    if FT == 2
        subplot(211)
        stairs(t, FTCS(FT).U(1, 1:end), 'b', 'LineWidth', 1.5)
        hold on
        xlabel('Time [min]'); ylabel('q_1 [l/min]'); grid on
        subplot(212)
        stairs(t, FTCS(FT).U(2, 1:end), 'b', 'LineWidth', 1.5)
        hold on
        xlabel('Time [min]'); ylabel('q_2 [l/min]'); grid on
	else
        subplot(211)
        stairs(t, FTCS(FT).U(1, 1:end), 'k-.', 'LineWidth', 1.5)
        plot(t, FTCS(FT).Umin(1)*ones(length(t)), 'r--')
        plot(t, FTCS(FT).Umax(1)*ones(length(t)), 'r--')
        hold off
        legend('MPC', 'FTMPC', 'Location', 'NorthWest');
        legend boxoff
        subplot(212)
        stairs(t, FTCS(FT).U(2, 1:end), 'k-.', 'LineWidth', 1.5)
        plot(t, FTCS(FT).Umin(2)*ones(length(t)), 'r--')
        plot(t, FTCS(FT).Umax(2)*ones(length(t)), 'r--')
        hold off
        print -dsvg figs/inputHE.svg
    end

%     % Failure inputs
%     figure(3)
%     if FTC == 0
%         subplot(211)
%         plot(t, FTCS(FTC).Umin(1, :), 'r--')
%         hold on
%         plot(t, FTCS(FTC).Umax(1, :), 'r--')
%         stairs(t, Ufail(1, :), 'b', 'LineWidth', 1.5)
%         xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
%         subplot(212)
%         plot(t, FTCS(FTC).Umin(2, :), 'r--')
%         hold on
%         plot(t, FTCS(FTC).Umax(2, :), 'r--')
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
    if FT == 2
        subplot(211)
        plot(t, FTCS(FT).RUIO(1).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        axis([0 inf 0 6.5])
        subplot(212)
        plot(t, FTCS(FT).RUIO(2).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        axis([0 inf 0 0.6])
	else
        subplot(211)
        plot(t, FTCS(FT).RUIO(1).error, 'k-.', 'LineWidth', 1.5)
        plot(t, threshold(1, :),  'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('|e_q|');
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
        legend boxoff
        subplot(212)
        plot(t, FTCS(FT).RUIO(2).error, 'k-.', 'LineWidth', 1.5)
        plot(t, threshold(2, :),  'r--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('|e_q|');
        hold off
        print -dsvg figs/RUIOerrorHE.svg
    end

	% UIOO error detection
    figure(5)
    if FT == 2
        subplot(211)
        plot(t, FTCS(FT).UIOO(1).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
        subplot(212)
        plot(t, UIOO(2).error, 'b', 'LineWidth', 1.5)
        hold on; grid on
    else
        subplot(211)
        plot(t, FTCS(FT).UIOO(1).error, 'k-.', 'LineWidth', 1.5)
        plot(t, threshold(3, :),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.04])
        xlabel('Time [min]'); ylabel('|e_x|');
        hold off
        legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthEast');
        legend boxoff
        subplot(212)
        plot(t, FTCS(FT).UIOO(2).error, 'k-.', 'LineWidth', 1.5)
        plot(t, threshold(4, :),  'r--', 'LineWidth', 1.5)
        axis([0 inf 0 0.03])
        xlabel('Time [min]'); ylabel('|e_x|');
        hold off
        print -dsvg figs/UIOOerrorHE.svg
    end

    % Fault estimation
    figure(6)
    if FT == 2
        subplot(211)
        stairs(t, FTCS(FT).RUIO(1).Fact, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, FTCS(FT).Ufail(1, :) - FTCS(FT).U(1, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_1 [l/min]'); grid on
        subplot(212)
        stairs(t, FTCS(FT).RUIO(2).Fact, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, FTCS(FT).Ufail(2, :) - FTCS(FT).U(2, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('Q_2 [l/min]'); grid on
    else
        subplot(211)
        stairs(t, FTCS(FT).RUIO(1).Fact, 'k-.', 'LineWidth', 1.5)
        hold off
        subplot(212)
        stairs(t, FTCS(FT).RUIO(2).Fact, 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/actuatorEstimationHE.svg
    end

    % Fault estimation
    figure(7)
    if FT == 2
        subplot(211)
        stairs(t, FTCS(FT).UIOO(1).Fsen, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, FTCS(FT).Yfail(1, :) - FTCS(FT).Y(1, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_1 [K]'); grid on
        subplot(212)
        stairs(t, FTCS(FT).UIOO(2).Fsen, 'b', 'LineWidth', 1.5)
        hold on
        stairs(t, FTCS(FT).Yfail(2, :) - FTCS(FT).Y(2, :), 'm--', 'LineWidth', 1.5)
        xlabel('Time [min]'); ylabel('\theta_2 [K]'); grid on
    else
        subplot(211)
        stairs(t, FTCS(FT).UIOO(1).Fsen, 'k-.', 'LineWidth', 1.5)
        hold off
        subplot(212)
        stairs(t, FTCS(FT).UIOO(2).Fsen, 'k-.', 'LineWidth', 1.5)
        hold off
        print -dsvg figs/sensorEstimationHE.svg
    end
    
%     % Objective
%     figure(8)
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
%     figure(9)
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
%     figure(10)
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
%     figure(11)
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
%     figure(12)
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
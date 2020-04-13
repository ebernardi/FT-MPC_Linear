%% Set Plots
clear all; close all; clc;

load FTCS

vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .7; .5];
tc = 0:Ts/10:Nsim/20;

disp('Plotting...')
FTC_OFF = 1; FTC_ON = 2;
%% Outputs
figure(1)

subplot(311)
plot(t, xsp(1)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
hold on; grid on
plot(t, FTCS(FTC_OFF).Y(1, :), 'b', t, FTCS(FTC_OFF).Yfail(1, :), 'y-.', 'LineWidth', 1.5);
plot(t, FTCS(FTC_ON).Y(1, :), 'k-.', t, FTCS(FTC_ON).Yfail(1, :), 'g:', 'LineWidth', 1.5);
axis([0 inf 98.5 102.5]); hold off
xlabel('Time [min]'); ylabel('Vol [l]');
leg = legend('x_s', 'Actual-MPC', 'Measured-MPC', 'Location', 'SouthWest');
set(leg, 'FontSize', 8, 'Orientation', 'horizontal');
leg.ItemTokenSize = [20, 15];

subplot(312)
plot(t, xsp(2)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
hold on; grid on
plot(t, FTCS(FTC_OFF).Y(2, :), 'b', t, FTCS(FTC_OFF).Yfail(2, :), 'y-.', 'LineWidth', 1.5);
h(1) = plot(t, FTCS(FTC_ON).Y(2, :), 'k-.', 'LineWidth', 1.5);
h(2) = plot(t, FTCS(FTC_ON).Yfail(2, :), 'g:', 'LineWidth', 1.5);
axis([0 inf 0.075 0.115]); hold off
xlabel('Time [min]'); ylabel('C_A [mol/l]');
leg = legend(h(:), 'Actual-FTMPC', 'Measured-FTMPC', 'Location', 'NorthWest');
set(leg, 'FontSize', 8, 'Orientation', 'horizontal');
leg.ItemTokenSize = [20, 15];

subplot(313)
plot(t, xsp(3)*ones(1, length(t)), 'r:', 'LineWidth', 1.5);
hold on; grid on
plot(t, FTCS(FTC_OFF).Y(3, :), 'b', t, FTCS(FTC_OFF).Yfail(3, :), 'y-.', 'LineWidth', 1.5);
plot(t, FTCS(FTC_ON).Y(3, :), 'k-.', t, FTCS(FTC_ON).Yfail(3, :), 'g:', 'LineWidth', 1.5);
axis([0 inf 436 444]);
hold off
xlabel('Time [min]'); ylabel('Temp [K]');

print -dsvg figs/outputCSTR.svg

%% Manipulated variables (with actuator faults)
figure(2)
subplot(211)
stairs(t, FTCS(FTC_OFF).Ufail(1, 1:end), 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Ufail(1, 1:end), 'k-.', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_OFF).Umin(1, 1:end), 'r--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_OFF).Umax(1, 1:end), 'r--', 'LineWidth', 1.5)
hold off
xlabel('Time [min]'); ylabel('q_s [l/min]');
legend('MPC', 'FTMPC', 'Location', 'NorthEast');
% axis([0 inf 93 100])

subplot(212)
stairs(t, FTCS(FTC_OFF).Ufail(2, 1:end), 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Ufail(2, 1:end), 'k-.', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_OFF).Umin(2, 1:end), 'r--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_OFF).Umax(2, 1:end), 'r--', 'LineWidth', 1.5)
% axis([0 inf 7.5 8.5])
hold off
xlabel('Time [min]'); ylabel('q_c [l/min]');

print -dsvg figs/inputCSTR.svg

%% RUIO error detection
figure(4)
subplot(211)
plot(t, FTCS(FTC_OFF).RUIO(1).error, 'b', 'LineWidth', 1.5)
hold on; grid on
plot(t, FTCS(FTC_ON).RUIO(1).error, 'k-.', 'LineWidth', 1.5)
plot(t, threshold(1, :),  'r--', 'LineWidth', 1.5)
hold off
axis([0 inf 0 0.8])
xlabel('Time [min]'); ylabel('|e_q|');

subplot(212)
plot(t, FTCS(FTC_OFF).RUIO(2).error, 'b', 'LineWidth', 1.5)
hold on; grid on
plot(t, FTCS(FTC_ON).RUIO(2).error, 'k-.', 'LineWidth', 1.5)
plot(t, threshold(2, :),  'r--', 'LineWidth', 1.5)
hold off
axis([0 inf 0 0.5])
xlabel('Time [min]'); ylabel('|e_q|');
legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');

print -dsvg figs/RUIOerrorCSTR.svg

%% UIOO error detection
figure(5)
subplot(211)
plot(t, FTCS(FTC_OFF).UIOO(1).error, 'b', 'LineWidth', 1.5)
hold on; grid on
plot(t, FTCS(FTC_ON).UIOO(1).error, 'k-.', 'LineWidth', 1.5)
plot(t, threshold(3, :),  'r--', 'LineWidth', 1.5)
hold off
axis([0 inf 0 1.5e-4])
xlabel('Time [min]'); ylabel('|e_x|');
% legend('MPC', 'FTMPC', 'Threshold', 'Location', 'NorthWest');

subplot(212)
plot(t, FTCS(FTC_OFF).UIOO(2).error, 'b', 'LineWidth', 1.5)
hold on; grid on
plot(t, FTCS(FTC_ON).UIOO(2).error, 'k-.', 'LineWidth', 1.5)
plot(t, threshold(4, :),  'r--', 'LineWidth', 1.5)
hold off
axis([0 inf 0 2e-1])
xlabel('Time [min]'); ylabel('|e_x|');

print -dsvg figs/UIOOerrorCSTR.svg

%% Actuator fault estimation
figure(6)
subplot(211)
stairs(t, FTCS(FTC_OFF).RUIO(1).Fact, 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Ufail(1, :) - FTCS(FTC_ON).U(1, :), 'm--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).RUIO(1).Fact, 'k-.', 'LineWidth', 1.5)
hold off
% axis([0 inf 0 6])
xlabel('Time [min]'); ylabel('Q_s [l/min]');
legend('MPC', 'Actual fault', 'FTMPC', 'Location', 'NorthEast');

subplot(212)
stairs(t, FTCS(FTC_OFF).RUIO(2).Fact, 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Ufail(2, :) - FTCS(FTC_ON).U(2, :), 'm--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).RUIO(2).Fact, 'k-.', 'LineWidth', 1.5)
hold off
% axis([0 inf -0.5 0])
xlabel('Time [min]'); ylabel('Q_c [l/min]');

print -dsvg figs/actuatorEstimationCSTR.svg

%% Sensor fault estimation
figure(7)
subplot(211)
stairs(t, FTCS(FTC_OFF).UIOO(1).Fsen, 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Yfail(1, :) - FTCS(FTC_ON).Y(1, :), 'm--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).UIOO(1).Fsen, 'k-.', 'LineWidth', 1.5)
hold off
axis([0 inf 0 2.5]);
xlabel('Time [min]'); ylabel('V [l]');
legend('MPC', 'Actual fault', 'FTMPC', 'Location', 'NorthWest');

subplot(212)
stairs(t, FTCS(FTC_OFF).UIOO(2).Fsen, 'b', 'LineWidth', 1.5)
hold on; grid on
stairs(t, FTCS(FTC_ON).Yfail(3, :) - FTCS(FTC_ON).Y(3, :), 'm--', 'LineWidth', 1.5)
stairs(t, FTCS(FTC_ON).UIOO(2).Fsen, 'k-.', 'LineWidth', 1.5)
hold off
axis([0 inf -3 0]);
xlabel('Time [min]'); ylabel('Temperature [K]');

print -dsvg figs/sensorEstimationCSTR.svg

% %%     % Objective
% %     figure(8)
% %     if FTC == 0
% %         plot(t, Obj, 'b', 'LineWidth', 1.5)
% %         hold on
% %         xlabel('Time [min]'); ylabel('Objective'); grid on
% %     else
% %         plot(t, Obj, 'k-.', 'LineWidth', 1.5)
% %         hold off
% %         axis([0 inf 0 400])
% %         print -dsvg figs/objectiveCSTR.svg
% %     end
% %     
% %% State evolution
% figure(9)
% %if FTC == 0
%     %plot3(x0(1), x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%     %hold on
%     %plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
%     %plot3(Y(1, :), Y(2, :), Y(3, :), 'b', 'LineWidth', 1.5)
%     %plot3(Y(1, end), Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
%     plot3(xsp(1), xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%     hold on
%     plot(Xpoly+X_lin, 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%     plot(Xs+X_lin, 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%     xlabel('V [l]'); ylabel('C_A [mol/l]'); zlabel('Temperature [K]'); grid on;
% %else
%    % plot3(Y_sim(1, :), Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)        
%     %plot3(Y(1, :), Y(2, :), Y(3, :), 'k-.', 'LineWidth', 1.5)
%     %plot3(Y(1, end), Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
%    % hold off
%     %print -dsvg figs/stateCSTR.svg
% %end
% %     
% %% State evolution
% Xx = Xpoly.projection(1:2).minHRep();
% Xxs = Xs.projection(1:2).minHRep();
% 
% figure(10)
% % if FTC == 0
% %     plot(x0(1), x0(2), 'g*', 'LineWidth', 1.5);
% %     hold on
% %     plot(Y_sim(1, :), Y_sim(2, :), 'y', 'LineWidth', 1.5)        
% %     plot(Y(1, :), Y(2, :), 'b.', 'LineWidth', 1.5)
% %     plot(Y(1, end), Y(2, end), 'mo', 'LineWidth', 1.5)
%     plot(xsp(1), xsp(2), 'rp', 'LineWidth', 1.5)
%     hold on
%     plot(Xxs+X_lin(1:2), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%     plot(Xx+X_lin(1:2), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%     xlabel('V [l]'); ylabel('C_A [mol/l]'); grid on;
% % else
% %     plot(Y_sim(1, :), Y_sim(2, :), 'g', 'LineWidth', 1.5)        
% %     plot(Y(1, :), Y(2, :), 'k.', 'LineWidth', 1.5)
% %     plot(Y(1, end), Y(2, end), 'ro', 'LineWidth', 1.5)
% %     hold off
% %     print -dsvg figs/state1-2CSTR.svg
% % end
% %     
% %% State evolution
% Xx = Xpoly.projection(1:2:3).minHRep();
% Xxs = Xs.projection(1:2:3).minHRep();
% 
% figure(11)
% % if FTC == 0
% %     hold on
% %     plot(Y_sim(1, :), Y_sim(3, :), 'b', 'LineWidth', 1.5)
% %     plot(Y(1, end), Y(3, end), 'mo', 'LineWidth', 1.5)
% % else
% %     plot(Y_sim(1, :), Y_sim(3, :), 'k--', 'LineWidth', 1.5)
% %     plot(Y(1, end), Y(3, end), 'y*', 'LineWidth', 1.5)
%     plot(xsp(1), xsp(3), 'ro', 'LineWidth', 1.5)
%     hold on
%     plot(x0(1), x0(3), 'gd', 'LineWidth', 1.5);
%     plot(Xxs+X_lin(1:2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'LineWidth', 1.5)
%     plot(Xx+X_lin(1:2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'LineWidth', 1.5)
%     xlabel('V [l]'); ylabel('Temperature [K]'); grid on;
% 
% %         hold off
% %         leg = legend('$x(0)$', '$x_{MPC}$', '$x_{MPC}(end)$', '$x_{FTMPC}$', '$x_{FTMPC}(end)$', '$x_s$', '$\bf{X}_s$', '$\bf{X}$', 'Location', 'SouthEast');
% %         set(leg, 'Interpreter', 'latex');
% %         print -dsvg figs/state1-3CSTR.svg
% %     end
% %     
% %% State evolution
% Xx = Xpoly.projection(2:3).minHRep();
% Xxs = Xs.projection(2:3).minHRep();
% 
% figure(12)
% % if FTC == 0
%     plot(x0(2), x0(3), 'g*', 'LineWidth', 1.5);
%     hold on
% %     plot(Y_sim(2, :), Y_sim(3, :), 'y', 'LineWidth', 1.5)
% %     plot(Y(2, :), Y(3, :), 'b.', 'LineWidth', 1.5)
% %     plot(Y(2, end), Y(3, end), 'mo', 'LineWidth', 1.5)
%     plot(xsp(2), xsp(3), 'rp', 'LineWidth', 1.5)
%     plot(Xxs+X_lin(2:3), 'Color', gris, 'Alpha', 0.2, 'edgecolor', gris, 'linestyle', '--', 'LineWidth', 1.5)
%     plot(Xx+X_lin(2:3), 'Color', vecrojo, 'Alpha', 0.05, 'edgecolor', vecrojo, 'linestyle', '--', 'LineWidth', 1.5)
%     xlabel('C_A [mol/l]'); ylabel('Temperature [K]'); grid on;
% %     else
% %         plot(Y_sim(2, :), Y_sim(3, :), 'g', 'LineWidth', 1.5)
% %         plot(Y(2, :), Y(3, :), 'k.', 'LineWidth', 1.5)
% %         plot(Y(2, end), Y(3, end), 'ro', 'LineWidth', 1.5)
% %         hold off
% %         print -dsvg figs/state2-3CSTR.svg    
% %     end
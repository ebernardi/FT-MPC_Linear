syms rho1 rho2 rhop Cp1 Cp2 Cpp Ar h1 h2 V1 V2 Vp theta1e theta2e q1 q2 theta1s theta2s thetap

% El fluido 2 es calefactor
% El fluido 1 es de proceso

%% Parámetros
Rho1 = 1;               % Densidad del fluido 1 (kg/l)
Rho2 = 1;               % Densidad del fluido 2 (kg/l)
Rhop = 7.874;        % Densidad de la pared (kg/l)
Cp_1 = 1000;         % Calor especifico del fluido 1 (cal/kg K)
Cp_2 = 1000;          % Calor especifico del fluido 2 (cal/kg K)
Cp_p = 1075.53;     % Calor especifico de la pared (cal/kg K)
a = 0.881;               % Area de intercambio (m^2)
h_1 = 32374;          % Coeficiente de transferencia de calor para fluido 1 (cal/min K m^2)
h_2 = 14716.6667; % Coeficiente de transferencia de calor para fluido 2 (cal/min K m^2)
V_1 = 16;                % Volumen de tubos (l)
V_2 = 2.11;             % Volumen de carcaza (l)
V_p = 1.19;             % Volumen de pared (l)
Theta_1e = 480;       % Temperatura de entrada de fluido 1 (K) (435 +/-10)
Theta_2e = 900;       % Temperatura de entrada de fluido 2 (K)

%% Modelo no lineal
system = [(q1*rho1*Cp1*(theta1e-theta1s)-Ar*h1*(theta1s-thetap))/(rho1*V1*Cp1);
                  (q2*rho2*Cp2*(theta2e-theta2s)+Ar*h2*(thetap-theta2s))/(rho2*V2*Cp2);
                  (Ar*h1*(theta1s-thetap)-Ar*h2*(thetap-theta2s))/(rhop*Cpp*Vp)];

inputs = [q1 q2];
outputs = [theta1s theta2s thetap];
states = [theta1s theta2s thetap];
nx = length(states); nu = length(inputs); ny = length(outputs);
C = eye(ny, nx);        % Matriz de salida
D = zeros(ny, nu);     % Matriz entrada/salida

%% Linealización
% Matrices simbólicas
A_sym = jacobian(system, states);
B_sym = jacobian(system, inputs);

% Obtengo la temperatura de la pared
Theta_p = (h_2*Theta_2s + h_1*Theta_1s) / (h_2 + h_1);

% Estados para la linealización
X_lin = [Theta_1s; Theta_2s; Theta_p];

% Obtengo el caudal de salida del fluido 1
Q1 = ( a*h_1*(Theta_1s-Theta_p) ) / ( Rho1*Cp_1*(Theta_1e - Theta_1s) );

% Obtengo el caudal de salida del fluido 2
Q2 = ( -a*h_2*(Theta_p - Theta_2s) ) / ( Rho2*Cp_2*(Theta_2e - Theta_2s) );

% Comprobación                  
% Theta_2s = ( (h_1+h_2)*Theta_p - h_1*Theta_1s ) / h_2
% Theta_2s = (Q2*Rho2*Cp_2*Theta_2e + (a*h_2*h_1*Theta_1s / (h_1+h_2)) )/ ...
%                       (Q2*Rho2*Cp_2 + a*h_2 - (a*h_2^2 / (h_1+h_2)))
% Theta_p = ( (Q1*Rho1*Cp_1 + a*h_1)*Theta_1s - Q1*Rho1*Cp_1*Theta_1e ) / (a*h_1)

% Manipuladas para la linealización
U_lin = [Q1; Q2];

% Matrices del sistema lineal
A = subs(A_sym, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
          {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, X_lin(1), X_lin(2), X_lin(3), U_lin(1), U_lin(2)});
B = subs(B_sym, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
          {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, X_lin(1), X_lin(2), X_lin(3), U_lin(1), U_lin(2)});
A = double(A);
B = double(B);

f = subs(system, {rho1, rho2, rhop, Cp1, Cp2, Cpp, Ar, h1, h2, V1, V2, Vp, theta1e, theta2e, theta1s, theta2s, thetap, q1, q2}, ...
          {Rho1, Rho2, Rhop, Cp_1, Cp_2, Cp_p, a, h_1, h_2, V_1, V_2, V_p, Theta_1e, Theta_2e, X_lin(1), X_lin(2), X_lin(3), U_lin(1), U_lin(2)});
f = double(f);

% Desviación del modelo lineal
delta = f - (A*X_lin+B*U_lin);
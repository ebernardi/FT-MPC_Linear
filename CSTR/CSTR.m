syms qs qc V CA T

%% Parameters
E_R = 1e4;          % [°K] (Activation energy term)
Te = 350;            % [°K] (Feed temperature)
Tce = 350;          % [°K] (Inlet coolant temperature)
dH = -2e5;          % [cal/mol] (Heat of reaction)
Cp = 1;               % [cal/g °K] (Specific heats)
rho = 1e3;			% [g/l] (Liquid densities)
CAe = 1;             % [mol/l] (Feed concentration)
ha = 7e5;			% [cal/min °K] (Heat transfer term)
k0 = 7.2e10;            % [l/min] (Reaction rate constant)
k1 = dH*k0/(rho*Cp);
k2 = rho*Cp/(rho*Cp);
k3 = ha/(rho*Cp);
k4 = 10;                % [l/min m^3/2] (Valve term)
q = 100;                % [l/min] (Feed flow rate)

%% Non-Linear model
system = [(q - qs);
          ((q/V)*(CAe - CA) - k0*CA*exp(-E_R/T));
		  ((q/V)*(Te - T) - k1*CA*exp(-E_R/T) + k2*(qc/V)*(1 - exp(-k3/qc))*(Tce - T))];

inputs = [qs qc];
outputs = [V CA T];
states = [V CA T];

nu = length(inputs); nx = length(states); ny = length(outputs);
C = eye(ny, nx);       % Output matrix
D = zeros(ny, nu);    % Input/Output matrix

%% Linealization
% Symbolic matrices
A_sym = jacobian(system, states);
B_sym = jacobian(system, inputs);

% Reactor temperature
% Tr = -E_R/log(-(q*(Ca-CAe))/(k0*Ca*Vr));

% Output Concentrarion
Ca = CAe/(1+(k0*(Vr/q)*exp(-E_R/Tr)));

% Linear states
X_lin = [Vr; Ca; Tr];

% Output flow rate
Qs = k4*sqrt(Vr);

% Coolant flow rate
Qc = double(solve((q/Vr)*(Te-Tr) - k1*Ca*exp(-E_R/Tr) + k2*(qc/Vr)*(1-exp(-k3/qc))*(Tce-Tr) == 0));

% Inputs
U_lin = [Qs; Qc];

% Linear systems matrices
A = subs(A_sym, {V, CA, T, qs, qc}, {Vr, Ca, Tr, Qs, Qc});
B = subs(B_sym, {V, CA, T, qs, qc}, {Vr, Ca, Tr, Qs, Qc});
f = subs(system, {V, CA, T, qs, qc}, {Vr, Ca, Tr, Qs, Qc});

A = double(A);
B = double(B);
f = double(f);

% Constant term
delta = f - (A*X_lin+B*U_lin);

% Continuous-time system
sys = ss(A, B, C, D);

% Euler discretization method
Ad = (A*Ts) + eye(nx); Bd = B*Ts; Cd = C; Dd = D; deltad = delta*Ts;
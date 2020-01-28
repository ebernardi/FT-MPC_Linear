clc; clear all; close all; yalmip('clear');

%% LTI model
Ts = 0.1;                       % Sample time [min]
run CSTR;

% Euler discretization method
Ad = (A*Ts) + eye(nx); Bd = B*Ts; Cd = C; Dd = D; deltad = delta*Ts;

%% Parameters
% Constraints
xmin = [95; 0.06; 350]-X_lin;
xmax = [105; 0.12; 500]-X_lin;
umin = [95; 95]-U_lin;
umax = [105; 105]-U_lin;

x0 = [102; 0.07; 445];             % Start-point
xsp = [100; 0.09; 440.7947];   % Set-point

%% Max reachable set
X = Polyhedron('lb', xmin, 'ub', xmax);     % State polyhedron
U = Polyhedron('lb', umin, 'ub', umax);    % Input polyhedron
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
	if Xs.contains(Xo)
        break;
	end
end

%% Set Plots
vecrojo = [0.7; 0; 0]; vecverde = [0; 0.8; 0]; vecazul = [0; 0; 0.6]; negro = [.1; .1; .1]; gris = [.5; .5; .5];
figure(1)
hold on
plot3(x0(1), x0(2), x0(3), 'ko');
plot3(xsp(1), xsp(2), xsp(3), 'bo');
plot(Xs+X_lin, 'Color', vecverde, 'Alpha', 0.1, 'edgecolor', vecverde, 'LineWidth', 1.5)
plot(X+X_lin, 'Alpha', 0.2, 'edgecolor', 'b', 'linestyle', '--')
hold off

figure(2)
plot(U+U_lin)
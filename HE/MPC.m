yalmip('clear')

%% Controller definition
x = sdpvar(nx*ones(1, Np+1), ones(1, Np+1));
u = sdpvar(nu*ones(1, Np), ones(1, Np));
xs = sdpvar(nx, 1);
uf = sdpvar(nu, 1);

% Artificial variables
xa = sdpvar(nx, 1); 
ua = sdpvar(nu, 1);

objective = 0; constraints = [];

% Stage constraints and objective
for k = 1:Np
    % Objective
    objective = objective + (x{k}-xa)'*Qx*(x{k}-xa); 
    objective = objective + (u{k}+uf-ua)'*Rx*(u{k}+uf-ua);

    % Dynamic constraint    
    constraints = [constraints, x{k+1} == Ad*x{k} + Bd*(u{k}+uf) + deltad];

    % Box-type constraint
    constraints = [constraints, umin <= u{k}+uf <= umax];
    constraints = [constraints, xmin <= x{k} <= xmax];
end

% Terminal constraints
objective = objective + (xa-xs)'*gamma*(xa-xs);                              % Terminal cost
constraints = [constraints, x{Np+1} == xa];                                    % Equality terminal constraint
constraints = [constraints, xa == Ad*xa + Bd*ua + deltad];            % Artificial variables equilibirum condition
% objective = objective + gamma*norm((xa-xs), 1);                         % Terminal cost
% constraints = [constraints, x{N+1} == xs];                                  % Strong constraints

% Defining the parameters in, and the solution
parameters = {x{1}, xs, uf};
solution = {u{1}, objective};

% Options for Optimizer  
% options = sdpsettings('solver', 'fmincon');
options = sdpsettings('solver', 'quadprog');     % quadprog works faster
options.verbose = 0;                                        % 0 to none, 2 to debug, 2+ for more

% Create prototype function
mpc = optimizer(constraints, objective, options, parameters, solution);
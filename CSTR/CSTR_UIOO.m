% Type of observer(N° observer). Sub-observer(N° sub-observer). Matrix 
UIOO = struct;

% Dimension of system matrices
n = 3; p = 2;

for j = 1:N
    UIOO(j).H = zeros(size(C));
    if j == 1
        UIOO(j).H(3, :) = C(3, :);
        UIOO(j).alpha = 0.3;
    else
        UIOO(j).H(1, :) = C(1, :);
        UIOO(j).alpha = 0.8;
    end
    UIOO(j).T2 = null(UIOO(j).H, 'r')';
    UIOO(j).J = UIOO(j).T2*C;

    yalmip('clear');
    X = sdpvar(n);
    S = sdpvar(n, p);

    const = [];
    const = [const, X >= 0];

    if j == 1
        UIOO(j).F = zeros(size(Bd(:, 1)));
    else
        UIOO(j).F = Bd(:, 1);
    end

    UIOO(j).W = sdpvar(n, p);

    UIOO(j).LMI = [2*UIOO(j).alpha*X, (Ad'*X+Ad'*UIOO(j).J'*S'-UIOO(j).J'*UIOO(j).W');
            (X*Ad+S*UIOO(j).J*Ad-UIOO(j).W*UIOO(j).J) 2*UIOO(j).alpha*X];

    const = [const, UIOO(j).LMI <= 0, (X+S*UIOO(j).J)*UIOO(j).F == 0];

    diagnostics = optimize(const);

    string = yalmiperror(diagnostics.problem);
    if diagnostics.problem == 0
        disp('Factible!')
    else
        disp('Infactible!')
        return
    end

    % Matrices to be determinated
    X = double(X);
    S = double(S);

    UIOO(j).E = X\S;
    UIOO(j).T1 = (eye(n) + UIOO(j).E*UIOO(j).J);

    UIOO(j).W = double(UIOO(j).W);

    UIOO(j).K = X\UIOO(j).W;
    UIOO(j).G = UIOO(j).T1*Bd;
    UIOO(j).Tg = UIOO(j).T1*deltad;
    UIOO(j).N = UIOO(j).T1*Ad - UIOO(j).K*UIOO(j).J;
    UIOO(j).L = UIOO(j).K - UIOO(j).N*UIOO(j).E;
end
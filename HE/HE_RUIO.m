% Type of observer(N° observer). Matrix 
RUIO = struct;

% Dimension of system matrices
n = 2; p = 3;

for j = 1:N
    
    if j == 1
        RUIO(j).N = [0 0; 0 1; 1 0];
        RUIO(j).Q = [0 0; 1 0; 0 1];
        RUIO(j).alpha = -0.2;
    else
        RUIO(j).N = [1 0; 0 0; 0 1];
        RUIO(j).Q = [0 1; 0 0; 1 0];
        RUIO(j).alpha = -0.2;
    end
%     if j == 1
%         RUIO(j).N = [0 0; 1 0; 0 1];
%         RUIO(j).Q = [0 0; 0 1; 1 0];
%         RUIO(j).alpha = -0.1;
%     else
%         RUIO(j).N = [0 0.1; 0 0; 1 0];
%         RUIO(j).Q = [1 0; 0 0; 0 1];
%         RUIO(j).alpha = -0.1;
%     end
    
    % Transformation matrix T for each model
    RUIO(j).D = Bd(:, j);
    RUIO(j).B = Bd;
    RUIO(j).T = [RUIO(j).N RUIO(j).D];
    RUIO(j).A_bar = RUIO(j).T\Ad*RUIO(j).T;
    RUIO(j).B_bar = RUIO(j).T\RUIO(j).B;
    RUIO(j).B_bar_1 = RUIO(j).B_bar(1:2, :);
    RUIO(j).B_bar_2 = RUIO(j).B_bar(3, :);
    RUIO(j).D_bar = RUIO(j).T\RUIO(j).D;
    RUIO(j).delta_bar = RUIO(j).T\deltad;
    RUIO(j).delta_bar_1 = RUIO(j).delta_bar(1:2, :);
    RUIO(j).delta_bar_2 = RUIO(j).delta_bar(3, :);    

    % Split state vector
    RUIO(j).A_bar_11 = RUIO(j).A_bar(1:2, 1:2);
    RUIO(j).A_bar_12 = RUIO(j).A_bar(1:2, 3);
    RUIO(j).A_bar_21 = RUIO(j).A_bar(3, 1:2);
    RUIO(j).A_bar_22 = RUIO(j).A_bar(3, 3);

    % Matrix U
    RUIO(j).rank_CD = rank(Cd*RUIO(j).D);
    RUIO(j).U = [Cd*RUIO(j).D RUIO(j).Q];
    RUIO(j).inv_U = inv(RUIO(j).U);
    RUIO(j).U_1 = RUIO(j).inv_U(1, :);
    RUIO(j).U_2 = RUIO(j).inv_U(2:3, :);

    % More matrices and check observability
    RUIO(j).A_tilde_1 = RUIO(j).A_bar_11 - RUIO(j).A_bar_12*RUIO(j).U_1*C*RUIO(j).N;
    RUIO(j).C_tilde_1 = Cd*RUIO(j).N;
    RUIO(j).E_1 = RUIO(j).A_bar_12*RUIO(j).U_1;

    RUIO(j).rank_C_tilde_1 = rank(RUIO(j).C_tilde_1);
    RUIO(j).O_M = [RUIO(j).C_tilde_1' (RUIO(j).C_tilde_1*RUIO(j).A_tilde_1)' (RUIO(j).C_tilde_1*RUIO(j).A_tilde_1*RUIO(j).A_tilde_1)']';
    RUIO(j).rank_Obs_M = rank(RUIO(j).O_M);
    
    % LMIs
    setlmis([]);
    X = lmivar(1, [n 1]);

    % LMI #1: X > 0
    lmiterm([-1 1 1 X], 1, 1);                  % LMI #1: X
    
	% LMI
    RUIO(j).W = lmivar(2, [n p]);

    % LMI #2: M < 0
    lmiterm([2 1 1 X], 2*RUIO(j).alpha, 1);  % LMI #1: 2*alpha*X; −left hand side
    lmiterm([2 2 1 RUIO(j).W], -1, RUIO(j).C_tilde_1, 's'); % LMI #1: -W*C_tilde_1; −left hand side
    lmiterm([2 2 1 X], 1, RUIO(j).A_tilde_1, 's'); % LMI #1: X*A_tilde_1; −left hand side
    lmiterm([2 2 2 X], 2*RUIO(j).alpha, 1);  % LMI #1: 2*alpha*X; −left hand side
    
    LMIs = getlmis;

    [~, xfeas] = feasp(LMIs);

    X = dec2mat(LMIs, xfeas, X);

    RUIO(j).W = dec2mat(LMIs, xfeas, RUIO(j).W);
    RUIO(j).L = X\RUIO(j).W;
    RUIO(j).K = RUIO(j).A_tilde_1 - RUIO(j).L*RUIO(j).C_tilde_1;
    RUIO(j).L_ast = RUIO(j).L + RUIO(j).E_1;

% Pole placement
%     % Compute gain matrices
%     if j == 1
%         RUIO(j).L = place(RUIO(j).A_tilde_1', RUIO(j).C_tilde_1', [-0.2 -0.2])';
%     else
%         RUIO(j).L = place(RUIO(j).A_tilde_1', RUIO(j).C_tilde_1', [-0.2 -0.2])';
%     end
%     RUIO(j).K = RUIO(j).A_tilde_1 - RUIO(j).L*RUIO(j).C_tilde_1;
%     RUIO(j).L_ast = RUIO(j).L + RUIO(j).E_1;
end
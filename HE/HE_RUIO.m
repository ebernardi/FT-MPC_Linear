% Type of observer(N° observer). Matrix 
RUIO = struct;

for j = 1:N
    
    if j == 1
        RUIO(j).N = [0 0; 0 1; 1 0];
        RUIO(j).Q = [0 0; 1 0; 0 1];
    else
        RUIO(j).N = [1 0; 0 0; 0 1];
        RUIO(j).Q = [0 1; 0 0; 1 0];
    end
    
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

    % Compute gain matrices
    if j == 1
        RUIO(j).L = place(RUIO(j).A_tilde_1', RUIO(j).C_tilde_1', [-0.2 -0.2])';
    else
        RUIO(j).L = place(RUIO(j).A_tilde_1', RUIO(j).C_tilde_1', [-0.2 -0.2])';
    end
	RUIO(j).K = RUIO(j).A_tilde_1 - RUIO(j).L*RUIO(j).C_tilde_1;
    RUIO(j).L_ast = RUIO(j).L + RUIO(j).E_1;
end


% %%  UIO 1 - Uso como entrada desconocida el caudal de ingreso al tanque Q
% D_1 = Bd(:, 1);
% 
% N_1 = [0 0; 0 1; 1 0];
% T_1 = [N_1 D_1];
% A_bar = T_1\Ad*T_1;
% B_bar = T_1\Bd;
% C_bar = Cd*T_1;
% D_bar = T_1\D_1;
% 
% % Separo el vector de estados no perturbados (n-q) del vector de
% % estados perturbados (q)
% A1_bar_11 = A_bar(1:2, 1:2);
% A1_bar_12 = A_bar(1:2, 3);
% A1_bar_21 = A_bar(3, 1:2);
% A1_bar_22 = A_bar(3, 3);
% 
% % Separo el vector de entradas no perturbadas (n-q) del vector de 
% % entradas perturbadas (q) 
% B1_bar_1 = B_bar(1:2, :);
% B1_bar_2 = B_bar(3, :);
% 
% % Armo la matriz no singular U
% Q1 = [0 0; 1 0; 0 1];
% U1 = [Cd*D_1 Q1];
% 
% % Obtengo inversa de U y la separo
% inv_U = inv(U1);
% U1_1 = inv_U(1, :);
% U2_1 = inv_U(2:3, :);
% 
% A_tilde_1 = A1_bar_11 - A1_bar_12*U1_1*Cd*N_1;
% E_1 = A1_bar_12*U1_1;
% C_tilde_1 = Cd*N_1;
% 
% % Verifico que el rango de C_tilde_1 sea (m-q)
% rank_C_tilde_1 = rank(C_tilde_1);
% 
% % Verifico Observabilidad
% O = [C_tilde_1' (C_tilde_1*A_tilde_1)' (C_tilde_1*A_tilde_1*A_tilde_1)']';
% rank_Obs = rank(O);
% 
% % Calculo la matriz de ganancias
% L = place(A_tilde_1', C_tilde_1', [-0.2 -0.2])';
% L_ast_1 = L + E_1;
% K_1 = A_tilde_1 - L*C_tilde_1;
% 
% %%  UIO 2 - Uso como entrada desconocida el caudal de ingreso a la camisa Qc
% D_2 = Bd(:, 2);
% 
% N_2 = [1 0; 0 0; 0 1];
% T_2 = [N_2 D_2];
% A_bar = T_2\Ad*T_2;
% B_bar = T_2\Bd;
% C_bar = Cd*T_2;
% D_bar = T_2\D_2;
% 
% % Separo el vector de estados no perturbados (n-q) del vector de
% % estados perturbados (q)
% A2_bar_11 = A_bar(1:2, 1:2);
% A2_bar_12 = A_bar(1:2, 3);
% A2_bar_21 = A_bar(3, 1:2);
% A2_bar_22 = A_bar(3, 3);
% 
% % Separo el vector de entradas no perturbadas (n-q) del vector de 
% % entradas perturbadas (q) 
% B2_bar_1 = B_bar(1:2, :);
% B2_bar_2 = B_bar(3, :);
% 
% % Armo la matriz no singular U
% Q2 = [0 1; 0 0; 1 0];
% U2 = [Cd*D_2 Q2];
% 
% % Obtengo la inversa de U y la separo
% inv_U = inv(U2);
% U1_2 = inv_U(1, :);
% U2_2 = inv_U(2:3, :);
% 
% A_tilde_1 = A2_bar_11-A2_bar_12*U1_2*Cd*N_2;
% E_1 = A2_bar_12*U1_2;
% C_tilde_1 = Cd*N_2;
% 
% % Verifico que el rango de C_tilde_1 sea (m-q)
% rank_C_tilde_1 = rank(C_tilde_1);
% 
% % Verifico Observabilidad
% O = [C_tilde_1' (C_tilde_1*A_tilde_1)' (C_tilde_1*A_tilde_1*A_tilde_1)']';
% rank_Obs = rank(O);
% 
% % Calculo la matriz de ganancias
% L = place(A_tilde_1', C_tilde_1', [-0.2 -0.2])';
% L_ast_2 = L + E_1;
% K_2 = A_tilde_1 - L*C_tilde_1;
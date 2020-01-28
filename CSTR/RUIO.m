% Verifico Observabilidad
O_cstr = [Cd' (Cd*Ad)' (Cd*Ad*Ad)']';
rank_Obs = rank(O_cstr);

%%  UIO 1 - Uso como entrada desconocida el caudal de ingreso al tanque Q
D_1 = Bd(:, 1);

% N_1 = [0 0; 0 -1; -1 0];
N_1 = [0 0; 0 1; 1 0];
T_1 = [N_1 D_1];
A_bar = T_1\Ad*T_1;
B_bar = T_1\Bd;
C_bar = Cd*T_1;
D_bar = T_1\D_1;

% Separo el vector de estados no perturbados (n-q) del vector de
% estados perturbados (q)
A1_bar_11 = A_bar(1:2, 1:2);
A1_bar_12 = A_bar(1:2, 3);
A1_bar_21 = A_bar(3, 1:2);
A1_bar_22 = A_bar(3, 3);

% Separo el vector de entradas no perturbadas (n-q) del vector de 
% entradas perturbadas (q) 
B1_bar_1 = B_bar(1:2, :);
B1_bar_2 = B_bar(3, :);

% Armo la matriz no singular U
Q1 = [0 0; 1 0; 0 1];
U1 = [Cd*D_1 Q1];

% Obtengo inversa de U y la separo
inv_U = inv(U1);
U1_1 = inv_U(1, :);
U2_1 = inv_U(2:3, :);

A_tilde_1 = A1_bar_11 - A1_bar_12*U1_1*Cd*N_1;
E_1 = A1_bar_12*U1_1;
C_tilde_1 = Cd*N_1;

% Verifico que el rango de C_tilde_1 sea (m-q)
rank_C_tilde_1 = rank(C_tilde_1);

% Verifico Observabilidad
O = [C_tilde_1' (C_tilde_1*A_tilde_1)' (C_tilde_1*A_tilde_1*A_tilde_1)']';
rank_Obs = rank(O);

% Calculo la matriz de ganancias
L = place(A_tilde_1', C_tilde_1', [-0.15 -0.2])';
L_ast_1 = L + E_1;
K_1 = A_tilde_1 - L*C_tilde_1;

%%  UIO 2 - Uso como entrada desconocida el caudal de ingreso a la camisa Qc
D_2 = Bd(:, 2);

% N_2 = [-0.9274 0; 0 -0.9274; 0 0]; %TODO: Ver como se propone%
N_2 = [1 0; 0 1; 0 0];
T_2 = [N_2 D_2];
A_bar = T_2\Ad*T_2;
B_bar = T_2\Bd;
C_bar = Cd*T_2;
D_bar = T_2\D_2;

% Separo el vector de estados no perturbados (n-q) del vector de
% estados perturbados (q)
A2_bar_11 = A_bar(1:2, 1:2);
A2_bar_12 = A_bar(1:2, 3);
A2_bar_21 = A_bar(3, 1:2);
A2_bar_22 = A_bar(3, 3);

% Separo el vector de entradas no perturbadas (n-q) del vector de 
% entradas perturbadas (q) 
B2_bar_1 = B_bar(1:2, :);
B2_bar_2 = B_bar(3, :);

% Armo la matriz no singular U
Q2 = [0 1; 1 0; 0 0];
U2 = [Cd*D_2 Q2];

% Obtengo la inversa de U y la separo
inv_U = inv(U2);
U1_2 = inv_U(1, :);
U2_2 = inv_U(2:3, :);

A_tilde_1 = A2_bar_11-A2_bar_12*U1_2*Cd*N_2;
E_1 = A2_bar_12*U1_2;
C_tilde_1 = Cd*N_2;

% Verifico que el rango de C_tilde_1 sea (m-q)
rank_C_tilde_1 = rank(C_tilde_1);

% Verifico Observabilidad
O = [C_tilde_1' (C_tilde_1*A_tilde_1)' (C_tilde_1*A_tilde_1*A_tilde_1)']';
rank_Obs = rank(O);

% Calculo la matriz de ganancias
L = place(A_tilde_1', C_tilde_1', [-0.2 -0.15])';
L_ast_2 = L + E_1;
K_2 = A_tilde_1 - L*C_tilde_1;
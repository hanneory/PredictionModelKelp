%% Simulation

% Forward Euler


h = 1;
N = numColsTime;
initialState = [0; 0; 0];
%x(1) = initialState;
y = [0; 0; 0];
% kelp(N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, 10, 10, 10, 10)
for n = 1:N
    T = Envdata.T(n);
    U = CurrentsDataset.Currentspeed(n);
    X = NutrientDataset.("Nutrient concentration")(n);
    I = IrradianceDataset.Irradiance(n); 
    y(n+1) = y(n) + kelp(N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, U, X, T, I)*h;
end
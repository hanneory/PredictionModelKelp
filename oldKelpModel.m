h = 0.0104;
Y_A = zeros(Nsample, NumberIterations);
Y_N = zeros(Nsample, NumberIterations);
Y_C = zeros(Nsample, NumberIterations);
dAdt = zeros(Nsample, NumberIterations);
dNdt = zeros(Nsample, NumberIterations);
dCdt = zeros(Nsample, NumberIterations);
C_content = zeros(Nsample, NumberIterations);
N_content = zeros(Nsample, NumberIterations);
netCarbonFixed = zeros(Nsample, NumberIterations);
grossFrond = zeros(Nsample, NumberIterations);
A_0 = A0_mu/100;
N_0 = 0.01;
C_0 = 0.05;
Y_A(:,1) = A_0;
% Y_A(:,1) = 37;
Y_N(:,1) = N_0;
Y_C(:,1) = C_0;
netCarbonFixed(:,1) = C_0;
grossFrond(:,1) = A0_mu/100;

for n = 1:(NumberIterations-1)
    A = Y_A(:,n);
    N = Y_N(:,n);
    C = Y_C(:,n);
    T = X_T(:, n);
    U = X_U(:,n);
    NO3 = X_XNO3(:,n);
    I = X_I(:,n);
    [A_dot, N_dot, C_dot, W_s, my, C_total, N_total, C_frac, N_frac] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, NO3, T, I, N, C, A);
    dAdt(:,n) = A_dot;
    dNdt(:,n) = N_dot;
    dCdt(:,n) = C_dot;
    netCarbonFixed(:, n+1) = netCarbonFixed(:, n) + netCarbon(W_s, C_dot).*h; 
    grossFrond(:, n+1) = grossFrond(:, n) + totalGross(my, A).*h; 
    C_content(:, n) = C_frac;
    N_content(:, n) = N_frac;
    Y_A(:, n+1) = Y_A(:, n) + A_dot.*h;
    Y_N(:, n+1) = Y_N(:, n) + N_dot*h;
    Y_C(:, n+1) = Y_C(:, n) + C_dot*h;
end

C_content(:, end) = C_content(:, end-1);
N_content(:, end) = N_content(:, end-1);

%% Plots

t = tiledlayout(2,2);

figure(1)
plot(time, Y_A);
grid on;
title('Area');

figure(2);
plot(time, Y_N);
grid on;
title('Nitrogen');

figure(3);
plot(time, Y_C);
grid on;
title('Carbon');

figure(4);
plot(time, C_content);
grid on;
title('Carbon content (fraction of dry weight');

figure(5);
plot(time, N_content);
grid on;
title('Nitrogen content (fraction of dry weight');

figure(6);
plot(time, grossFrond);
grid on;
title('Gross Frond area');

figure(7);
plot(time, netCarbonFixed);
grid on;
title('Gross Carbon');



%% Net Carbon Fixed
function C = netCarbon (W_s, dCdt)
C = W_s .* dCdt;
end

%% Total gross frond area
function A_tot = totalGross (my, A)
A_tot = my.*A;
end

%% Effect of temperature
function f_temp = EffectTemp (T_k)
T = T_k;
f_temp = zeros(size(T,1), 1);
for n = 1:size(T,1)
if (T(n)>=-1.8) && (T(n) < 10)
    f_temp(n) = 0.08*T(n) + 0.2;
elseif (T(n)>= 10) && (T(n)<= 15)
    f_temp(n) = 1;
elseif (T(n)>15) && (T(n)<=19)
    f_temp(n) = 19/4 - T(n)/4;
elseif (T(n)>19)
    f_temp(n) = 0;
end
end
end

%% Gross photosynthesis
function P = grossPhotosynthesis ( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I)
syms B;
T_PL = 271;
T_PH = 296;
P_max = (P_1.*exp(T_AP/T_R1-T_AP./T)) ./ (1+exp(T_APL./T - T_APL/T_PL) + exp(T_APH/T_PH - T_APH./T));
B_vec = zeros(size(T,1), 1);
B_0 = 1*10^(-9);
for i = 1:size(T,1)
    for n = 1:10
        fun = @(B) - (alpha*I_sat) / (log(1+alpha/B)) * (alpha/(alpha+B)) * ((B/(alpha+B))^(B/alpha)) - P_max(i);
        B_new = fminsearch(fun, B_0);
        B_0 = B_new;
    end
    B_vec(i) = B_new*100;
end
P_S = (alpha*I_sat)./(log(1+alpha./B_vec));
P = P_S .* (1-exp(-(alpha.*I)./P_S)) .* exp(- (B_vec.*I)./P_S);
end

%% Specific Growth Rate
function my = specificGrowthRate (f_area, f_photo, f_temp, N_min, C_min, N, C)
my = zeros(size(N,1), 1);
for n = 1:size(N,1)
if ((1-N_min/N(n)) <= (1-C_min/C(n)))
    my(n) = f_area(n) * f_photo * f_temp(n) * (1-N_min/N(n));
else
    my(n) = f_area(n) * f_photo * f_temp(n) * (1-C_min/C(n));
end
end
end

%% Model
function [A_dot, N_dot, C_dot, W_s, my, C_total, N_total, C_frac, N_frac] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, NO3, T, I, N, C, A)
% Effect of size on growth rate
f_area = m_1.*exp(-(A./A_O).^2) + m_2;

%Effect of temperature on growth rate
f_temp = EffectTemp(T);

% Seasonal influence on growth rate
%f_photo = a_1*(1+sign(lambda)*abs(lambda)^0.5) + a_2;
f_photo = 1.5;

% Frond erosion
ny = (10.^(-6).*exp(epsilon.*A)) ./ (1+10.^(-6).*(exp(epsilon.*A)-1));

% Nitrate uptake rate
J = J_max.* (NO3./(K_X + NO3)) .* ((N_max-N)./(N_max-N_min)) .* (1-exp(-U/U_065));

% Gross photosynthesis
P = grossPhotosynthesis(alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T+274.15, I);

% Temperature dependent respiration
R = R_1 .* exp((T_AR/T_R1)- T_AR./(T+274.15));

% Carbon exudation
E = 1 - exp(gamma.*(C_min-C));

% Specific growth rate
my = specificGrowthRate(f_area, f_photo, f_temp, N_min, C_min, N, C);

% Amount of frond area lost
A_lost = A.*(C_min - C) ./ C_struct;

% Structural weight
W_s = k_A.*A;

% Total dry weight
W_d = k_A .* (1+ k_N.*(N-N_min) + N_min + k_C.*(C-C_min) + C_min).*A;

% Total wet weight
W_w = k_A .* (k_dw^(-1) + k_N.*(N-N_min) + N_min + k_C.*(C-C_min) + C_min).*A;

% Total carbon content
C_total = (C + C_struct) .* W_s;

% Total nitrogen content
N_total = (N + N_struct) .* W_s;

% Carbon content (fraction of dry weight)
C_frac = C_total ./ W_d;

% Nitrogen content (fraction of dry weight)
N_frac = N_total ./W_d;

% Rate of change of frond area
A_dot = (my-ny).*A;

% Rate of change in nitrogen reserves
N_dot = k_A^(-1)*J-my.*(N+N_struct);

% Rate of change in carbon reserves
C_dot = k_A^(-1).*(P.*(1-E)-R) - (C+C_struct).*my;
% state_dot = [(my-ny)*A;
%     k_A^(-1)*J-my*(N+N_struct);
%     k_A^(-1)*(P*(1-E-R)) - (C+C_struct)*my];

end

%% KelpModel with values from dataset withour offset
% for n = 1:(NumberIterations-1)
%     A = Y(1,n);
%     N = Y(2,n);
%     C = Y(3,n);
%     [state_dot, A_lost, W_s, W_d, W_w, C_total, N_total, my] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, X(2, n), X(3, n), X(1, n), X(4, n), N, C, A);
% %     netCarbonFixed(n+1) = netCarbon(W_s, state_dot(3), n, n+1) + netCarbonFixed(n);
% %     grossFrond(n+1) = totalGross(my, A, n, n+1) + grossFrond(n);
%     Y(:, n+1) = Y(:,n) + state_dot*h;
% end5*10^(-5))/(1.0000 * eps^-9))
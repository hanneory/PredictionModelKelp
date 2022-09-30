%Parameters
A_O = 6;
alpha = 3.75*10^(-5);
C_min = 0.01;
C_struct = 0.20;
gamma = 0.5; 
epsilon = 0.22;
I_sat = 200;
J_max = 1.4*10^(-4);
k_A = 0.6;
k_dw = 0.0785;
k_C = 2.1213;
k_N = 2.72;
m_1 = 0.1085;
m_2 = 0.03;
my_max = 0.18;
N_min = 0.01;
N_max = 0.022;
N_struct = 0.01;
P_1 = 1.22*10^(-3);
P_2 = 1.44*10^(-3);
a_1 = 0.85;
a_2 = 0.3;
R_1 = 2.785*10^(-4);
R_2 = 5.429*10^(-4);
T_R1 = 285;
T_R2 = 290;
T_AP = 1694.4;
T_APH = 25924;
T_APL = 27774;
T_AR = 11033;
U_065 = 0.03;
K_X = 4;


% state = [Area;
%          Nitrogen;
%          Carbon];

%% Forward Euler
h = 1;
NumberIterations = 17520;
Y = zeros(3, NumberIterations);
netCarbonFixed = zeros(1, NumberIterations-1);
grossFrond = zeros(1, NumberIterations-1);
A_0 = 0.01;
C_0 = 0.01;
N_0 = 0.011;
Y(:,1) = [A_0;N_0;C_0];
for n = 1:NumberIterations-1
    T = TemperatureDataset.Temperature(n);
    U = CurrentsDataset.Currentspeed(n);
    X = NutrientDataset.("Nutrient concentration")(n);
    I = IrradianceDataset.Irradiance(n); 
    A = Y(1,n);
    N = Y(2,n);
    C = Y(3,n);
    [state_dot, A_lost, W_s, W_d, W_w, C_total, N_total, my] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, X, T, I, N, C, A);
    %netCarbonFixed(n) = netCarbon(W_s, state_dot(3), n, n+1);
    %grossFrond(n) = totalGross(my, A, n, n+1);
    Y(:, n+1) = Y(:,n) + state_dot*h;
end
X = 1:1:17520;
figure(1)
plot(X, Y(1,:));

figure(2)
plot(X, Y(2,:));

figure(3)
plot(X, Y(3,:));
% 
% figure (4)
% plot(17519, netCarbonFixed);
% 
% figure (5)
% plot(X-1, grossFrond);

%% Net Carbon Fixed
function C = netCarbon (W_s, dCdt, tmin, tmax)
    syms t
    fun = W_s * dCdt;
    C = int(fun, t, tmin, tmax);
end

%% Total gross frond area
function A_tot = totalGross (my, A, tmin, tmax)
    syms t
    fun = my*A;
    A_tot = int(fun, t, tmin, tmax);
end

%% Effect of temperature
function f_temp = EffectTemp (T_k)
    T = T_k - 272.15;
    if (T>=-1.8) && (T < 10)
        f_temp = 0.08*T + 0.2;
    elseif (T>= 10) && (T<= 15)
        f_temp = 1;
    elseif (T>15) && (T<=19)
        f_temp = 19/4 - T/4;
    elseif (T>19)
        f_temp = 0;
    end
end

%% Gross photosynthesis
function P = grossPhotosynthesis ( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I)
    syms B;
    T_PL = 271;
    T_PH = 296;
    P_max = (P_1*exp(T_AP/T_R1-T_AP/T)) / (1+exp(T_APL/T- T_APL/T_PL) + exp(T_APH/T_PH - T_PH/T));

    fun = @(B) - (alpha*I_sat) / (log(1+alpha/B)) * (alpha/(alpha+B)) * ((B/(alpha+B))^(B/alpha)) - P_max;
    B_0 = 10^(-9);
    B = fminsearch(fun, B_0);
    
    P_S = (alpha*I_sat)/(log(1+alpha/B));
    P = P_S * (1-exp(-(alpha*I)/P_S)) * exp(- (B*I)/P_S);
end

%% Specific Growth Rate
function my = specificGrowthRate (f_area, f_photo, f_temp, N_min, C_min, N, C)
   if ((1-N_min/N) <= (1-C_min/C))
       my = f_area * f_photo * f_temp * (1-N_min/N);
   else
       my = f_area * f_photo * f_temp * (1-C_min/C);
   end 
end

%% Model
function [state_dot, A_lost, W_s, W_d, W_w, C_total, N_total, my] = kelp(N_struct, C_struct, k_A, k_N, k_C, k_dw, N_min, C_min, m_1, m_2, A_O, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, alpha, I_sat, P_1, T_AP, T_APL,T_APH, U, X, T, I, N, C, A)

    f_area = m_1*exp(-(A/A_O)^2) + m_2;

    f_temp = EffectTemp(T);

    %f_photo = a_1*(1+sign(lambda)*abs(lambda)^0.5) + a_2;
    f_photo = 1;

    ny = (10^(-6)*exp(epsilon*A)) / (1+10^(-6)*(exp(epsilon*A)-1));
    
    J = J_max* (X/(K_X + X)) * ((N_max-N)/(N_max-N_min)) * (1-exp(-U/U_065));

    P = grossPhotosynthesis(alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I);

    R = R_1 * exp((T_AR/T_R1)-(T_AR/T));

    E = 1 - exp(gamma*(C_min-C));

    my = specificGrowthRate(f_area, f_photo, f_temp, N_min, C_min, N, C);
    
    A_lost = A*(C_min - C) / C_struct;

    W_s = k_A*A;

    W_d = k_A * (1+ k_N*(N-N_min) + N_min + k_C*(C-C_min) + C_min)*A;

    W_w = k_A * (k_dw^(-1) + k_N*(N-N_min) + N_min + k_C*(C-C_min) + C_min)*A;

    C_total = (C + C_struct) * W_s;

    N_total = (N + N_struct) * W_s;

    state_dot = [(my-ny)*A;
        k_A^(-1)*J-my*(N+N_struct);
        k_A^(-1)*(P*(1-E-R)) - (C+C_struct)*my];
end






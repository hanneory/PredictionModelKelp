%Parameters
A_0 = 6;
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
Tspan = 5:1:20;
Nspan = 0.01:0.002:0.022;
Cspan = 0.01:0.01:0.3;
%linspace

%%Moode variables
% my, B, P_S

%% State vairables
% state = [A; N; C]

%% Environmental variables
% T, I, U, X

%% Biomass variables
% W_s, W_w, W_d

% Initial state
% state = [Area;
%          Nitrogen;
%          Carbon];


%kelp (N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, 0.5, 0.5, 300, gamma);
%[t, y] = ode45(@(t,y) kelp (N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, a_1, a_2, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, 0.5, 0.5, 10, gamma), [0 20], state);
%grossPhotosynthesis( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, 278, 200)
%ForwardEuler (N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma);

%% Forward Euler
h = 1;
NumberIterations = 17520;
Y = zeros(3, NumberIterations+1);
A = 0.01;
C = 0.02;
N = 0.02;
Y(:,1) = [A;N;C];
for n = 1:NumberIterations
    T = TemperatureDataset.Temperature(n);
    U = CurrentsDataset.Currentspeed(n);
    X = NutrientDataset.("Nutrient concentration")(n);
    I = IrradianceDataset.Irradiance(n); 
    Y(:, n+1) = Y(:,n) + kelp(N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, U, X, T, I, N, C, A)*h;
end

%%
function f_temp = EffectTemp (T_k)
    T = T_k - 273.15;
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

function P = grossPhotosynthesis ( alpha, I_sat, P_1, T_AP, T_R1, T_APL,T_APH, T, I)
    syms B;
    T_PL = 271;
    T_PH = 296;
    P_max = (P_1*exp(T_AP/T_R1-T_AP/T)) / (1+exp(T_AP/T- T_APL/T_PL) + exp(T_APH/T_PH - T_PH/T));


    P(B) = (alpha*I_sat) / (log(1+alpha/B)) * (alpha/(alpha+B)) * ((B/(alpha+B))^(B/alpha));
    pretty(P(B));
    Pd(B) = diff(P(B));
    pretty(Pd(B))
    B0 = 1*10^-9;
    nmax = 10;
    B = B0;

    for i = 1:nmax
        B_next = B - P(B)/Pd(B);
        B = B_next;
    end
    
    P_max2 = double(P(B));
    P_S = (alpha*I_sat)/(log(1+alpha/B));
    P = P_S * (1-exp(-(alpha*I)/P_S)) * exp(- (B*I)/P_S);

end

function my = specificGrowthRate (f_area, f_photo, f_temp, N_min, C_min, N, C)
   if ((1-N_min/N) <= (1-C_min/C))
       my = f_area * f_photo * f_temp * (1-N_min/N);
   else
       my = f_area * f_photo * f_temp * (1-C_min/C);
   end 
end

function [state_dot] = kelp(N_struct, C_struct, k_A, N_min, C_min, m_1, m_2, A_0, epsilon, K_X, N_max, J_max, U_065, R_1, T_AR, T_R1, gamma, U, X, T, I, N, C, A)

    f_area = m_1*exp(-(A/A_0)^2) + m_2;

    f_temp = EffectTemp(T);

    %lambda = 1;

    %f_photo = a_1*(1+sign(lambda)*abs(lambda)^0.5) + a_2;
    f_photo = 1;

    ny = (10^(-6)*exp(epsilon*A))/(1+10^(-6)*(exp(epsilon*A)-1));
    
    J = J_max*(X/(K_X + X))*((N_max-N)/(N_max-N_min))*(1-exp(-U/U_065));

    P = 0.5;

    R = R_1 * exp((T_AR/T_R1)-(T_AR/T));

    E = 1 - exp(gamma*(C_min-C));


    %my = f_area * f_photo * f_temp * min([1-N_min/N, 1-C_min/C]);
    %my = f_area * f_photo * f_temp * (1-N_min/N);
    my = specificGrowthRate(f_area, f_photo, f_temp, N_min, C_min, N, C);

    state_dot = [(my-ny)*A;
        k_A^(-1)*J-my*(N+N_struct);
        k_A^(-1)*(P*(1-E-R)) - (C+C_struct)*my];
end





